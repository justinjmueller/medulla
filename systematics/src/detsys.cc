/**
 * @file detsys.cc
 * @brief Implementation of the DetsysCalculator class.
 * @details This file contains the implementation of the DetsysCalculator class.
 * The DetsysCalculator class is used to calculate the weights for the detector
 * systematics using a spline interpolation of the ratio of the nominal and
 * sample histograms. Configuration of the class is done using a TOML-based
 * configuration file.
 * @see cfg::ConfigurationTable 
 * @author mueller@fnal.gov
 */
#include <map>
#include <vector>
#include <random>
#include <algorithm>

#include "detsys.h"
#include "configuration.h"
#include "utilities.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TSpline.h"
#include "TFile.h"

std::string sys::detsys::DetsysCalculator::make_hist_key(const std::string & variable_name, const std::string & variation_name) const
{
    if(variables.size() == 1)
        return variable_name;
    return variation_name + "__" + variable_name;
}

std::string sys::detsys::DetsysCalculator::make_spline_key(const std::string & detsys_name, const std::string & variable_name) const
{
    if(variables.size() == 1)
        return detsys_name;
    return detsys_name + "__" + variable_name;
}

// Constructor for the DetsysCalculator class that initializes the class using
// the configuration table, the output file, and the input file. 
sys::detsys::DetsysCalculator::DetsysCalculator(cfg::ConfigurationTable & table, TFile * output, TFile * input)
{
    // Roll random z-scores to create a set of universes for later.
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0, 1);
    nuniverses = table.get_int_field("variations.nuniverses");
    for(size_t n(0); n < nuniverses; ++n)
        random_zscores.push_back(dist(gen));

    // Create the output directories to store the histograms and splines.
    // Also, load a few configuration details.
    histogram_directory = create_directory(output, table.get_string_field("output.histogram_destination"));
    result_directory = create_directory(output, table.get_string_field("variations.result_destination"));
    variations = table.get_string_vector("variations.keys");

    if(table.has_field("variations.variables"))
        variables = table.get_string_vector("variations.variables");
    else if(table.has_field("variations.variable"))
    {
        try
        {
            variables = table.get_string_vector("variations.variable");
        }
        catch(const cfg::ConfigurationError &)
        {
            variables = {table.get_string_field("variations.variable")};
        }
    }
    else
        throw cfg::ConfigurationError("The field variations.variable or variations.variables must be present in the configuration file.");

    if(variables.empty())
        throw cfg::ConfigurationError("The field variations.variables must contain at least one variable.");

    variable = variables.front();

    std::vector<std::vector<double>> bins_by_variable;
    if(table.has_field("variations.bins_list"))
    {
        bins_by_variable = table.get_double_matrix("variations.bins_list");
        if(bins_by_variable.size() != variables.size())
            throw cfg::ConfigurationError("variations.bins_list size must match variations.variables size.");
    }
    else if(table.has_field("variations.bins"))
    {
        try
        {
            bins_by_variable = table.get_double_matrix("variations.bins");
        }
        catch(const cfg::ConfigurationError &)
        {
            std::vector<double> shared_bins = table.get_double_vector("variations.bins");
            bins_by_variable.assign(variables.size(), shared_bins);
        }

        if(bins_by_variable.size() != variables.size())
            throw cfg::ConfigurationError("variations.bins size must match variations.variable size when using multi-variable mode.");
    }
    else
        throw cfg::ConfigurationError("The field variations.bins or variations.bins_list must be present in the configuration file.");

    bool variable_length = table.get_bool_field("variations.variable_length", false);

    // Loop over the variations and create the histograms. A "variation" is a
    // single sample that implements some change w.r.t. the nominal sample in
    // the detector model.
    for(std::string variation : variations)
    {
        double pot(0);
        // Check if the variation has an exposure tree instead of a histogram.
        std::string exp_tree_name = table.get_string_field("variations.origin") + variation + "/" + table.get_string_field("variations.tree") + "_exposure";
        if(input->Get(exp_tree_name.c_str()))
        {
            // If the exposure tree exists, use it to calculate the POT.
            // The tree has a branch "pot" that we need to sum over.
            TTree * exp_tree = input->Get<TTree>(exp_tree_name.c_str());
            double pot_value;
            exp_tree->SetBranchAddress("pot", &pot_value);
            for(int i(0); i < exp_tree->GetEntries(); ++i)
            {
                exp_tree->GetEntry(i);
                pot += pot_value; // Sum the POT values
            }
            pot /= 1e18; // Convert to 1e18 POT
        }
        else
        {
            // If the exposure tree does not exist, use the POT histogram.
            std::string pot_name = table.get_string_field("variations.origin") + variation + '/' + "POT";
            TH1D * h = (TH1D *) input->Get(pot_name.c_str());
            pot = h->GetBinContent(1) / 1e18; // Convert to 1e18 POT
        }
        std::cout << "Variation " << variation << " has " << pot << "e18 POT." << std::endl;

        std::string name = table.get_string_field("variations.origin") + variation + '/' + table.get_string_field("variations.tree");
        TTree * t = input->Get<TTree>(name.c_str());
        for(size_t ivar(0); ivar < variables.size(); ++ivar)
        {
            const std::string & configured_variable = variables[ivar];
            std::string hist_key = make_hist_key(configured_variable, variation);
            const std::vector<double> & bins = bins_by_variable[ivar];

            double value;
            t->SetBranchAddress(configured_variable.c_str(), &value);
            if(!variable_length && bins.size() == 3)
                histograms[hist_key] = new TH1D(hist_key.c_str(), configured_variable.c_str(), bins[0], bins[1], bins[2]);
            else
                histograms[hist_key] = new TH1D(hist_key.c_str(), configured_variable.c_str(), bins.size() - 1, bins.data());

            for(int i(0); i < t->GetEntries(); ++i)
            {
                t->GetEntry(i);
                histograms[hist_key]->Fill(value, 1.0 / pot);
            }
        }
    }

    // Loop over the detector systematics and create the splines. A single
    // entry in the "sys" block of the configuration file with type "variation"
    // specifies a single detector systematic parameter, but in general may
    // consist of multiple variations spanning the range of the parameter.
    for(cfg::ConfigurationTable & t : table.get_subtables("sys"))
    {
        // Skip non-variation detector systematics.
        if(t.get_string_field("type") != "variation")
            continue;
        
        // The "points" field specifies the variations to be used in the spline
        // construction. The "scale" field specifies the scale factors for each
        // variation, which is really just a way to weight the variations in the
        // spline construction. 
        std::vector<std::string> points = t.get_string_vector("points");
        std::vector<double> scale_factors = t.get_double_vector("scale");

        // Create a dummy histogram to store useful metadata for the spline
        // construction. Finally, the "zscores" field specifies the 
        // corresponding z-scores for each variation (how many standard
        // deviations the systematic parameter is from the nominal value).
        std::string name = t.get_string_field("name");
        zscores.insert(std::make_pair(name, t.get_double_vector("nsigma")));
        for(const std::string & configured_variable : variables)
        {
            std::string spline_key = make_spline_key(name, configured_variable);
            TH1D * base = histograms[make_hist_key(configured_variable, points[0])];
            int nbins = base->GetXaxis()->GetNbins();
            const double * xedges = base->GetXaxis()->GetXbins()->GetArray();
            hdummies.insert(std::make_pair(spline_key, new TH1D("hdummy", "hdummy", nbins, xedges)));

            // This block creates a TH2D that will be used to store the input for
            // the spline construction. The TH2D is filled with the ratio of the
            // variations to the nominal sample (across the range of the variable)
            // for each of the spline points. Each spline point is adjusted by the
            // scale factor configured in the "detsys" block.
            double ylow = *std::min_element(zscores[name].begin(), zscores[name].end());
            double yup = *std::max_element(zscores[name].begin(), zscores[name].end());
            TH2D * h = new TH2D("tmp", "tmp", nbins, xedges, points.size(), ylow, yup);
            for(size_t i(0); i < points.size(); ++i)
            {
                std::cout << "Creating spline input for systematic " << name << " variation " << points[i] << " with scale factor " << scale_factors[i] << " and z-score " << zscores[name][i] << std::endl;
                std::cout << "Variation histogram: " << make_hist_key(configured_variable, points[i]) << ", entries: " << histograms[make_hist_key(configured_variable, points[i])]->GetEntries() << std::endl;
                std::cout << "Nominal histogram: " << make_hist_key(configured_variable, t.get_string_field("ordinate")) << ", entries: " << histograms[make_hist_key(configured_variable, t.get_string_field("ordinate"))]->GetEntries() << std::endl;

                hdummies[spline_key]->Divide(histograms[make_hist_key(configured_variable, points[i])], histograms[make_hist_key(configured_variable, t.get_string_field("ordinate"))]);
                for(int j(0); j < hdummies[spline_key]->GetNbinsX(); ++j)
                {
                    if(hdummies[spline_key]->GetBinContent(j + 1) != 0)
                        h->SetBinContent(j + 1, i + 1,  1 + scale_factors[i] * (hdummies[spline_key]->GetBinContent(j + 1) - 1));
                    else
                        h->SetBinContent(j + 1, i + 1, 1);
                }
            }

            // This block creates the splines for the detector systematic parameter.
            // The splines are created for each bin of the dummy histogram.
            splines.insert(std::make_pair(spline_key, std::vector<TSpline3 *>()));
            for(int j(0); j < hdummies[spline_key]->GetNbinsX(); ++j)
            {
                std::vector<double> x, y;
                for(size_t i(0); i < points.size(); ++i)
                {
                    x.push_back(zscores[name][i]);
                    y.push_back(h->GetBinContent(j + 1, i + 1));
                }
                splines[spline_key].push_back(new TSpline3("spline", x.data(), y.data(), x.size()));
            }
        }
    }
}

// Add a variable to the list of result histograms.
void sys::detsys::DetsysCalculator::add_variable(SysVariable & variable)
{
    // Create the TH1D and TH2D that will be used to store the results of
    // the detector systematic universes.
    for(auto & [key, value] : zscores)
    {
        if(variables.size() == 1)
        {
            std::string name = variable.name + "_" + key;
            detsys_results1D[name] = new TH1D(name.c_str(), name.c_str(), 1000, -0.25, 0.25);
            detsys_results2D[name] = new TH2D(name.c_str(), name.c_str(), variable.nbins, variable.min, variable.max, nuniverses, 0, nuniverses);
        }
        else
        {
            for(const std::string & configured_variable : variables)
            {
                std::string name = variable.name + "_" + key + "_" + configured_variable;
                detsys_results1D[name] = new TH1D(name.c_str(), name.c_str(), 1000, -0.25, 0.25);
                detsys_results2D[name] = new TH2D(name.c_str(), name.c_str(), variable.nbins, variable.min, variable.max, nuniverses, 0, nuniverses);
            }
        }
    }
}

// Default constructor for the DetsysCalculator class.
sys::detsys::DetsysCalculator::DetsysCalculator()
    : initialized(false)
    {}

// Accessor method for the initialized flag.
bool sys::detsys::DetsysCalculator::is_initialized()
{
    return initialized;
}

// Accessor method for the variable name.
std::string sys::detsys::DetsysCalculator::get_variable()
{
    return variable;
}

std::vector<std::string> sys::detsys::DetsysCalculator::get_variables()
{
    return variables;
}

// Accessor method for the number of universes.
size_t sys::detsys::DetsysCalculator::get_nuniverses()
{
    return nuniverses;
}

// Increment the nominal count by the specified value.
void sys::detsys::DetsysCalculator::increment_nominal_count(double value)
{
    nominal_count += value;
}

// Method to get the zscores for a given detector systematic parameter.
std::vector<double> sys::detsys::DetsysCalculator::get_zscores(std::string name)
{
    return zscores[name];
}

// Write the histogram of each configured variation and the splines for each
// detector systematic parameter to the output file.
void sys::detsys::DetsysCalculator::write()
{
    histogram_directory->cd();
    for(auto & [key, value] : histograms)
    {
        value->GetYaxis()->SetTitle("Signal Candidates / 1e18 POT");
        result_directory->WriteObject(value, key.c_str());
    }
    for(auto & [key, value] : splines)
    {
        TDirectory * tmp = result_directory->mkdir((key+"_splines").c_str());
        tmp->cd();
        for(size_t i(0); i < value.size(); ++i)
            tmp->WriteObject(value[i], (std::string("bin") + std::to_string(i)).c_str());
        histogram_directory->cd();
    }
}

// Write the result histograms for each detector systematic parameter to the
// output file.
void sys::detsys::DetsysCalculator::write_results()
{
    result_directory->cd();
    for(auto & [key, value] : detsys_results2D)
    {
        std::string name = key + "_2D";
        histogram_directory->WriteObject(value, name.c_str());
        for(int i(0); i < value->GetNbinsY(); ++i)
        {
            double sum(0);
            for(int j(0); j < value->GetNbinsX(); ++j)
                sum += value->GetBinContent(j+1, i+1);
            detsys_results1D[key]->Fill((sum - nominal_count) / nominal_count);
        }
    }
    for(auto & [key, value] : detsys_results1D)
    {
        std::string name = key + "_1D";
        histogram_directory->WriteObject(value, name.c_str());
    }
    
}

// Accessor method for the histograms.
TH1D * sys::detsys::DetsysCalculator::operator[](std::string key)
{
    return histograms[key];
}

// Method to get the weight for a given detector systematic parameter, value
// of the binning variable, and z-score.
double sys::detsys::DetsysCalculator::get_weight(std::string name, double value, double zscore)
{
    std::string spline_key = make_spline_key(name, variable);
    if(value < hdummies[spline_key]->GetXaxis()->GetXmin() || value > hdummies[spline_key]->GetXaxis()->GetXmax())
        return 1;
    else
    {
        int bin = hdummies[spline_key]->FindBin(value);
        return splines[spline_key][bin-1]->Eval(zscore);
    }
}

// Method to add a value to the detector systematic parameter histogram
// for all pre-roll z-scores (universes).
void sys::detsys::DetsysCalculator::add_value(std::string varname, double binvar, std::string detsysname, const std::map<std::string, double> & values)
{
    if(variables.size() == 1)
    {
        auto it = values.find(variable);
        if(it == values.end())
            return;

        std::string name = varname + "_" + detsysname;
        for(size_t i(0); i < nuniverses; ++i)
            detsys_results2D[name]->Fill(binvar, i, get_weight(detsysname, it->second, random_zscores[i]));
        return;
    }

    for(const std::string & configured_variable : variables)
    {
        auto value_it = values.find(configured_variable);
        if(value_it == values.end())
            continue;

        std::string spline_key = make_spline_key(detsysname, configured_variable);
        double configured_value = value_it->second;
        if(configured_value < hdummies[spline_key]->GetXaxis()->GetXmin() || configured_value > hdummies[spline_key]->GetXaxis()->GetXmax())
            continue;

        int bin = hdummies[spline_key]->FindBin(configured_value);
        std::string name = varname + "_" + detsysname + "_" + configured_variable;
        for(size_t i(0); i < nuniverses; ++i)
            detsys_results2D[name]->Fill(binvar, i, splines[spline_key][bin-1]->Eval(random_zscores[i]));
    }
}
