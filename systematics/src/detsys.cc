/**
 * @file detsys.cc
 * @brief Implementation of the DetsysCalculator class.
 * @details This file contains the implementation of the DetsysCalculator class.
 * The DetsysCalculator class is used to calculate the weights for the detector
 * systematics using a spline interpolation of the ratio of the nominal and
 * sample histograms. Configuration of the class is done using a TOML-based
 * configuration file.
 * @see sys::cfg::ConfigurationTable 
 * @author mueller@fnal.gov
 */
#include <map>
#include <vector>
#include <random>

#include "detsys.h"
#include "configuration.h"
#include "utilities.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TSpline.h"
#include "TFile.h"

// Constructor for the DetsysCalculator class that initializes the class using
// the configuration table, the output file, and the input file. 
sys::detsys::DetsysCalculator::DetsysCalculator(sys::cfg::ConfigurationTable & table, TFile * output, TFile * input)
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
    variable = table.get_string_field("variations.variable");

    // Loop over the variations and create the histograms. A "variation" is a
    // single sample that implements some change w.r.t. the nominal sample in
    // the detector model.
    for(std::string variation : variations)
    {
        std::string pot_name = table.get_string_field("variations.origin") + variation + '/' + "POT";
        TH1D * h = (TH1D *) input->Get(pot_name.c_str());
        double pot = h->GetBinContent(1) / 1e18;

        std::string name = table.get_string_field("variations.origin") + variation + '/' + table.get_string_field("variations.tree");
        TTree * t = input->Get<TTree>(name.c_str());
        double value;
        t->SetBranchAddress(variable.c_str(), &value);
        std::vector<double> bins = table.get_double_vector("variations.bins");
        histograms[variation] = new TH1D(variation.c_str(), variation.c_str(), bins[0], bins[1], bins[2]);
        for(int i(0); i < t->GetEntries(); ++i)
        {
            t->GetEntry(i);
            histograms[variation]->Fill(value, 1.0 / pot);
        }
    }

    // Loop over the detector systematics and create the splines. A single
    // entry in the "sys" block of the configuration file with type "variation"
    // specifies a single detector systematic parameter, but in general may
    // consist of multiple variations spanning the range of the parameter.
    for(sys::cfg::ConfigurationTable & t : table.get_subtables("sys"))
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
        hdummies.insert(std::make_pair(name, new TH1D("hdummy", "hdummy", histograms[points[0]]->GetNbinsX(), histograms[points[0]]->GetXaxis()->GetXmin(), histograms[points[0]]->GetXaxis()->GetXmax())));

        // This block creates a TH2D that will be used to store the input for
        // the spline construction. The TH2D is filled with the ratio of the
        // variations to the nominal sample (across the range of the variable)
        // for each of the spline points. Each spline point is adjusted by the
        // scale factor configured in the "detsys" block.
        TH2D * h = new TH2D("tmp", "tmp", hdummies[name]->GetNbinsX(), hdummies[name]->GetXaxis()->GetXmin(), hdummies[name]->GetXaxis()->GetXmax(), points.size(), zscores[name][0], zscores[name].back());
        for(size_t i(0); i < points.size(); ++i)
        {
            hdummies[name]->Divide(histograms[points[i]], histograms[t.get_string_field("ordinate")]);
            for(int j(0); j < hdummies[name]->GetNbinsX(); ++j)
            {
                if(hdummies[name]->GetBinContent(j + 1) != 0)
                    h->SetBinContent(j + 1, i + 1,  1 + scale_factors[i] * (hdummies[name]->GetBinContent(j + 1) - 1));
                else
                    h->SetBinContent(j + 1, i + 1, 1);
            }
        }

        // This block creates the splines for the detector systematic parameter.
        // The splines are created for each bin of the dummy histogram. The
        splines.insert(std::make_pair(name, std::vector<TSpline3 *>()));
        for(int j(0); j < hdummies[name]->GetNbinsX(); ++j)
        {
            std::vector<double> x, y;
            for(size_t i(0); i < points.size(); ++i)
            {
                x.push_back(zscores[name][i]);
                y.push_back(h->GetBinContent(j + 1, i + 1));
            }
            splines[name].push_back(new TSpline3("spline", x.data(), y.data(), x.size()));
        }
    }
}

// Add a variable to the list of result histograms.
void sys::detsys::DetsysCalculator::add_variable(SysVariable & variable)
{
    // Create the TH1D and TH2D that will be used to store the results of
    // the detector systematic universes.
    for(auto & [key, value] : hdummies)
    {
        std::string name = variable.name + "_" + key;
        detsys_results1D[name] = new TH1D(name.c_str(), name.c_str(), 1000, -0.25, 0.25);
        detsys_results2D[name] = new TH2D(name.c_str(), name.c_str(), variable.nbins, variable.min, variable.max, nuniverses, 0, nuniverses);
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
    if(value < hdummies[name]->GetXaxis()->GetXmin() || value > hdummies[name]->GetXaxis()->GetXmax())
        return 1;
    else
    {
        int bin = hdummies[name]->FindBin(value);
        return splines[name][bin-1]->Eval(zscore);
    }
}

// Method to add a value to the detector systematic parameter histogram
// for all pre-roll z-scores (universes).
void sys::detsys::DetsysCalculator::add_value(std::string varname, double binvar, std::string detsysname, double value)
{
    std::string name = varname + "_" + detsysname;
    for(size_t i(0); i < nuniverses; ++i)
        detsys_results2D[name]->Fill(binvar, i, get_weight(detsysname, value, random_zscores[i]));
}
