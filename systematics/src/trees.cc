/**
 * @file trees.cc
 * @brief Implementation file for the trees namespace.
 * @details This file contains the implementation of the functions that read
 * and interface with the TTrees produced by the CAFAna analysis framework.
 * Different "copying" actions can be performed on the TTrees, such as a simple
 * copy or adding systematics to the output file based on the selected signal
 * candidates and the configured systematics.
 * @author mueller@fnal.gov
 */
#include <cmath>
#include <iostream>

#include "trees.h"
#include "detsys.h"
#include "utilities.h"
#include "configuration.h"
#include "genie_record.h"
#include "systematic.h"
#include "weight_reader.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

// Copy the input TTree to the output TTree.
void sys::trees::copy_tree(cfg::ConfigurationTable & table, TFile * output, TFile * input)
{
    /**
     * @brief Create the output subdirectory following the nesting outlined
     * in the configuration file.
     */
    TDirectory * directory = (TDirectory *) output;
    directory = create_directory(directory, table.get_string_field("destination").c_str());
    directory->cd();

    /**
     * @brief Check if the exposure information ("POT", "Livetime") has
     * alread been copied and saved. If not, copy the exposure information
     * to the output TTree.
     */
    if(!directory->GetListOfKeys()->Contains("POT"))
    {
        std::cout << "Copying POT and Livetime histograms." << std::endl;
        TDirectory * parent = (TDirectory *) input;
        parent = get_parent_directory(parent, table.get_string_field("origin").c_str());
        TH1D * pot = (TH1D *) parent->Get("POT");
        TH1D * livetime = (TH1D *) parent->Get("Livetime");
        directory->WriteObject(pot, "POT");
        directory->WriteObject(livetime, "Livetime");
    }
    
    /**
     * @brief Create the output TTree with the name specified in the
     * configuration file.
     */
    TTree * output_tree = new TTree(table.get_string_field("name").c_str(), table.get_string_field("name").c_str());

    /**
     * @brief Connect to the input TTree and associated branches.
     * @details Three are N+3 branches in the input TTree, where N is the
     * number of branches in the input TTree of type double. The three
     * other branches (Run, Subrun, Evt) are of type int. This block uses
     * a single array to store the values of the double branches and three
     * separate variables to store the values of the int branches.
     */
    TTree * input_tree = (TTree *) input->Get(table.get_string_field("origin").c_str());
    int run, subrun, event;
    double br[input_tree->GetNbranches()-3];
    for (int i = 0; i < input_tree->GetNbranches()-3; i++)
        input_tree->SetBranchAddress(input_tree->GetListOfBranches()->At(i)->GetName(), br+i);
    input_tree->SetBranchAddress("Run", &run);
    input_tree->SetBranchAddress("Subrun", &subrun);
    input_tree->SetBranchAddress("Evt", &event);

    /**
     * @brief Create the branches in the output TTree following the same
     * structure as the input TTree.
     * @details The same array (type double) and three variables (type int)
     * as used above with the input TTree are used to create the branches
     * in the output TTree. The branches are created in the same order as
     * the input TTree. This process streamlines the copying of the values
     * from the input TTree to the output TTree.
     */
    for (int i = 0; i < input_tree->GetNbranches()-3; i++)
        output_tree->Branch(input_tree->GetListOfBranches()->At(i)->GetName(), br+i);
    output_tree->Branch("Run", &run);
    output_tree->Branch("Subrun", &subrun);
    output_tree->Branch("Evt", &event);

    /**
     * @brief Loop over the input TTree and copy the values to the output
     * TTree.
     */
    for(int i(0); i < input_tree->GetEntries(); ++i)
    {
        input_tree->GetEntry(i);
        output_tree->Fill();
    }

    /**
     * @brief Write the output TTree to the output ROOT file.
     */
    directory->WriteObject(output_tree, table.get_string_field("name").c_str());
}

// Add reweightable systematics to the output TTree.
void sys::trees::copy_with_weight_systematics(cfg::ConfigurationTable & config, cfg::ConfigurationTable & table, TFile * output, TFile * input, sys::detsys::DetsysCalculator & calc)
{
    /**
     * @brief Create the output subdirectory following the nesting outlined
     * in the configuration file.
     */
    TDirectory * directory = (TDirectory *) output;
    directory = create_directory(directory, table.get_string_field("destination").c_str());
    directory->cd();
    
    /**
     * @brief Check if the exposure information ("POT", "Livetime") has
     * alread been copied and saved. If not, copy the exposure information
     * to the output TTree.
     */
    if(!directory->GetListOfKeys()->Contains("POT"))
    {
        std::cout << "Copying POT and Livetime histograms." << std::endl;
        TDirectory * parent = (TDirectory *) input;
        parent = get_parent_directory(parent, table.get_string_field("origin").c_str());
        TH1D * pot = (TH1D *) parent->Get("POT");
        TH1D * livetime = (TH1D *) parent->Get("Livetime");
        directory->WriteObject(pot, "POT");
        directory->WriteObject(livetime, "Livetime");
    }

    /**
     * @brief Connect to the input TTree and associated branches.
     * @details Three are N+3 branches in the input TTree, where N is the
     * number of branches in the input TTree of type double. The three
     * other branches (Run, Subrun, Evt) are of type int. This block uses
     * a single array to store the values of the double branches and three
     * separate variables to store the values of the int branches. There is
     * one quirk, however, as we would also like to have access to the
     * "true_neutrino_id" branch in the input TTree directly.
     */
    TTree * input_tree = (TTree *) input->Get(table.get_string_field("origin").c_str());
    std::map<std::string, double> brs;
    double nu_id;
    Int_t run, subrun, event;
    for(int i(0); i < input_tree->GetNbranches()-3; ++i)
    {
        std::string brname = input_tree->GetListOfBranches()->At(i)->GetName();
        
        // We explicitly handle this branch, so we skip it in this loop.
        if(brname == "true_neutrino_id")
            continue;

        // Initialize the branch value to 0 and set the branch address.
        brs[brname] = 0;
        input_tree->SetBranchAddress(brname.c_str(), &brs[brname]);
    }
    input_tree->SetBranchAddress("true_neutrino_id", &nu_id);
    input_tree->SetBranchAddress("Run", &run);
    input_tree->SetBranchAddress("Subrun", &subrun);
    input_tree->SetBranchAddress("Evt", &event);

    /**
     * @brief Check whether the GENIE event records have been requested, and
     * whether the inputs can actually supply them.
     * @details The GENIE event records are only copied to the output file if the
     * "general.store_genie_evt_rec" field is set in the configuration file. The
     * feature rests on three foundations that are not present in every input:
     * the "true_neutrino_id" branch of the selected candidate tree, which
     * indexes "rec.mc.nu"; the "rec.mc.nu.genie_evtrec_idx" branch of the CAF
     * files, which indexes the GENIE event record tree; and the
     * "GenieEvtRecTree" itself. The first is checked here, and the remaining two
     * are checked below once the CAF files have been opened. If any of them is
     * missing the request is downgraded to a warning rather than an error, so
     * that a configuration file shared across samples does not fail on the
     * samples that happen to lack the records.
     */
    bool store_genie = config.get_bool_field("general.store_genie_evt_rec", false);
    if(store_genie && input_tree->GetBranch("true_neutrino_id") == nullptr)
    {
        std::cerr << "Warning: 'general.store_genie_evt_rec' is set, but the tree "
                  << table.get_string_field("origin") << " has no 'true_neutrino_id' branch. "
                  << "The GENIE event records will not be stored." << std::endl;
        store_genie = false;
    }

    /**
     * @brief Create the output TTree with the name specified in the
     * configuration file.
     * @details The output TTree is created with the same branches as the
     * input TTree, plus the Run, Subrun, and Evt branches. The same array
     * (type double) and three variables (type int) as used above with the
     * input TTree are used to create the branches in the output TTree. The
     * branches are created in the same order as the input TTree. This
     * process streamlines the copying of the values from the input TTree to
     * the output TTree.
     */
    TTree * output_tree = new TTree(table.get_string_field("name").c_str(), table.get_string_field("name").c_str());
    for(auto & br : brs)
        output_tree->Branch(br.first.c_str(), &br.second);
    output_tree->Branch("true_neutrino_id", &nu_id);
    output_tree->Branch("Run", &run);
    output_tree->Branch("Subrun", &subrun);
    output_tree->Branch("Evt", &event);

    /**
     * @brief Create an output TTree for non-matched signal candidates.
     * @details This differs from the output TTree above in that it is meant to
     * capture signal candidates that do not have a matching neutrino in the
     * input CAF files (e.g., cosmics or failed truth matching). This TTree is
     * only created if the "output.save_nonmatched" field in the configuration
     * file is set to true.
     */
    TTree * nonmatched_tree = nullptr;
    if(config.get_bool_field("output.nonmatched", true))
    {
        nonmatched_tree = new TTree((table.get_string_field("name") + "_nonmatched").c_str(), (table.get_string_field("name") + "_nonmatched").c_str());
        for(auto & br : brs)
            nonmatched_tree->Branch(br.first.c_str(), &br.second);
        nonmatched_tree->Branch("Run", &run);
        nonmatched_tree->Branch("Subrun", &subrun);
        nonmatched_tree->Branch("Evt", &event);
    }
    
    /**
     * @brief Create the map of selected signal candidates.
     * @details This block creates a map of selected signal candidates. The
     * map is built by looping over the input TTree and storing an index of the
     * run, subrun, event, nu_id, and nu_energy branches as the key. The
     * value is the index of the entry in the input TTree.
     */
    std::map<index_t, size_t> candidates;
    bool use_additional_hash = config.get_bool_field("input.use_additional_hash", false);
    for(int i(0); i < input_tree->GetEntries(); ++i)
    {
        input_tree->GetEntry(i);
        // Only consider entries with a valid neutrino ID. That is, cosmics and
        // failed truth matching will not be included in the candidates map and
        // will be copied to the non-matched TTree if it has been created.
        if(nu_id >= 0)
        {
            if(!use_additional_hash)
                candidates.insert(std::make_pair<index_t, size_t>(std::make_tuple(run, subrun, event, nu_id, 0), i));
            else
                candidates.insert(std::make_pair<index_t, size_t>(std::make_tuple(run, subrun, event, nu_id, brs["true_neutrino_energy"]), i));
        }
    }

    /**
     * @brief Configure the weight-based systematics.
     * @details This block configures the weight-based systematics. The
     * systematics are split (by type) into separate TTrees, which is
     * enforced by the "type" field in the configuration block for each
     * systematic. Because we do not wish to loop over the selected signal
     * candidates multiple times, we must store the systematic information
     * in such a way that we can easily accomodate this scheme. The variable
     * "systematics" is a map of Systematic objects keyed by the name of the
     * systematic parameter. Each Systematic object contains metadata about
     * the systematic parameter (name, index, type, etc.), some
     * configuration information, and a pointer to the output TTree,
     * weights vector, and zscores vector.
     */
    std::map<std::string, Systematic *> systematics;
    std::map<std::string, TTree *> systrees;

    /**
     * @brief Create histograms for storing the systematic results as a 
     * function of a collection of variables.
     * @details This block creates histograms for storing the systematic
     * weights / selected ratios as a function of a the variables specified
     * in the configuration file. The histograms are stored in a map with
     * the key being a pair of the variable name and the systematic name.
     * The 1D histogram contains a single entry per universe with a fill
     * value corresponding to the ratio of the selected signal candidates
     * with the universe weight to the nominal count. The 2D histogram
     * contains a 2D histogram with the variable on the x-axis and the
     * universe index on the y-axis. The fill value is the universe weight.
     * The 1D histograms can be easily inspected to see the one-bin effect
     * (uncertainty) of the systematic on the selected signal candidates.
     * The 2D histograms contain similar information, but can additionally
     * be used to inspect the effect of the systematic as a function of the
     * variable or calculate a covariance matrix.
     */
    std::vector<SysVariable> sysvariables;
    std::map<syst_t, TH2D *> results2d;
    std::map<syst_t, TH1D *> results1d;
    for(cfg::ConfigurationTable & t : config.get_subtables("sysvar"))
    {
        sysvariables.push_back(SysVariable(t));
        calc.add_variable(sysvariables.back());
    }

    /**
     * @brief Loop over the systematic types in the configuration file.
     * @details This block loops over the systematic types in the
     * configuration file. The "table" field in the configuration file
     * specifies the name of the table lists in the configuration file that
     * contain the exact definition of the systematics. Principally, this
     * loop is used to load and configure the systematics of each type in
     * sequential order.
     */
    std::vector<std::string> table_types = table.get_string_vector("table_types");
    for(const std::string & s : table_types)
    {
        std::string tname = table.get_string_field("name") + '_' + s;
        systrees[tname] = new TTree(
            (tname + "Tree").c_str(),
            (tname + "Tree").c_str());
        systrees[tname]->SetDirectory(nullptr);
        systrees[tname]->Branch("Run", &run);
        systrees[tname]->Branch("Subrun", &subrun);
        systrees[tname]->Branch("Evt", &event);
        systrees[tname]->SetDirectory(directory);
        systrees[tname]->SetAutoFlush(1000);
    }

    std::vector<double> default_clip = config.has_field("general.weight_clip")
        ? config.get_double_vector("general.weight_clip")
        : std::vector<double>{};

    for(cfg::ConfigurationTable & t : config.get_subtables("sys"))
    {
        std::string tname = table.get_string_field("name") + '_' + t.get_string_field("type");
        systematics.insert(std::make_pair<std::string, Systematic *>(t.get_string_field("name"), new Systematic(t, systrees[tname], default_clip)));
        Systematic * tmp = systematics[t.get_string_field("name")];
        tmp->get_tree()->Branch(t.get_string_field("name").c_str(), &systematics[t.get_string_field("name")]->get_weights());
        if(tmp->get_nsigma()->size() > 0)
        {
            tmp->get_tree()->Branch((t.get_string_field("name") + "_sigma").c_str(), &systematics[t.get_string_field("name")]->get_nsigma());
        }
    }

    sys::WeightReader reader(config.get_string_field("input.weights"));

    /**
     * @brief Configure the copying of the GENIE event records.
     * @details This block completes the checks begun above, now that the CAF
     * files have been opened: the reader reports whether both the
     * "rec.mc.nu.genie_evtrec_idx" branch and the "GenieEvtRecTree" that it
     * indexes were found. The records are written to a tree that is filled in
     * lockstep with the output TTree, in the manner of the systematic trees
     * above, so that entry N of the record tree belongs to entry N of the
     * selected candidate tree.
     */
    if(store_genie && !reader.has_genie_evtrec())
    {
        std::cerr << "Warning: 'general.store_genie_evt_rec' is set, but the input CAF files "
                  << "do not carry both 'rec.mc.nu.genie_evtrec_idx' and 'GenieEvtRecTree'. "
                  << "The GENIE event records will not be stored." << std::endl;
        store_genie = false;
    }

    /**
     * @brief Confirm that the GENIE classes can be written before copying any
     * record.
     * @details Writing a record requires ROOT to stream a TObject-derived class.
     * Without a compiled dictionary ROOT emulates the class well enough to read
     * it, but writing it dispatches a virtual call against an object that has no
     * vtable, which is a segmentation fault rather than a recoverable error. The
     * check is therefore made up front, and the feature is disabled rather than
     * risking a crash partway through a long job.
     */
    std::string reason;
    if(store_genie && !sys::records_are_writable(reader.get_genie_tree(), reason))
    {
        std::cerr << "Warning: 'general.store_genie_evt_rec' is set, but the GENIE event records "
                  << "cannot be written: " << reason << ". "
                  << "The GENIE event records will not be stored." << std::endl;
        store_genie = false;
    }
    sys::GenieRecordWriter genie_writer(table.get_string_field("name") + "_genieTree", directory);

    std::vector<index_t> saved_indices;
    double nominal_count(0);
    while(reader.next())
    {
        /**
         * @brief Loop over the neutrinos in the CAF input files.
         * @details This block loops over the neutrinos in the CAF input
         * files. The loop is used to populate the output TTree with the
         * selected signal candidates and the universe weights for matched
         * neutrinos. The loop also retrieves the selected signal candidate
         * that has been matched with the parent neutrino and copies the
         * values to the output TTree.
         */
        for(size_t idn(0); idn < reader.get_nnu(); ++idn)
        {
            index_t index;
            if(!use_additional_hash)
                index = std::make_tuple(reader.get_run(), reader.get_subrun(), reader.get_event(), idn, 0);
            else
                index = std::make_tuple(reader.get_run(), reader.get_subrun(), reader.get_event(), idn, (double)reader.get_energy(idn));
            if(candidates.find(index) != candidates.end())
            {
                /**
                 * @brief Retrieve the selected signal candidate and copy
                 * the values to the output TTree.
                 * @details This block retrieves the selected signal
                 * candidate that has been matched with the parent neutrino
                 * and copies the values to the output TTree.
                 */
                input_tree->GetEntry(candidates[index]);
                run = reader.get_run();
                subrun = reader.get_subrun();
                event = reader.get_event();
                calc.increment_nominal_count(1.0);
                nominal_count += 1.0;
                output_tree->Fill();

                /**
                 * @brief Store the GENIE event record for the parent neutrino.
                 * @details The loop index "idn" is the position of the parent
                 * neutrino within "rec.mc.nu", which is exactly what the
                 * "true_neutrino_id" branch of the selected candidate tree holds
                 * and what the candidates map was keyed on. The record index
                 * that it yields is relative to the CAF file that is currently
                 * loaded, so it is resolved against that file's own record tree.
                 * Exactly one entry is appended per selected candidate, so the
                 * record tree stays aligned with the output TTree even for
                 * candidates whose record is missing.
                 */
                if(store_genie)
                    genie_writer.fill(reader.get_genie_tree(), reader.get_file_index(), reader.get_genie_evtrec_idx(idn));

                /**
                 * @brief Store the universe weights in the output TTree.
                 * @details This block stores the universe weights in the
                 * output TTree for each of the configured systematics.  
                 */
                for(auto & [key, value] : systematics)
                {
                    value->get_weights()->clear();
                    if(value->get_type() == Type::kMULTISIM || value->get_type() == Type::kMULTISIGMA)
                    {
                        for(SysVariable & sv : sysvariables)
                        {
                            syst_t syskey = std::make_pair(sv.name, value->get_index());
                            reader.set(value->get_index());
                            if(results1d.find(syskey) == results1d.end())
                            {
                                results1d[syskey] = new TH1D((sv.name + "_" + key + "_1d").c_str(), (sv.name + "_" + key + "_1d").c_str(), 1000, -0.25, 0.25);
                                results1d[syskey]->SetDirectory(nullptr);
                                results2d[syskey] = new TH2D((sv.name + "_" + key + "_2d").c_str(), (sv.name + "_" + key + "_2d").c_str(), sv.nbins, sv.min, sv.max, reader.get_nuniv(idn), 0, reader.get_nuniv(idn));
                                results2d[syskey]->SetDirectory(nullptr);
                            }
                            for(size_t u(0); u < reader.get_nuniv(idn); ++u)
                            {
                                double w = value->clip(reader.get_weight(idn, u));
                                value->get_weights()->push_back(w);
                                results2d[syskey]->Fill(brs[sv.name], u, w);
                            }
                        }
                    }
                    else
                    {
                        for(double & z : calc.get_zscores(key))
                            value->get_weights()->push_back(value->clip(calc.get_weight(key, brs[calc.get_variable()], z)));
                        for(SysVariable & sv : sysvariables)
                            calc.add_value(sv.name, brs[sv.name], key, brs[calc.get_variable()]);
                    }
                } // End of loop over the configured systematics.

                /**
                 * @brief Fill the systematic TTrees.
                 * @details This block fills the systematic TTrees with
                 * the universe weights for the parent neutrino. Each
                 * configured systematic should have its weights vector
                 * populated by the above loop.
                 */
                for(auto & [key, value] : systrees)
                    value->Fill();

                // Save the index of the matched signal candidate.
                saved_indices.push_back(index);
            } // End of block for matched signal candidates.
        }
    }

    // Fill the non-matched TTree if it has been created.
    if(nonmatched_tree)
    {
        // Non-matched signal candidates that are actually neutrinos can happen
        // due to file mismatching. Though the user is expected to ensure that
        // the input TTree and the input weights file correspond to the same
        // set of events, this is not enforced by the code. Therefore, we copy
        // any neutrino entries that do not have a match in the input weights
        // file to the non-matched TTree as a "audible" sign that something is
        // amiss.
        for(auto & [key, value] : candidates)
        {
            if(std::find(saved_indices.begin(), saved_indices.end(), key) != saved_indices.end())
                continue;
            input_tree->GetEntry(value);
            run = std::get<0>(key);
            subrun = std::get<1>(key);
            event = std::get<2>(key);
            nonmatched_tree->Fill();
        }

        // The primary use case for the non-matched TTree is to capture cosmics
        // and failed truth matching. We explicitly write all entries of the
        // input tree that have a neutrino ID less than 0 to the non-matched
        // TTree to capture these cases. An unmatched reco interaction (no
        // corresponding entry in sr->dlp_true -- the common case for a
        // cosmic in an overlay sample) is given true_neutrino_id = NaN by
        // the selection framework's kNoMatchValue convention, not a
        // negative number, so it must be checked for explicitly: NaN
        // compares false against both "< 0" and ">= 0" (IEEE 754), so
        // without this check these entries would silently be dropped from
        // both this tree and the matched candidates map above.
        for(int i(0); i < input_tree->GetEntries(); ++i)
        {
            input_tree->GetEntry(i);
            if(nu_id < 0 || std::isnan(nu_id))
            {
                run = reader.get_run();
                subrun = reader.get_subrun();
                event = reader.get_event();
                nonmatched_tree->Fill();
            }
        }

        directory->WriteObject(nonmatched_tree, nonmatched_tree->GetName());
        delete nonmatched_tree;
    }

    // Write the output TTree to the output file.
    directory->WriteObject(output_tree, table.get_string_field("name").c_str());
    for(auto & [key, value] : systrees)
        directory->WriteObject(value, (key+"Tree").c_str());

    // Write the GENIE event records to the output file.
    if(store_genie)
        genie_writer.write();
    
    // Write the systematic histograms to the output file.
    std::string destination = config.get_string_field("output.histogram_destination", "");
    if(destination != "")
    {
        TDirectory * histogram_directory = create_directory(output, destination.c_str());
        for(auto & [key, value] : results2d)
        {
            std::string name = value->GetName();
            histogram_directory->WriteObject(value, name.c_str());
            for(int i(0); i < value->GetNbinsY(); ++i)
            {
                double sum(0);
                for(int j(0); j < value->GetNbinsX(); ++j)
                    sum += value->GetBinContent(j+1, i+1);
                results1d[key]->Fill((sum - nominal_count) / nominal_count);
            }
            delete value;
        }
        for(auto & [key, value] : results1d)
        {
            std::string name = value->GetName();
            histogram_directory->WriteObject(value, name.c_str());
            delete value;
        }
    }

    // Write detector systematic histograms to the output file.
    if(calc.is_initialized())
        calc.write_results();
}