/**
 * @file trees.h
 * @brief Header file for the trees namespace.
 * @details This file contains the header for the trees namespace. The trees
 * namespace contains functions that read and interface with the TTrees
 * produced by the CAFAna analysis framework. Different "copying" actions can
 * be performed on the TTrees, such as a simple copy or adding systematics to
 * the output file based on the selected signal candidates and the configured
 * systematics.
 * @author mueller@fnal.gov
 */
#ifndef TREES_H
#define TREES_H
#include <iostream>

#include "detsys.h"
#include "utilities.h"
#include "configuration.h"
#include "systematic.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"

/**
 * @namespace sys::trees
 * @brief Namespace for functions that read and interface with the TTrees
 * produced by the CAFAna analysis framework.
 * @details This namespace contains functions that read and interface with the
 * TTrees produced by the CAFAna analysis framework. Different "copying" 
 * actions can be performed on the TTrees, such as a simple copy or the
 * addition of systematics to the output file based on the selected signal
 * candidates and the configured systematics.
 */
namespace sys::trees
{
    /**
     * @brief Type definitions for the selected signal candidates indexing, the
     * universe weights, and systematic indexing (variable name and index).
     */
    typedef std::tuple<Double_t, Double_t, Double_t, Double_t> index_t;
    typedef std::map<index_t, size_t> map_t;
    typedef std::pair<std::string, int64_t> syst_t;

    /**
     * @brief Copy the input TTree to the output TTree.
     * @details This function copies the input TTree to the output TTree. The
     * function loops over the input TTree and copies the values of the branches
     * to the output TTree. The output TTree is created with the same branches
     * as the input TTree.
     * @param table The table that contains the configuration for the tree.
     * @param output The output TFile.
     * @param input The input TFile.
     * @return void
     */
    void copy_tree(sys::cfg::ConfigurationTable & table, TFile * output, TFile * input)
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

    /**
     * @brief Add reweightable systematics to the output TTree.
     * @details This function adds reweightable systematics to the output
     * TTree. The function loops over the input TTree to build a map for the
     * selected signal candidates to their index in the input TTree. The
     * function then loops over the neutrinos in the CAF input files and
     * populates the output TTree with the selected signal candidates and the
     * universe weights for matched neutrinos.
     * @param table The table that contains the configuration for the tree.
     * @param output The output TFile.
     * @param input The input TFile.
     * @return void
     */
    void copy_with_weight_systematics(sys::cfg::ConfigurationTable & config, sys::cfg::ConfigurationTable & table, TFile * output, TFile * input, sys::detsys::DetsysCalculator & calc)
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
         * "nu_id" branch in the input TTree directly.
         */
        TTree * input_tree = (TTree *) input->Get(table.get_string_field("origin").c_str());
        std::map<std::string, double> brs;
        double nu_id;
        Int_t run, subrun, event;
        for(int i(0); i < input_tree->GetNbranches()-3; ++i)
        {
            std::string brname = input_tree->GetListOfBranches()->At(i)->GetName();
            brs[brname] = 0;
            input_tree->SetBranchAddress(brname.c_str(), &brs[brname]);
        }
        input_tree->SetBranchAddress("nu_id", &nu_id);
        input_tree->SetBranchAddress("Run", &run);
        input_tree->SetBranchAddress("Subrun", &subrun);
        input_tree->SetBranchAddress("Evt", &event);


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
        output_tree->Branch("Run", &run);
        output_tree->Branch("Subrun", &subrun);
        output_tree->Branch("Evt", &event);
        
        /**
         * @brief Create the map of selected signal candidates.
         * @details This block creates a map of selected signal candidates. The
         * map is built by looping over the input TTree and storing the index
         * (run, subrun, event, nu_id) of the selected signal candidates in the
         * map. The index is used to match the selected signal candidates with
         * the universe weights for the parent neutrino.
         */
        map_t candidates;
        for(int i(0); i < input_tree->GetEntries(); ++i)
        {
            input_tree->GetEntry(i);
            candidates.insert(std::make_pair<index_t, size_t>(std::make_tuple(run, subrun, event, nu_id), i));
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
        for(sys::cfg::ConfigurationTable & t : config.get_subtables("sysvar"))
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
            systrees[s] = new TTree((s+"Tree").c_str(), (s+"Tree").c_str());
            systrees[s]->SetDirectory(nullptr);
            systrees[s]->Branch("Run", &run);
            systrees[s]->Branch("Subrun", &subrun);
            systrees[s]->Branch("Evt", &event);
            systrees[s]->SetDirectory(directory);
            systrees[s]->SetAutoFlush(1000);
        }

        for(sys::cfg::ConfigurationTable & t : config.get_subtables("sys"))
        {
            systematics.insert(std::make_pair<std::string, Systematic *>(t.get_string_field("name"), new Systematic(t, systrees[t.get_string_field("type")])));
            systematics[t.get_string_field("name")]->get_tree()->Branch(t.get_string_field("name").c_str(), &systematics[t.get_string_field("name")]->get_weights());
        }

        /**
         * @brief Load the input CAF files.
         * @details This block loads the input CAF files. The input CAF files
         * are specified in a file list that is read from the configuration
         * file. The file list is read line-by-line and the input CAF files are
         * stored in a vector.
         */
        std::vector<std::string> input_files;
        std::ifstream file_list(config.get_string_field("input.caflist"));
        std::string line;
        while(std::getline(file_list, line))
            input_files.push_back(line);
        file_list.close();

        /**
         * @brief Loop over the input CAF files.
         * @details This block loops over the input CAF files. This loop begins
         * the process of matching the selected signal candidates with the universe
         * weights for parent neutrino. 
         */
        size_t nprocessed(0);
        double nominal_count(0);
        for(std::string input_file : input_files)
        {
            if(nprocessed % 100 == 0)
                std::cout << "Processed " << nprocessed << " files." << std::endl;
            /**
             * @brief Open and validate the input CAF file.
             * @details This block opens the input CAF file and checks that the
             * files is not a zombie (corrupted) and that it contains the TTree
             * "recTree". If the file is a zombie or does not contain the TTree,
             * the code prints an error message and continues to the next file.
             * Practically, skipping a file means that signal events from that file
             * will not be included in the output ROOT file, though a situation
             * where this occurs is not expected.
             */
            TFile * caf = TFile::Open(input_file.c_str(), "READ");
            if(caf->IsZombie() || !caf->GetListOfKeys()->Contains("recTree"))
            {
                std::cerr << "Error: File" << input_file << " does not exist." << std::endl;
                continue;
            }

            /**
             * @brief Connect to the input CAF file fields.
             * @details This block connects to a minimal set of the input CAF file
             * fields. We need the run, subrun, event, and the true interaction
             * (parent neutrino) information to match the selected signal
             * candidates and retrieve the universe weights.
             */
            TTreeReader reader("recTree", caf);
            TTreeReaderValue<uint32_t> rrun(reader, "rec.hdr.run");
            TTreeReaderValue<uint32_t> rsubrun(reader, "rec.hdr.subrun");
            TTreeReaderValue<uint32_t> revt(reader, "rec.hdr.evt");
            TTreeReaderArray<caf::SRTrueInteraction> mc(reader, "rec.mc.nu");

            /**
             * @brief Loop over the events in the input CAF file.
             * @details This block loops over the events in the input CAF file. At
             * each event, the code checks if a neutrino interaction from the input
             * CAF file matches a selected signal candidate. If a match is found,
             * the code retrieves the universe weights for the parent neutrino and
             * stores them in the output TTree.
             */
            while(reader.Next())
            {
                for(const caf::SRTrueInteraction & nu : mc)
                {
                    index_t index(*rrun, *rsubrun, *revt, nu.index);
                    if(candidates.find(index) != candidates.end())
                    {
                        calc.increment_nominal_count(1.0);
                        nominal_count += 1.0;
                        /**
                         * @brief Retrieve the selected signal candidate and copy
                         * the values to the output TTree.
                         * @details This block retrieves the selected signal
                         * candidate that has been matched with the parent neutrino
                         * and copies the values to the output TTree.
                         */
                        input_tree->GetEntry(candidates[index]);
                        run = *rrun;
                        subrun = *rsubrun;
                        event = *revt;
                        output_tree->Fill();
                        
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
                                    if(results1d.find(syskey) == results1d.end())
                                    {
                                        results1d[syskey] = new TH1D((sv.name + "_" + key + "_1d").c_str(), (sv.name + "_" + key + "_1d").c_str(), 1000, -0.25, 0.25);
                                        results1d[syskey]->SetDirectory(nullptr);
                                        results2d[syskey] = new TH2D((sv.name + "_" + key + "_2d").c_str(), (sv.name + "_" + key + "_2d").c_str(), sv.nbins, sv.min, sv.max, nu.wgt[value->get_index()].univ.size(), 0, nu.wgt[value->get_index()].univ.size());
                                        results2d[syskey]->SetDirectory(nullptr);
                                    }
                                    for(size_t u(0); u < nu.wgt[value->get_index()].univ.size(); ++u)
                                    {
                                        value->get_weights()->push_back(nu.wgt[value->get_index()].univ[u]);
                                        results2d[syskey]->Fill(brs[sv.name], u, nu.wgt[value->get_index()].univ[u]);
                                    }
                                }
                            }
                            else
                            {
                                for(double & z : calc.get_zscores(key))
                                    value->get_weights()->push_back(calc.get_weight(key, brs[calc.get_variable()], z));
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

                    } // End of block for matched signal candidates.
                } // End of loop over the neutrino interactions in the input CAF file.
            } // End of loop over the events in the input CAF file.
            caf->Close();
            nprocessed++;
        } // End of loop over the input CAF files.
        directory->WriteObject(output_tree, table.get_string_field("name").c_str());
        for(auto & [key, value] : systrees)
            directory->WriteObject(value, (key+"Tree").c_str());
        
        // Write the systematic histograms to the output file.
        TDirectory * histogram_directory = create_directory(output, config.get_string_field("output.histogram_destination"));
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

        // Write detector systematic histograms to the output file.
        calc.write_results();
    }
}
#endif