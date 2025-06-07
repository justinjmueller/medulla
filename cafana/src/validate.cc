/**
 * @file validate.cc
 * @brief Validation execution file for the SPINE analysis framework.
 * @details This file contains the main function for testing the SPINE
 * analysis framework. It generates a ROOT file with a TTree containing
 * interactions and particles, which is used to test the framework's
 * logic and functionality.
 * @author mueller@fnal.gov
 */
#include <iostream>

#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"
#include "sbnanaobj/StandardRecord/SRParticleDLP.h"
#include "sbnanaobj/StandardRecord/SRParticleTruthDLP.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

#include "test.h"

/**
 * @brief Main function for the validation code.
 * @details This function serves two purposes: generating the structured CAF
 * file input for the framework testing and validating the output of the
 * framework against the expected results. These modes are toggled by
 * the presence of different command line flags.
 * @param argc The number of command line arguments.
 * @param argv The command line arguments. The options are:
 * - `--generate`: Generate the structured CAF file input for the framework
 *                 testing.
 * - `--validate`: Validate the output of the framework against the expected
 *                 results.
 * @return int The exit code of the program. Returns 0 on success, non-zero
 * on failure.
 */
int main(int argc, char * argv[])
{
    // Check if the command line arguments are valid.
    if(argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " --generate | --validate" << std::endl;
        return 1;
    }

    // Check the command line arguments for the mode.
    std::string mode = argv[1];
    if(mode != "--generate" && mode != "--validate")
    {
        std::cerr << "Invalid mode: " << mode << ". Use --generate or --validate." << std::endl;
        return 1;
    }

    // If the mode is generate, we create a ROOT file with a TTree containing
    // interactions and particles, which is used to test the framework's logic
    // and functionality.
    if(mode == "--generate")
    {
        TFile f("validation.root", "RECREATE");

        TH1F * pot = new TH1F("TotalPOT", "TotalPOT", 1, 0, 1);
        TH1F * nevt = new TH1F("TotalEvents", "TotalEvents", 1, 0, 1);

        TTree * t = new TTree("recTree", "Standard Record Tree");

        caf::StandardRecord * rec = new caf::StandardRecord();
        t->Branch("rec", &rec);

        // Basic interaction with particles (No interaction matches).
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, {2, 2, 2, 2, 2}));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, {2, 2, 2, 2, 2}));
        write_event(rec, 1, 1, 0, pot, nevt, t);

        // Basic interaction with particles (Reco -> True match only).
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(1, 0, {2, 2, 2, 2, 2}));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(1, 0, {2, 2, 2, 2, 2}));
        pair(rec->dlp[1], rec->dlp_true[1]);
        write_event(rec, 1, 1, 1, pot, nevt, t);

        // Write the tree and histograms to the file.
        t->Write();
        pot->Write();
        nevt->Write();
        f.Close();

        // Clean up the allocated memory.
        delete rec;

        return 0;
    }

    // If the mode is validate, we run the validation logic.
    if(mode == "--validate")
    {
        // First, we load connect to the ROOT file containing the
        // results of the validation.
        TFile f("test.root", "READ");
        if(!f.IsOpen())
        {
            std::cerr << "Error: Could not open the file 'test.root'." << std::endl;
            return 1;
        }

        // Retrieve the TTree from the file.
        TTree * t = static_cast<TTree *>(f.Get("events/test/test_reco"));
        if(!t)
        {
            std::cerr << "Error: Could not find the TTree 'test_reco' in the file." << std::endl;
            f.Close();
            return 1;
        }

        // Connect to the branches of the TTree.
        Int_t run, subrun, event_num;
        double id, true_vtxx, reco_vtxx;
        t->SetBranchAddress("Run", &run);
        t->SetBranchAddress("Subrun", &subrun);
        t->SetBranchAddress("Evt", &event_num);
        t->SetBranchAddress("reco_interaction_id", &id);
        t->SetBranchAddress("true_vertex_x", &true_vtxx);
        t->SetBranchAddress("reco_vertex_x", &reco_vtxx);

        // Fill the rows vector with the data from the TTree.
        std::vector<row_t> rows;
        for(int i = 0; i < t->GetEntries(); ++i)
        {
            t->GetEntry(i);
            row_t row;
            row["run"] = run;
            row["subrun"] = subrun;
            row["event_num"] = event_num;
            row["interaction_id"] = id;
            row["true_vertex_x"] = true_vtxx;
            row["reco_vertex_x"] = reco_vtxx;
            rows.push_back(row);
        }

        // Expected results for validation.
        std::vector<condition_t> expected_results = {
            {"Condition #0", {{"run", 1}, {"subrun", 1}, {"event_num", 0}, {"interaction_id", 0}, {"true_vertex_x", -210.0}, {"reco_vertex_x", -210.0}}},
            {"Condition #1", {{"run", 1}, {"subrun", 1}, {"event_num", 1}, {"interaction_id", 1}, {"true_vertex_x", -210.0}, {"reco_vertex_x", -210.0}}},
        };
        auto match_metadata = [](const row_t & row, const condition_t & condition) {
            for(const auto & field : condition.second)
            {
                if(field.first == "run" || field.first == "subrun" || field.first == "event_num")
                {
                    if(row.at(field.first) != field.second)
                        return false; // Metadata mismatch.
                }
            }
            return true;
        };

        // Check if each condition_t entry is present in the rows vector.
        std::cout << "\033[1m--- Running validation ---\033[0m" << std::endl;
        for(const auto & condition : expected_results)
        {
            bool found = false;
            for(const auto & row : rows)
            {
                if(match_metadata(row, condition) && row == condition.second)
                {
                    std::cout << "\033[32mValidation passed:\033[0m   " << condition.first << "." << std::endl;
                    found = true;
                    break;
                }
                if(match_metadata(row, condition) && row != condition.second)
                {
                    found = true;
                    std::cout << "\033[33mValidation mismatch:\033[0m " << condition.first << "." << std::endl;
                    
                    // Print the fields that are mismatched.
                    for(const auto & field : condition.second)
                    {
                        if(row.at(field.first) != field.second)
                        {
                            std::cout << "    " << field.first
                                      << " - expected: " << field.second
                                      << ", got: " << row.at(field.first) << std::endl;
                        }
                    }
                }
            }
            if(!found)
                std::cout << "\033[31mValidation failed:\033[0m   " << condition.first << "." << std::endl;
        }
        std::cout << "\033[1m---        DONE        ---\033[0m" << std::endl;

        f.Close();
    }

    return 0;
}