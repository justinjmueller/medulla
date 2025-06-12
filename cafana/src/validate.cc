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
        // Common variables for both simulation-like and data-like events.
        TH1F * pot, * nevt;
        TTree * t;
        caf::StandardRecord * rec;

        /**
         * @brief Generate some basic "simulation-like" events for the
         * validation.
         * @details These events are used to test the framework's response
         * to the "simulation-like" events, which mimic the structure of
         * simulated events. The events are structured to test various logic
         * paths in the framework, including cases where interactions are not
         * matched (truth <--> reco) and cases where the interactions do
         * actually have matches. Additionally, the proxy for "is the event
         * selected" is the presence of a valid flash match.
         * 
         * - ES00: This represents the case where a reco interaction and a
         *   truth interaction are present, but are not matched. Both have
         *   valid flash matches, so they are selected.
         * 
         * - ES01: This represents the case where a reco interaction and a
         *   truth interaction are present, but are not matched. Both have
         *   no valid flash matches, so they are not selected.
         * 
         * - ES02: This represents the case where a reco interaction and a
         *   truth interaction are present, and they are matched. Both have
         *   valid flash matches, so they are selected.
         * 
         * - ES03: This represents the case where a reco interaction and a
         *   truth interaction are present, and they are matched. Both have
         *   no valid flash matches, so they are not selected.
         */
        // Open the file and initialize everything.
        TFile sim("validation_simlike.root", "RECREATE");
        pot = new TH1F("TotalPOT", "TotalPOT", 1, 0, 1);
        nevt = new TH1F("TotalEvents", "TotalEvents", 1, 0, 1);
        t = new TTree("recTree", "Standard Record Tree");
        rec = new caf::StandardRecord();
        t->Branch("rec", &rec);

        // ES00
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, {2, 2, 2, 2, 2}));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, {2, 2, 2, 2, 2}));
        write_event(rec, 1, 1, 0, pot, nevt, t);

        // ES01
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, {2, 2, 2, 2, 2}, false));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, {2, 2, 2, 2, 2}, false));
        write_event(rec, 1, 1, 1, pot, nevt, t);

        // ES02
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, {2, 2, 2, 2, 2}));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, {2, 2, 2, 2, 2}));
        pair(rec->dlp[0], rec->dlp_true[0]);
        write_event(rec, 1, 1, 2, pot, nevt, t);        

        // ES03
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, {2, 2, 2, 2, 2}, false));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, {2, 2, 2, 2, 2}, false));
        pair(rec->dlp[0], rec->dlp_true[0]);
        write_event(rec, 1, 1, 3, pot, nevt, t);

        // Write the tree and histograms to the file.
        t->Write();
        pot->Write();
        nevt->Write();
        sim.Close();

        // Clean up the allocated memory.
        delete rec;

        /**
         * @brief Generate some basic "data-like" events for the validation.
         * @details These events are used to test the framework's response to
         * the "data-like" events, which mimic the structure of reconstructed
         * events in data. These events have no truth information, but this
         * does mean that we need to test scenarios where we ask for truth
         * information in the branches (as one might do when reusing the same
         * selection logic for both simulation and data).
         * 
         * - ED00: This represents a reco interaction with a valid flash match,
         *   which should be selected.
         * 
         * - ED01: This represents a reco interaction with no valid flash match,
         *   which should not be selected.
         */
        // Open the file and initialize everything.
        TFile data("validation_datalike.root", "RECREATE");
        pot = new TH1F("TotalPOT", "TotalPOT", 1, 0, 1);
        nevt = new TH1F("TotalEvents", "TotalEvents", 1, 0, 1);
        t = new TTree("recTree", "Standard Record Tree");
        rec = new caf::StandardRecord();
        t->Branch("rec", &rec);

        // ED00
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, {2, 2, 2, 2, 2}));
        write_event(rec, 1, 1, 0, pot, nevt, t);

        // ED01
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, {2, 2, 2, 2, 2}, false));
        write_event(rec, 1, 1, 1, pot, nevt, t);

        // Write the tree and histograms to the file.
        t->Write();
        pot->Write();
        nevt->Write();
        data.Close();

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

        std::cout << "\033[1m--- Running validation ---\033[0m" << std::endl;
        /**
         * @brief The first set of events to validate is the "sim-like" events
         * and the response of the framework when run over them in "reco" mode.
         * @details This set of events is used to validate the framework's
         * response to the "sim-like" events, which mimic the structure of
         * simulated events. There are a few different conditions that we are 
         * looking for in the validation:
         * 
         * - SR00: This represents a reco interaction with no valid truth match
         *   under an additional truth cut that would otherwise pass the reco-
         *   only selection. This does not pass the selection.
         * 
         * - SR01: This represents a reco interaction with no valid truth match
         *   under an additional truth cut that would also not pass the
         *   reco-only selection. This does not pass the selection.
         * 
         * - SR02: This represents a reco interaction with no valid truth match
         *   under no additional truth cut and that passes the reco-only
         *   selection. This passes the selection with a valid reco-var.
         * 
         * - SR03: This represents a reco interaction with no valid truth match
         *   under no additional truth cut and that passes the reco-only
         *   selection. This passes the selection with a NaN truth-var.
         * 
         * - SR04: This represents a reco interaction with no valid truth match
         *   under no additional truth cut and that does not pass the reco-only
         *   selection. This does not pass the selection.
         * 
         * - SR05: This represents a reco interaction with a valid truth match
         *   under an additional truth cut that would also pass the reco-only
         *   selection. This passes the selection with a valid reco-var.
         * 
         * - SR06: This represents a reco interaction with a valid truth match
         *   under an additional truth cut that would also pass the reco-only
         *   selection. This passes the selection with a valid truth-var.
         * 
         * - SR07: This represents a reco interaction with a valid truth match
         *   under and additional truth cut that would also not pass the
         *   reco-only selection. This does not pass the selection.
         * 
         * - SR08: This represents a reco interaction with a valid truth match
         *   under no additional truth cut and that passes the reco-only
         *   selection. This passes the selection with a valid reco-var.
         * 
         * - SR09: This represents a reco interaction with a valid truth match
         *   under no additional truth cut and that passes the reco-only
         *   selection. This passes the selection with a valid truth-var.
         * 
         * - SR10: This represents a reco interaction with a valid truth match
         *   under no additional truth cut and that does not pass the reco-only
         *   selection. This does not pass the selection.
         */
        std::cout << "\n\033[1mSimulation-like events with mode == 'reco' \033[0m" << std::endl;

        //  Read the event data from the TTree in the ROOT file.
        std::vector<row_t> rows = read_event_data("events/test_simlike/test_reco");

        // Expected results for validation.
        std::vector<condition_t> conditions = {
            {"SR02", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}, {"reco_vertex_x", -210.0}}},
            {"SR03", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}, {"true_vertex_x", kNaN}}},
            {"!SR04", {{"Run", 1}, {"Subrun", 1}, {"Evt", 1}}},
            {"SR08", {{"Run", 1}, {"Subrun", 1}, {"Evt", 2}, {"reco_vertex_x", -210.0}}},
            {"SR09", {{"Run", 1}, {"Subrun", 1}, {"Evt", 2}, {"true_vertex_x", -210.0}}},
            {"!SR10", {{"Run", 1}, {"Subrun", 1}, {"Evt", 3}}},
        };

        // Check if each condition_t entry is present in the rows vector.
        match_conditions(rows, conditions);

        //  Read the event data from the TTree in the ROOT file.
        rows = read_event_data("events/test_simlike/test_reco_with_truth_cut");

        // Expected results for validation.
        conditions = {
            {"!SR00", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}}},
            {"!SR01", {{"Run", 1}, {"Subrun", 1}, {"Evt", 1}}},
            {"SR05", {{"Run", 1}, {"Subrun", 1}, {"Evt", 2}, {"reco_vertex_x", -210.0}}},
            {"SR06", {{"Run", 1}, {"Subrun", 1}, {"Evt", 2}, {"true_vertex_x", -210.0}}},
            {"!SR07", {{"Run", 1}, {"Subrun", 1}, {"Evt", 3}}},
        };

        // Check if each condition_t entry is present in the rows vector.
        match_conditions(rows, conditions);

        /**
         * @brief The second set of events to validate is the "data-like" events
         * and the response of the framework when run over them in "reco" mode.
         * @details This set of events is used to validate the framework's
         * response to the "data-like" events, which mimic the structure of
         * reconstructed events in data. These events have no truth information,
         * but this does mean that we need to test scenarios where we ask for
         * truth information in the branches.
         * 
         * - DR00: This represents a reco interaction with both a truth cut and
         *   the usual selection cut placed. The interaction should pass the
         *   normal selection cut, and the fact that the sample is labeled as
         *   data should bypass the truth cut. Thus, this should have a valid
         *   reco-var.
         * 
         * - DR01: This represents a reco interaction with both a truth cut and
         *   the usual selection cut placed. The interaction should pass the
         *   normal selection cut, and the fact that the sample is labeled as
         *   data should bypass the truth cut. Thus, this should have a NaN
         *   truth-var.
         * 
         * - DR02: This represents a reco interaction with both a truth cut and
         *   the usual selection cut placed. The interaction should not pass
         *   the normal selection cut, and therefore should not pass the
         *   selection.
         * 
         * - DR03: This represents a reco interaction with no truth cut and the
         *   usual selection cut placed. The interaction should pass the normal
         *   selection cut, and therefore should have a valid reco-var.
         * 
         * - DR04: This represents a reco interaction with no truth cut and the
         *   usual selection cut placed. The interaction should pass the normal
         *   selection cut, and therefore should have a NaN truth-var.  
         */
        std::cout << "\n\033[1mData-like events with mode == 'reco' \033[0m" << std::endl;

        // Read the event data from the TTree in the ROOT file.
        rows = read_event_data("events/test_datalike/test_reco");
        
        // Expected results for validation.
        conditions = {
            {"DR00", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}, {"reco_vertex_x", -210.0}}},
            {"DR01", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}, {"true_vertex_x", kNaN}}},
            {"!DR02", {{"Run", 1}, {"Subrun", 1}, {"Evt", 1}}},
            {"DR03", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}, {"reco_vertex_x", -210.0}}},
            {"DR04", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}, {"true_vertex_x", kNaN}}},
            {"!DR05", {{"Run", 1}, {"Subrun", 1}, {"Evt", 1}}},
        };

        // Check if each condition_t entry is present in the rows vector.
        match_conditions(rows, conditions);

        // Finished!
        std::cout << "\n\033[1m---        DONE        ---\033[0m" << std::endl;
        f.Close();
    }
    return 0;
}