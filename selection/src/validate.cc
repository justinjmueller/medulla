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
         *   valid flash matches, so they are selected. ES00A also has a
         *   trigger time of 2000, so we can test the global trigger time
         *   event-level cut.
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

        // Default multiplicity for the final state particles (1 photon)
        const multiplicity_t fs = {1, 0, 0, 0, 0}; 

        // ES00A
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, fs));
        write_event(rec, 1, 0, 0, pot, nevt, t, 2000);

        // ES00B
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, fs));
        mark_contained(&rec->dlp[0], &rec->dlp_true[0]);
        write_event(rec, 1, 1, 0, pot, nevt, t);

        // ES00C
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, fs));
        pair(rec->dlp[0].particles[0], rec->dlp_true[0].particles[0]);
        write_event(rec, 1, 2, 0, pot, nevt, t);

        // ES00D
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, fs));
        pair(rec->dlp[0].particles[0], rec->dlp_true[0].particles[0]);
        mark_contained(&rec->dlp[0], &rec->dlp_true[0]);
        write_event(rec, 1, 3, 0, pot, nevt, t);

        // ES01A
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs, false));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, fs, false));
        write_event(rec, 1, 0, 1, pot, nevt, t);

        // ES01B
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs, false));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, fs, false));
        mark_contained(&rec->dlp[0], &rec->dlp_true[0]);
        write_event(rec, 1, 1, 1, pot, nevt, t);

        // ES01C
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs, false));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, fs, false));
        pair(rec->dlp[0].particles[0], rec->dlp_true[0].particles[0]);
        write_event(rec, 1, 2, 1, pot, nevt, t);

        // ES01D
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs, false));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, fs, false));
        pair(rec->dlp[0].particles[0], rec->dlp_true[0].particles[0]);
        mark_contained(&rec->dlp[0], &rec->dlp_true[0]);
        write_event(rec, 1, 3, 1, pot, nevt, t);

        // ES02A
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, fs));
        pair(rec->dlp[0], rec->dlp_true[0]);
        write_event(rec, 1, 0, 2, pot, nevt, t);

        // ES02B
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, fs));
        pair(rec->dlp[0], rec->dlp_true[0]);
        mark_contained(&rec->dlp[0], &rec->dlp_true[0]);
        write_event(rec, 1, 1, 2, pot, nevt, t);

        // ES02C
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, fs));
        pair(rec->dlp[0], rec->dlp_true[0]);
        pair(rec->dlp[0].particles[0], rec->dlp_true[0].particles[0]);
        write_event(rec, 1, 2, 2, pot, nevt, t);

        // ES02D
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, fs));
        pair(rec->dlp[0], rec->dlp_true[0]);
        pair(rec->dlp[0].particles[0], rec->dlp_true[0].particles[0]);
        mark_contained(&rec->dlp[0], &rec->dlp_true[0]);
        write_event(rec, 1, 3, 2, pot, nevt, t);

        // ES03A
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs, false));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, fs, false));
        pair(rec->dlp[0], rec->dlp_true[0]);
        write_event(rec, 1, 0, 3, pot, nevt, t);

        // ES03B
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs, false));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, fs, false));
        pair(rec->dlp[0], rec->dlp_true[0]);
        mark_contained(&rec->dlp[0], &rec->dlp_true[0]);
        write_event(rec, 1, 1, 3, pot, nevt, t);
        
        // ES03C
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs, false));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, fs, false));
        pair(rec->dlp[0], rec->dlp_true[0]);
        pair(rec->dlp[0].particles[0], rec->dlp_true[0].particles[0]);
        write_event(rec, 1, 2, 3, pot, nevt, t);
        
        // ES03D
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs, false));
        rec->dlp_true.push_back(generate_interaction<caf::SRInteractionTruthDLP>(0, 0, fs, false));
        pair(rec->dlp[0], rec->dlp_true[0]);
        pair(rec->dlp[0].particles[0], rec->dlp_true[0].particles[0]);
        mark_contained(&rec->dlp[0], &rec->dlp_true[0]);
        write_event(rec, 1, 3, 3, pot, nevt, t);

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

        // ED00A
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs));
        write_event(rec, 1, 0, 0, pot, nevt, t);

        // ED00B
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs));
        mark_contained(&rec->dlp[0]);
        write_event(rec, 1, 1, 0, pot, nevt, t);

        // ED01A
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs, false));
        write_event(rec, 1, 0, 1, pot, nevt, t);

        // ED01B
        rec->dlp.push_back(generate_interaction<caf::SRInteractionDLP>(0, 0, fs, false));
        mark_contained(&rec->dlp[0]);
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
            {"SR02", {{"Run", 1}, {"Subrun", 0}, {"Evt", 0}, {"reco_vertex_x", -210.0}}},
            {"SR03", {{"Run", 1}, {"Subrun", 0}, {"Evt", 0}, {"true_vertex_x", kNaN}}},
            {"!SR04", {{"Run", 1}, {"Subrun", 0}, {"Evt", 1}}},
            {"SR08", {{"Run", 1}, {"Subrun", 0}, {"Evt", 2}, {"reco_vertex_x", -210.0}}},
            {"SR09", {{"Run", 1}, {"Subrun", 0}, {"Evt", 2}, {"true_vertex_x", -210.0}}},
            {"!SR10", {{"Run", 1}, {"Subrun", 0}, {"Evt", 3}}},
        };

        // Check if each condition_t entry is present in the rows vector.
        match_conditions(rows, conditions);

        //  Read the event data from the TTree in the ROOT file.
        rows = read_event_data("events/test_simlike/test_reco_with_truth_cut");

        // Expected results for validation.
        conditions = {
            {"!SR00", {{"Run", 1}, {"Subrun", 0}, {"Evt", 0}}},
            {"!SR01", {{"Run", 1}, {"Subrun", 0}, {"Evt", 1}}},
            {"SR05", {{"Run", 1}, {"Subrun", 0}, {"Evt", 2}, {"reco_vertex_x", -210.0}}},
            {"SR06", {{"Run", 1}, {"Subrun", 0}, {"Evt", 2}, {"true_vertex_x", -210.0}}},
            {"!SR07", {{"Run", 1}, {"Subrun", 0}, {"Evt", 3}}},
        };

        // Check if each condition_t entry is present in the rows vector.
        match_conditions(rows, conditions);

        /**
         * @brief The second set of events to validate is the "sim-like" events
         * and the response of the framework when run over them in "truth" mode.
         * @details This set of events is used to validate the framework's
         * response to the "sim-like" events, which mimic the structure of
         * simulated events. These events have truth information, and
         * therefore we can test the framework's response to the truth
         * information in the branches.
         * 
         * - ST00: This represents a truth interaction with no valid reco match
         *   under an additional reco cut that would otherwise pass the truth-
         *   only selection. This does not pass the selection.
         * 
         * - ST01: This represents a truth interaction with no valid reco match
         *   under an additional reco cut that would also not pass the truth-
         *   only selection. This does not pass the selection.
         * 
         * - ST02: This represents a truth interaction with no valid reco match
         *   under no additional reco cut and that passes the truth-only
         *   selection. This passes the selection with a valid truth-var.
         * 
         * - ST03: This represents a truth interaction with no valid reco match
         *   under no additional reco cut and that passes the truth-only
         *   selection. This passes the selection with a NaN reco-var.
         * 
         * - ST04: This represents a truth interaction with no valid reco match
         *   under no additional reco cut and that does not pass the truth-only
         *   selection. This does not pass the selection.
         * 
         * - ST05: This represents a truth interaction with a valid reco match
         *   under an additional reco cut that would also pass the truth-only
         *   selection. This passes the selection with a valid truth-var.
         * 
         * - ST06: This represents a truth interaction with a valid reco match
         *   under an additional reco cut that would also pass the truth-only
         *   selection. This passes the selection with a valid reco-var.
         * 
         * - ST07: This represents a truth interaction with a valid reco match
         *   under and additional reco cut that would also not pass the truth-
         *   only selection. This does not pass the selection.
         * 
         * - ST08: This represents a truth interaction with a valid reco match
         *   under no additional reco cut and that passes the truth-only
         *   selection. This passes the selection with a valid truth-var.
         * 
         * - ST09: This represents a truth interaction with a valid reco match
         *   under no additional reco cut and that passes the truth-only
         *   selection. This passes the selection with a valid reco-var.
         * 
         * - ST10: This represents a truth interaction with a valid reco match
         *   under no additional reco cut and that does not pass the truth-only
         *   selection. This does not pass the selection.
         */
        std::cout << "\n\033[1mSimulation-like events with mode == 'truth' \033[0m" << std::endl;
        
        // Read the event data from the TTree in the ROOT file.
        rows = read_event_data("events/test_simlike/test_truth");
        
        // Expected results for validation.
        conditions = {
            {"ST02", {{"Run", 1}, {"Subrun", 0}, {"Evt", 0}, {"true_vertex_x", -210.0}}},
            {"ST03", {{"Run", 1}, {"Subrun", 0}, {"Evt", 0}, {"reco_vertex_x", kNaN}}},
            {"!ST04", {{"Run", 1}, {"Subrun", 0}, {"Evt", 1}}},
            {"ST08", {{"Run", 1}, {"Subrun", 0}, {"Evt", 2}, {"true_vertex_x", -210.0}}},
            {"ST09", {{"Run", 1}, {"Subrun", 0}, {"Evt", 2}, {"reco_vertex_x", -210.0}}},
            {"!ST10", {{"Run", 1}, {"Subrun", 0}, {"Evt", 3}}},
        };

        // Check if each condition_t entry is present in the rows vector.
        match_conditions(rows, conditions);

        // Read the event data from the TTree in the ROOT file.
        rows = read_event_data("events/test_simlike/test_truth_with_reco_cut");

        // Expected results for validation.
        conditions = {
            {"!ST00", {{"Run", 1}, {"Subrun", 0}, {"Evt", 0}}},
            {"!ST01", {{"Run", 1}, {"Subrun", 0}, {"Evt", 1}}},
            {"ST05", {{"Run", 1}, {"Subrun", 0}, {"Evt", 2}, {"true_vertex_x", -210.0}}},
            {"ST06", {{"Run", 1}, {"Subrun", 0}, {"Evt", 2}, {"reco_vertex_x", -210.0}}},
            {"!ST07", {{"Run", 1}, {"Subrun", 0}, {"Evt", 3}}},
        };

        // Check if each condition_t entry is present in the rows vector.
        match_conditions(rows, conditions);

        /**
         * @brief The third set of events to validate is the "data-like" events
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
            {"DR00", {{"Run", 1}, {"Subrun", 0}, {"Evt", 0}, {"reco_vertex_x", -210.0}}},
            {"DR01", {{"Run", 1}, {"Subrun", 0}, {"Evt", 0}, {"true_vertex_x", kNaN}}},
            {"!DR02", {{"Run", 1}, {"Subrun", 0}, {"Evt", 1}}},
            {"DR03", {{"Run", 1}, {"Subrun", 0}, {"Evt", 0}, {"reco_vertex_x", -210.0}}},
            {"DR04", {{"Run", 1}, {"Subrun", 0}, {"Evt", 0}, {"true_vertex_x", kNaN}}},
            {"!DR05", {{"Run", 1}, {"Subrun", 0}, {"Evt", 1}}},
        };

        // Check if each condition_t entry is present in the rows vector.
        match_conditions(rows, conditions);

        /**
         * @brief The fourth set of events to validate is the "sim-like" events
         * and the response of the framework when run over them in "reco" mode
         * with particle-level variables.
         * @details This set of events effectively tests the framework's
         * behavior when run in a mode where the selection logic is applied
         * to the particles of the interactions, rather than the interactions
         * themselves. Effectively, the framework implements an additional
         * nested loop over the particles (of the same type) for each
         * interaction, and applies some selection logic that is similar to the
         * one used for the interactions.
         * 
         * - SPR00: This represents a reco interaction with no valid truth match
         *   under no additional truth cut that would otherwise not pass the
         *   reco-selection. Furthermore, the particle in this interaction has
         *   no truth match and does not pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR01: This represents a reco interaction with no valid truth match
         *   under no additional truth cut that would otherwise not pass the
         *   reco-selection. Furthermore, the particle in this interaction has
         *   no truth match and does pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR02: This represents a reco interaction with no valid truth match
         *   under no additional truth cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has no
         *   truth match and does not pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR03: This represents a reco interaction with no valid truth match
         *   under no additional truth cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has no
         *   truth match and does pass the reco cut. Passes the selection with
         *   a valid reco-var.
         * 
         * - SPR04: This represents a reco interaction with no valid truth match
         *   under no additional truth cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has no
         *   truth match and does pass the reco cut. Passes the selection with
         *   a NaN truth-var.
         * 
         * - SPR05: This represents a reco interaction with no valid truth match
         *   under an additional truth cut that would otherwise not pass the
         *   reco-selection. Furthermore, the particle in this interaction has
         *   no truth match and does not pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR06: This represents a reco interaction with no valid truth match
         *   under an additional truth cut that would otherwise not pass the
         *   reco-selection. Furthermore, the particle in this interaction has
         *   no truth match and does pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR07: This represents a reco interaction with no valid truth match
         *   under an additional truth cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has no
         *   truth match and does not pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR08: This represents a reco interaction with no valid truth match
         *   under an additional truth cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has no
         *   truth match and does pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR09: This represents a reco interaction with no valid truth match
         *   under no additional truth cut that would otherwise not pass the
         *   reco-selection. Furthermore, the particle in this interaction has a
         *   truth match and does not pass the reco cut. Does not pass the
         *   selection.
         *  
         * - SPR10: This represents a reco interaction with no valid truth match
         *   under no additional truth cut that would otherwise not pass the
         *   reco-selection. Furthermore, the particle in this interaction has a
         *   truth match and does pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR11: This represents a reco interaction with no valid truth match
         *   under no additional truth cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has a
         *   truth match and does not pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR12: This represents a reco interaction with no valid truth match
         *   under no additional truth cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has a
         *   truth match and does pass the reco cut. Passes the selection with
         *   a valid reco-var.
         * 
         * - SPR13: This represents a reco interaction with no valid truth match
         *   under no additional truth cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has a
         *   truth match and does pass the reco cut. Passes the selection with
         *   a valid truth-var.
         * 
         * - SPR14: This represents a reco interaction with no valid truth match
         *   under an additional truth cut that would otherwise not pass the
         *   reco-selection. Furthermore, the particle in this interaction has a
         *   truth match and does not pass the reco cut. Does not pass the
         *   selection.
         *  
         * - SPR15: This represents a reco interaction with no valid truth match
         *   under an additional truth cut that would otherwise not pass the
         *   reco-selection. Furthermore, the particle in this interaction has a
         *   truth match and does pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR16: This represents a reco interaction with no valid truth match
         *   under an additional truth cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has a
         *   truth match and does not pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR17: This represents a reco interaction with no valid truth match
         *   under an additional truth cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has a
         *   truth match and does pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR18: This represents a reco interaction with a valid truth match
         *   under no additional truth cut and that does not pass the reco-
         *   selection. Furthermore, the particle in this interaction has no
         *   truth match and does not pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR19: This represents a reco interaction with a valid truth match
         *   under no additional truth cut and that does not pass the reco-
         *   selection. Furthermore, the particle in this interaction has no
         *   truth match and does pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR20: This represents a reco interaction with a valid truth match
         *   under no additional truth cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has no
         *   truth match and does not pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR21: This represents a reco interaction with a valid truth match
         *   under no additional truth cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has no
         *   truth match and does pass the reco cut. Passes the selection with
         *   a valid reco-var.
         * 
         * - SPR22: This represents a reco interaction with a valid truth match
         *   under no additional truth cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has no
         *   truth match and does pass the reco cut. Passes the selection with
         *   a NaN truth-var.
         * 
         * - SPR23: This represents a reco interaction with a valid truth match
         *   under an additional truth cut and that does not pass the reco-
         *   selection. Furthermore, the particle in this interaction has no
         *   truth match and does not pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR24: This represents a reco interaction with a valid truth match
         *   under an additional truth cut and that does not pass the reco-
         *   selection. Furthermore, the particle in this interaction has no
         *   truth match and does pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR25: This represents a reco interaction with a valid truth match
         *   under an additional truth cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has no
         *   truth match and does not pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR26: This represents a reco interaction with a valid truth match
         *   under an additional truth cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has no
         *   truth match and does pass the reco cut. Passes the selection with
         *   a valid reco-var.
         * 
         * - SPR27: This represents a reco interaction with a valid truth match
         *   under an additional truth cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has no
         *   truth match and does pass the reco cut. Passes the selection with
         *   a NaN truth-var.
         * 
         * - SPR28: This represents a reco interaction with a valid truth match
         *   under no additional truth cut and that does not pass the reco-
         *   selection. Furthermore, the particle in this interaction has a
         *   truth match and does not pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR29: This represents a reco interaction with a valid truth match
         *   under no additional truth cut and that does not pass the reco-
         *   selection. Furthermore, the particle in this interaction has a
         *   truth match and does pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR30: This represents a reco interaction with a valid truth match
         *   under no additional truth cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has a
         *   truth match and does not pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR31: This represents a reco interaction with a valid truth match
         *   under no additional truth cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has a
         *   truth match and does pass the reco cut. Passes the selection with
         *   a valid reco-var.
         * 
         * - SPR32: This represents a reco interaction with a valid truth match
         *   under no additional truth cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has a
         *   truth match and does pass the reco cut. Passes the selection with
         *   a valid truth-var.
         * 
         * - SPR33: This represents a reco interaction with a valid truth match
         *   under an additional truth cut and that does not pass the reco-
         *   selection. Furthermore, the particle in this interaction has a
         *   truth match and does not pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR34: This represents a reco interaction with a valid truth match
         *   under an additional truth cut and that does not pass the reco-
         *   selection. Furthermore, the particle in this interaction has a
         *   truth match and does pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR35: This represents a reco interaction with a valid truth match
         *   under an additional truth cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has a
         *   truth match and does not pass the reco cut. Does not pass the
         *   selection.
         * 
         * - SPR36: This represents a reco interaction with a valid truth match
         *   under an additional truth cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has a
         *   truth match and does pass the reco cut. Passes the selection with
         *   a valid reco-var.
         * 
         * - SPR37: This represents a reco interaction with a valid truth match
         *   under an additional truth cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has a
         *   truth match and does pass the reco cut. Passes the selection with
         *   a valid truth-var.
         */
        std::cout << "\n\033[1mSimulation-like events with mode == 'reco' and particle-level variables \033[0m" << std::endl;
        
        // Read the event data from the TTree in the ROOT file.
        rows = read_event_data("events/test_simlike/test_reco_particles");

        // Expected results for validation.
        conditions = {
            {"!SPR00", {{"Run", 1}, {"Subrun", 0}, {"Evt", 1}}},
            {"!SPR01", {{"Run", 1}, {"Subrun", 1}, {"Evt", 1}}},
            {"!SPR02", {{"Run", 1}, {"Subrun", 0}, {"Evt", 0}}},
            {"SPR03", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}, {"reco_particle_ke", 200.0}}},
            {"SPR04", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}, {"true_particle_ke", kNaN}}},
            {"!SPR09", {{"Run", 1}, {"Subrun", 2}, {"Evt", 1}}},
            {"!SPR10", {{"Run", 1}, {"Subrun", 3}, {"Evt", 1}}},
            {"!SPR11", {{"Run", 1}, {"Subrun", 2}, {"Evt", 0}}},
            {"SPR12", {{"Run", 1}, {"Subrun", 3}, {"Evt", 0}, {"reco_particle_ke", 200.0}}},
            {"SPR13", {{"Run", 1}, {"Subrun", 3}, {"Evt", 0}, {"true_particle_ke", 200.0}}},
            {"!SPR18", {{"Run", 1}, {"Subrun", 0}, {"Evt", 3}}},
            {"!SPR19", {{"Run", 1}, {"Subrun", 1}, {"Evt", 3}}},
            {"!SPR20", {{"Run", 1}, {"Subrun", 0}, {"Evt", 2}}},
            {"SPR21", {{"Run", 1}, {"Subrun", 1}, {"Evt", 2}, {"reco_particle_ke", 200.0}}},
            {"SPR22", {{"Run", 1}, {"Subrun", 1}, {"Evt", 2}, {"true_particle_ke", kNaN}}},
            {"!SPR28", {{"Run", 1}, {"Subrun", 2}, {"Evt", 3}}},
            {"!SPR29", {{"Run", 1}, {"Subrun", 3}, {"Evt", 3}}},
            {"!SPR30", {{"Run", 1}, {"Subrun", 2}, {"Evt", 2}}},
            {"SPR31", {{"Run", 1}, {"Subrun", 3}, {"Evt", 2}, {"reco_particle_ke", 200.0}}},
            {"SPR32", {{"Run", 1}, {"Subrun", 3}, {"Evt", 2}, {"true_particle_ke", 200.0}}},
        };

        // Check if each condition_t entry is present in the rows vector.
        match_conditions(rows, conditions);

        // Read the event data from the TTree in the ROOT file.
        rows = read_event_data("events/test_simlike/test_reco_particles_with_truth_cut");
        
        // Expected results for validation.
        conditions = {
            {"!SPR05", {{"Run", 1}, {"Subrun", 0}, {"Evt", 1}}},
            {"!SPR06", {{"Run", 1}, {"Subrun", 1}, {"Evt", 1}}},
            {"!SPR07", {{"Run", 1}, {"Subrun", 0}, {"Evt", 0}}},
            {"!SPR08", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}}},
            {"!SPR14", {{"Run", 1}, {"Subrun", 2}, {"Evt", 1}}},
            {"!SPR15", {{"Run", 1}, {"Subrun", 3}, {"Evt", 1}}},      
            {"!SPR16", {{"Run", 1}, {"Subrun", 2}, {"Evt", 0}}},
            {"!SPR17", {{"Run", 1}, {"Subrun", 3}, {"Evt", 0}}},
            {"!SPR23", {{"Run", 1}, {"Subrun", 0}, {"Evt", 3}}},
            {"!SPR24", {{"Run", 1}, {"Subrun", 1}, {"Evt", 3}}},
            {"!SPR25", {{"Run", 1}, {"Subrun", 0}, {"Evt", 2}}},
            {"SPR26", {{"Run", 1}, {"Subrun", 1}, {"Evt", 2}, {"reco_particle_ke", 200.0}}},
            {"SPR27", {{"Run", 1}, {"Subrun", 1}, {"Evt", 2}, {"true_particle_ke", kNaN}}},
            {"!SPR33", {{"Run", 1}, {"Subrun", 2}, {"Evt", 3}}},
            {"!SPR34", {{"Run", 1}, {"Subrun", 3}, {"Evt", 3}}},
            {"!SPR35", {{"Run", 1}, {"Subrun", 2}, {"Evt", 2}}},
            {"SPR36", {{"Run", 1}, {"Subrun", 3}, {"Evt", 2}, {"reco_particle_ke", 200.0}}},
            {"SPR37", {{"Run", 1}, {"Subrun", 3}, {"Evt", 2}, {"true_particle_ke", 200.0}}},
        };

        // Check if each condition_t entry is present in the rows vector.
        match_conditions(rows, conditions);

        /**
         * @brief The fifth set of events to validate is the "sim-like" events
         * and the response of the framework when run over them in "true" mode
         * with particle-level variables.
         * @details This set of events effectively tests the framework's
         * behavior when run in a mode where the selection logic is applied
         * to the particles of the interactions, rather than the interactions
         * themselves. Effectively, the framework implements an additional
         * nested loop over the particles (of the same type) for each
         * interaction, and applies some selection logic that is similar to the
         * one used for the interactions.
         * 
         * - SPT00: This represents a true interaction with no valid reco match
         *   under no additional reco cut that would otherwise not pass the
         *   reco-selection. Furthermore, the particle in this interaction has
         *   no reco match and does not pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT01: This represents a true interaction with no valid reco match
         *   under no additional reco cut that would otherwise not pass the
         *   reco-selection. Furthermore, the particle in this interaction has
         *   no reco match and does pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT02: This represents a true interaction with no valid reco match
         *   under no additional reco cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has no
         *   reco match and does not pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT03: This represents a true interaction with no valid reco match
         *   under no additional reco cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has no
         *   reco match and does pass the truth cut. Passes the selection with
         *   a valid truth-var.
         * 
         * - SPT04: This represents a true interaction with no valid reco match
         *   under no additional reco cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has no
         *   reco match and does pass the truth cut. Passes the selection with
         *   a NaN reco-var.
         * 
         * - SPT05: This represents a true interaction with no valid reco match
         *   under an additional truth cut that would otherwise not pass the
         *   reco-selection. Furthermore, the particle in this interaction has
         *   no reco match and does not pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT06: This represents a true interaction with no valid reco match
         *   under an additional truth cut that would otherwise not pass the
         *   reco-selection. Furthermore, the particle in this interaction has
         *   no reco match and does pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT07: This represents a true interaction with no valid reco match
         *   under an additional truth cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has no
         *   reco match and does not pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT08: This represents a true interaction with no valid reco match
         *   under an additional truth cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has no
         *   reco match and does pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT09: This represents a true interaction with no valid reco match
         *   under no additional reco cut that would otherwise not pass the
         *   reco-selection. Furthermore, the particle in this interaction has a
         *   reco match and does not pass the truth cut. Does not pass the
         *   selection.
         *  
         * - SPT10: This represents a true interaction with no valid reco match
         *   under no additional reco cut that would otherwise not pass the
         *   reco-selection. Furthermore, the particle in this interaction has a
         *   reco match and does pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT11: This represents a true interaction with no valid reco match
         *   under no additional reco cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has a
         *   reco match and does not pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT12: This represents a true interaction with no valid reco match
         *   under no additional reco cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has a
         *   reco match and does pass the truth cut. Passes the selection with
         *   a valid truth-var.
         * 
         * - SPT13: This represents a true interaction with no valid reco match
         *   under no additional reco cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has a
         *   reco match and does pass the truth cut. Passes the selection with
         *   a valid reco-var.
         * 
         * - SPT14: This represents a true interaction with no valid reco match
         *   under an additional truth cut that would otherwise not pass the
         *   reco-selection. Furthermore, the particle in this interaction has a
         *   reco match and does not pass the truth cut. Does not pass the
         *   selection.
         *  
         * - SPT15: This represents a true interaction with no valid reco match
         *   under an additional truth cut that would otherwise not pass the
         *   reco-selection. Furthermore, the particle in this interaction has a
         *   reco match and does pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT16: This represents a true interaction with no valid reco match
         *   under an additional truth cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has a
         *   reco match and does not pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT17: This represents a true interaction with no valid reco match
         *   under an additional truth cut that would otherwise pass the reco-
         *   selection. Furthermore, the particle in this interaction has a
         *   reco match and does pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT18: This represents a true interaction with a valid truth match
         *   under no additional reco cut and that does not pass the reco-
         *   selection. Furthermore, the particle in this interaction has no
         *   reco match and does not pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT19: This represents a true interaction with a valid truth match
         *   under no additional reco cut and that does not pass the reco-
         *   selection. Furthermore, the particle in this interaction has no
         *   reco match and does pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT20: This represents a true interaction with a valid truth match
         *   under no additional reco cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has no
         *   reco match and does not pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT21: This represents a true interaction with a valid truth match
         *   under no additional reco cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has no
         *   reco match and does pass the truth cut. Passes the selection with
         *   a valid truth-var.
         * 
         * - SPT22: This represents a true interaction with a valid truth match
         *   under no additional reco cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has no
         *   reco match and does pass the truth cut. Passes the selection with
         *   a NaN reco-var.
         * 
         * - SPT23: This represents a true interaction with a valid truth match
         *   under an additional truth cut and that does not pass the reco-
         *   selection. Furthermore, the particle in this interaction has no
         *   reco match and does not pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT24: This represents a true interaction with a valid truth match
         *   under an additional truth cut and that does not pass the reco-
         *   selection. Furthermore, the particle in this interaction has no
         *   reco match and does pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT25: This represents a true interaction with a valid truth match
         *   under an additional truth cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has no
         *   reco match and does not pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT26: This represents a true interaction with a valid truth match
         *   under an additional truth cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has no
         *   reco match and does pass the truth cut. Passes the selection with
         *   a valid truth-var.
         * 
         * - SPT27: This represents a true interaction with a valid truth match
         *   under an additional truth cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has no
         *   reco match and does pass the truth cut. Passes the selection with
         *   a NaN reco-var.
         * 
         * - SPT28: This represents a true interaction with a valid truth match
         *   under no additional reco cut and that does not pass the reco-
         *   selection. Furthermore, the particle in this interaction has a
         *   reco match and does not pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT29: This represents a true interaction with a valid truth match
         *   under no additional reco cut and that does not pass the reco-
         *   selection. Furthermore, the particle in this interaction has a
         *   reco match and does pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT30: This represents a true interaction with a valid truth match
         *   under no additional reco cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has a
         *   reco match and does not pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT31: This represents a true interaction with a valid truth match
         *   under no additional reco cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has a
         *   reco match and does pass the truth cut. Passes the selection with
         *   a valid truth-var.
         * 
         * - SPT32: This represents a true interaction with a valid truth match
         *   under no additional reco cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has a
         *   reco match and does pass the truth cut. Passes the selection with
         *   a valid reco-var.
         * 
         * - SPT33: This represents a true interaction with a valid truth match
         *   under an additional truth cut and that does not pass the reco-
         *   selection. Furthermore, the particle in this interaction has a
         *   reco match and does not pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT34: This represents a true interaction with a valid truth match
         *   under an additional truth cut and that does not pass the reco-
         *   selection. Furthermore, the particle in this interaction has a
         *   reco match and does pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT35: This represents a true interaction with a valid truth match
         *   under an additional truth cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has a
         *   reco match and does not pass the truth cut. Does not pass the
         *   selection.
         * 
         * - SPT36: This represents a true interaction with a valid truth match
         *   under an additional truth cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has a
         *   reco match and does pass the truth cut. Passes the selection with
         *   a valid truth-var.
         * 
         * - SPT37: This represents a true interaction with a valid truth match
         *   under an additional truth cut that would also pass the reco-only
         *   selection. Furthermore, the particle in this interaction has a
         *   reco match and does pass the truth cut. Passes the selection with
         *   a valid reco-var.
         */
        std::cout << "\n\033[1mSimulation-like events with mode == 'true' and particle-level variables \033[0m" << std::endl;
        
        // Read the event data from the TTree in the ROOT file.
        rows = read_event_data("events/test_simlike/test_truth_particles");

        // Expected results for validation.
        conditions = {
            {"!SPT00", {{"Run", 1}, {"Subrun", 0}, {"Evt", 1}}},
            {"!SPT01", {{"Run", 1}, {"Subrun", 1}, {"Evt", 1}}},
            {"!SPT02", {{"Run", 1}, {"Subrun", 0}, {"Evt", 0}}},
            {"SPT03", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}, {"true_particle_ke", 200.0}}},
            {"SPT04", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}, {"reco_particle_ke", kNaN}}},
            {"!SPT09", {{"Run", 1}, {"Subrun", 2}, {"Evt", 1}}},
            {"!SPT10", {{"Run", 1}, {"Subrun", 3}, {"Evt", 1}}},
            {"!SPT11", {{"Run", 1}, {"Subrun", 2}, {"Evt", 0}}},
            {"SPT12", {{"Run", 1}, {"Subrun", 3}, {"Evt", 0}, {"true_particle_ke", 200.0}}},
            {"SPT13", {{"Run", 1}, {"Subrun", 3}, {"Evt", 0}, {"reco_particle_ke", 200.0}}},
            {"!SPT18", {{"Run", 1}, {"Subrun", 0}, {"Evt", 3}}},
            {"!SPT19", {{"Run", 1}, {"Subrun", 1}, {"Evt", 3}}},
            {"!SPT20", {{"Run", 1}, {"Subrun", 0}, {"Evt", 2}}},
            {"SPT21", {{"Run", 1}, {"Subrun", 1}, {"Evt", 2}, {"true_particle_ke", 200.0}}},
            {"SPT22", {{"Run", 1}, {"Subrun", 1}, {"Evt", 2}, {"reco_particle_ke", kNaN}}},
            {"!SPT28", {{"Run", 1}, {"Subrun", 2}, {"Evt", 3}}},
            {"!SPT29", {{"Run", 1}, {"Subrun", 3}, {"Evt", 3}}},
            {"!SPT30", {{"Run", 1}, {"Subrun", 2}, {"Evt", 2}}},
            {"SPT31", {{"Run", 1}, {"Subrun", 3}, {"Evt", 2}, {"true_particle_ke", 200.0}}},
            {"SPT32", {{"Run", 1}, {"Subrun", 3}, {"Evt", 2}, {"reco_particle_ke", 200.0}}},
        };

        // Check if each condition_t entry is present in the rows vector.
        match_conditions(rows, conditions);

        // Read the event data from the TTree in the ROOT file.
        rows = read_event_data("events/test_simlike/test_truth_particles_with_reco_cut");
        
        // Expected results for validation.
        conditions = {
            {"!SPT05", {{"Run", 1}, {"Subrun", 0}, {"Evt", 1}}},
            {"!SPT06", {{"Run", 1}, {"Subrun", 1}, {"Evt", 1}}},
            {"!SPT07", {{"Run", 1}, {"Subrun", 0}, {"Evt", 0}}},
            {"!SPT08", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}}},
            {"!SPT14", {{"Run", 1}, {"Subrun", 2}, {"Evt", 1}}},
            {"!SPT15", {{"Run", 1}, {"Subrun", 3}, {"Evt", 1}}},      
            {"!SPT16", {{"Run", 1}, {"Subrun", 2}, {"Evt", 0}}},
            {"!SPT17", {{"Run", 1}, {"Subrun", 3}, {"Evt", 0}}},
            {"!SPT23", {{"Run", 1}, {"Subrun", 0}, {"Evt", 3}}},
            {"!SPT24", {{"Run", 1}, {"Subrun", 1}, {"Evt", 3}}},
            {"!SPT25", {{"Run", 1}, {"Subrun", 0}, {"Evt", 2}}},
            {"SPT26", {{"Run", 1}, {"Subrun", 1}, {"Evt", 2}, {"true_particle_ke", 200.0}}},
            {"SPT27", {{"Run", 1}, {"Subrun", 1}, {"Evt", 2}, {"reco_particle_ke", kNaN}}},
            {"!SPT33", {{"Run", 1}, {"Subrun", 2}, {"Evt", 3}}},
            {"!SPT34", {{"Run", 1}, {"Subrun", 3}, {"Evt", 3}}},
            {"!SPT35", {{"Run", 1}, {"Subrun", 2}, {"Evt", 2}}},
            {"SPT36", {{"Run", 1}, {"Subrun", 3}, {"Evt", 2}, {"true_particle_ke", 200.0}}},
            {"SPT37", {{"Run", 1}, {"Subrun", 3}, {"Evt", 2}, {"reco_particle_ke", 200.0}}},
        };

        // Check if each condition_t entry is present in the rows vector.
        match_conditions(rows, conditions);

        /**
         * @brief The sixth set of events to validate is the "data-like" events
         * and the response of the framework when run over them in "reco" mode
         * with particle-level variables.
         * @details This set of events effectively tests the framework's
         * behavior when run in a mode where the selection logic is applied
         * to the particles of the interactions, rather than the interactions
         * themselves. Effectively, the framework implements an additional
         * nested loop over the particles (of the same type) for each
         * interaction, and applies some selection logic that is similar to the
         * one used for the interactions.
         * 
         * - DPR00: This represents a reco interaction that does not pass the
         *   reco-selection and with a particle that does not pass the reco
         *   cut. Does not pass the selection.
         * 
         * - DPR01: This represents a reco interaction that does not pass the
         *   reco-selection and with a particle that does pass the reco cut.
         *   Does not pass the selection.
         * 
         * - DPR02: This represents a reco interaction that passes the reco-
         *   selection and with a particle that does not pass the reco cut.
         *   Does not pass the selection.
         * 
         * - DPR03: This represents a reco interaction that passes the reco-
         *   selection and with a particle that does pass the reco cut. Passes
         *   the selection with a valid reco-var.
         * 
         * - DPR04: This represents a reco interaction that passes the reco-
         *   selection and with a particle that does pass the reco cut. Passes
         *   the selection with a NaN truth-var.
         * 
         * - DPR05: This represents a reco interaction that does not pass the
         *   reco-selection and with a particle that does not pass the reco
         *   cut. There is an additional truth cut at the interaction level.
         *   Does not pass the selection.
         * 
         * - DPR06: This represents a reco interaction that does not pass the
         *   reco-selection and with a particle that does pass the reco cut.
         *   There is an additional truth cut at the interaction level. Does
         *   not pass the selection.
         * 
         * - DPR07: This represents a reco interaction that passes the reco-
         *   selection and with a particle that does not pass the reco cut.
         *   There is an additional truth cut at the interaction level. Does
         *   not pass the selection.
         *  
         * - DPR08: This represents a reco interaction that passes the reco-
         *   selection and with a particle that does pass the reco cut. There
         *   is an additional truth cut at the interaction level. Passes the
         *   selection with a valid reco-var.
         * 
         * - DPR09: This represents a reco interaction that passes the reco-
         *   selection and with a particle that does pass the reco cut. There
         *   is an additional truth cut at the interaction level. Passes the
         *   selection with a NaN truth-var.
         */
        std::cout << "\n\033[1mData-like events with mode == 'reco' and particle-level variables \033[0m" << std::endl;

        // Read the event data from the TTree in the ROOT file.
        rows = read_event_data("events/test_datalike/test_reco_particles");
        
        // Expected results for validation.
        conditions = {
            {"!DPR00", {{"Run", 1}, {"Subrun", 0}, {"Evt", 1}}},
            {"!DPR01", {{"Run", 1}, {"Subrun", 1}, {"Evt", 1}}},
            {"!DPR02", {{"Run", 1}, {"Subrun", 0}, {"Evt", 0}}},
            {"DPR03", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}, {"reco_particle_ke", 200.0}}},
            {"DPR04", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}, {"true_particle_ke", kNaN}}},
            {"!DPR05", {{"Run", 1}, {"Subrun", 0}, {"Evt", 1}}},
            {"!DPR06", {{"Run", 1}, {"Subrun", 1}, {"Evt", 1}}},
            {"!DPR07", {{"Run", 1}, {"Subrun", 0}, {"Evt", 0}}},
            {"DPR08", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}, {"reco_particle_ke", 200.0}}},
            {"DPR09", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}, {"true_particle_ke", kNaN}}},
        };

        // Check if each condition_t entry is present in the rows vector.
        match_conditions(rows, conditions);

        /**
         * @brief The seventh set of events to validate is the "sim-like"
         * events and the response of the framework when run over them in
         * "event" mode (i.e., extraction of event-level variables).
         * @details This set of events effectively tests the framework's
         * behavior when run in a mode where the selection logic is applied to
         * each event as a whole, rather than to the interactions or particles
         * within the interactions.
         * 
         * - SEV00: This represents an event with a trigger time that passes
         *   the cut on the trigger time.
         * 
         * - SEV01: This represents an event with a trigger time that does not
         *   pass the cut on the trigger time.
         */
        std::cout << "\n\033[1mSimulation-like events with mode == 'event' \033[0m" << std::endl;

        // Read the event data from the TTree in the ROOT file.
        rows = read_event_data("events/test_simlike/test_event");

        // Expected results for validation.
        conditions = {
            {"SEV00", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}, {"event_ntrue", 1.0}}},
            {"!SEV01", {{"Run", 1}, {"Subrun", 0}, {"Evt", 0}}},
        };

        // Check if each condition_t entry is present in the rows vector.
        match_conditions(rows, conditions);

        /**
         * @brief The eight set of events to validate is the "sim-like" events
         * and the response of the framework when run over them in "reco" mode
         * with an event-level cut.
         * @details This set of events effectively tests the framework's
         * behavior when run in a mode where the selection logic is applied to
         * reco interactions with an additional event-level cut.
         * 
         * - SER00: This represents a reco interaction that passes the
         *   reco-selection and with a parent event that passes the
         *   event-level cut.
         * 
         * - SER01: This represents a reco interaction that passes the
         *   reco-selection and with a parent event that does not pass the
         *   event-level cut.
         */
        std::cout << "\n\033[1mSimulation-like events with mode == 'reco' and event-level cut \033[0m" << std::endl;

        // Read the event data from the TTree in the ROOT file.
        rows = read_event_data("events/test_simlike/test_reco_with_event_cut");

        // Expected results for validation.
        conditions = {
            {"SER00", {{"Run", 1}, {"Subrun", 1}, {"Evt", 0}, {"reco_vertex_x", -210.0}}},
            {"!SER01", {{"Run", 1}, {"Subrun", 0}, {"Evt", 0}}},
        };

        // Check if each condition_t entry is present in the rows vector.
        match_conditions(rows, conditions);

        // Finished!
        std::cout << "\n\033[1m---        DONE        ---\033[0m" << std::endl;
        f.Close();
    }
    return 0;
}