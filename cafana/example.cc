/**
 * @file example.cc
 * @brief Example of how to run a selection analysis using the
 * spine_anaplot_tools package.
 * @details This file contains an example of how to run a full (if basic)
 * neutrino selection analysis using the spine_anaplot_tools package. This
 * includes the declaration of TTrees with branches implemented by a
 * combination of a selection (cut) function and a variable function.
 * @author mueller@fnal.gov
 */
#define PLACEHOLDERVALUE std::numeric_limits<double>::quiet_NaN()
#define PRIMARYFUNC pvars::custom_primary_classification
#define PIDFUNC pvars::pid
#define PROTON_BINDING_ENERGY 30.9 // MeV
#define BEAM_IS_NUMI false
#define WRITE_PURITY_TREES false

#include "include/mctruth.h"
#include "include/variables.h"
#include "include/muon2024/variables_muon2024.h"
#include "include/cuts.h"
#include "include/muon2024/cuts_muon2024.h"
#include "include/spinevar.h"
#include "include/analysis.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"

int main()
{
    /**
     * @brief Create an instance of the Analysis class.
     * @details This creates an instance of the Analysis class, which is used
     * to run the analysis on the specified samples. The name of the analysis,
     * and therefore the name of the output file, is specified as an argument
     * to the constructor.
     */
    ana::Analysis analysis("nuexample");

    /**
     * @brief Add a sample to the analysis.
     * @details This adds a sample to the analysis by creating a SpectrumLoader
     * object and adding it to the Analysis class. The SpectrumLoader object
     * represents the sample in the analysis, and is used to load the data from
     * the ROOT file and apply the cuts and variables. The name passed to the
     * AddLoader function is used to create a directory in the output ROOT file
     * to store the results of the analysis.
     */
    ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/production/simulation/nominal/flat/input*.flat.root");
    analysis.AddLoader("mc", &mc, true);
    
    /**
     * @brief Add a set of variables for selected interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    #define CUT cuts::fiducial_containment_flash_cut
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_selected_nu;
    vars_selected_nu.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_nu.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_nu.insert({"pdg", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_pdg, &CUT, &TCUT)});
    vars_selected_nu.insert({"cc", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_cc, &CUT, &TCUT)});
    vars_selected_nu.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_selected_nu.insert({"interaction_type", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_type, &CUT, &TCUT)});
    vars_selected_nu.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::muon2024::category, &CUT, &TCUT)});
    vars_selected_nu.insert({"true_energy", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_energy, &CUT, &TCUT)});
    vars_selected_nu.insert({"true_edep", SpineVar<TTYPE,RTYPE>(&vars::visible_energy, &CUT, &TCUT)});
    vars_selected_nu.insert({"reco_edep", SpineVar<RTYPE,RTYPE>(&vars::visible_energy, &CUT, &TCUT)});
    vars_selected_nu.insert({"true_edep_calosub", SpineVar<TTYPE,RTYPE>(&vars::visible_energy_calosub, &CUT, &TCUT)});
    vars_selected_nu.insert({"reco_edep_calosub", SpineVar<RTYPE,RTYPE>(&vars::visible_energy_calosub, &CUT, &TCUT)});
    vars_selected_nu.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &TCUT)});
    vars_selected_nu.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &CUT, &TCUT)});
    vars_selected_nu.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &CUT, &TCUT)});

    analysis.AddTree("selectedNu", vars_selected_nu, false);

    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_selected_cos;
    vars_selected_cos.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_cos.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_cos.insert({"pdg", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_pdg, &CUT, &TCUT)});
    vars_selected_cos.insert({"cc", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_cc, &CUT, &TCUT)});
    vars_selected_cos.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_selected_cos.insert({"interaction_type", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_type, &CUT, &TCUT)});
    vars_selected_cos.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::muon2024::category, &CUT, &TCUT)});
    vars_selected_cos.insert({"true_energy", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_energy, &CUT, &TCUT)});
    vars_selected_cos.insert({"true_edep", SpineVar<TTYPE,RTYPE>(&vars::visible_energy, &CUT, &TCUT)});
    vars_selected_cos.insert({"reco_edep", SpineVar<RTYPE,RTYPE>(&vars::visible_energy, &CUT, &TCUT)});
    vars_selected_cos.insert({"true_edep_calosub", SpineVar<TTYPE,RTYPE>(&vars::visible_energy_calosub, &CUT, &TCUT)});
    vars_selected_cos.insert({"reco_edep_calosub", SpineVar<RTYPE,RTYPE>(&vars::visible_energy_calosub, &CUT, &TCUT)});
    vars_selected_cos.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &TCUT)});
    vars_selected_cos.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &CUT, &TCUT)});
    vars_selected_cos.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &CUT, &TCUT)});
    
    analysis.AddTree("selectedCos", vars_selected_cos, false);
    
    #undef TCUT
    #define TCUT cuts::neutrino
    #undef CUT
    #define CUT cuts::no_cut
    std::map<std::string, ana::SpillMultiVar> vars_purity_nu;
    vars_purity_nu.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_purity_nu.insert({"containment_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::containment_cut), &CUT, &TCUT)});
    vars_purity_nu.insert({"flash_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::flash_cut), &CUT, &TCUT)});
    vars_purity_nu.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &TCUT)});
    vars_purity_nu.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_purity_nu.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_purity_nu.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});

    analysis.AddTree("purityNu", vars_purity_nu, false);

    /**
     * @brief Add a set of variables for signal interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    #define SIGCUT cuts::fiducial_containment_neutrino_cut
    std::map<std::string, ana::SpillMultiVar> vars_signal;
    vars_signal.insert({"nu_id", SpineVar<TTYPE,TTYPE>(&vars::neutrino_id, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"baseline", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_baseline, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"pdg", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_pdg, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"cc", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_cc, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"interaction_mode", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_mode, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"interaction_type", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_type, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"category", SpineVar<TTYPE,TTYPE>(&vars::muon2024::category, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_energy", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_energy, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_edep", SpineVar<TTYPE,TTYPE>(&vars::visible_energy, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_edep_calosub", SpineVar<TTYPE,TTYPE>(&vars::visible_energy_calosub, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"fiducial_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::fiducial_cut), &SIGCUT, &SIGCUT)});
    vars_signal.insert({"containment_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::containment_cut), &SIGCUT, &SIGCUT)});
    vars_signal.insert({"flash_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::flash_cut), &SIGCUT, &SIGCUT)});

    analysis.AddTree("signalNu", vars_signal, true);

    /**
     * @brief Run the analysis.
     * @details This runs the analysis on the samples specified by the
     * SpectrumLoaders and variables added to the Analysis class. It loops over
     * each sample (here only one), applies the cuts and variables to the data,
     * and stores the results in a TFile.
     */
    analysis.Go();

    return 0;
}