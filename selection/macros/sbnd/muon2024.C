/**
 * @file muon2024.C
 * @brief The main analysis macro for the muon2024 analysis on SBND Monte
 * Carlo simulation.
 * @details This macro drives the analysis by configuring the variables, cuts,
 * and samples to be used in the analysis. This is accomplished through the use
 * of the Analysis class, which containerizes the configuration of the analysis
 * and reduces the amount of boilerplate code needed to run the analysis.
 * @author mueller@fnal.gov
*/

/**
 * @brief Block of preprocessor definitions for the analysis.
 * @details This block of preprocessor definitions is used to configure the
 * analysis. The definitions control the behavior of the analysis, such as
 * which beam is used, which cuts are applied, and which trees are created.
 */
#define PLACEHOLDERVALUE std::numeric_limits<double>::quiet_NaN()
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

void muon2024()
{
    /**
     * @brief Create an instance of the Analysis class.
     * @details This creates an instance of the Analysis class, which is used
     * to run the analysis on the specified samples. The name of the analysis,
     * and therefore the name of the output file, is specified as an argument
     * to the constructor.
     */
    ana::Analysis analysis("muon2024_sbnd");

    /**
     * @brief Add samples to the analysis.
     * @details This adds samples to the analysis by creating SpectrumLoader
     * objects and adding them to the Analysis class. The SpectrumLoader object
     * represents the sample in the analysis, and is used to load the data from
     * the ROOT file and apply the cuts and variables. The name passed to the
     * AddLoader function is used to create a directory in the output ROOT file
     * to store the results of the analysis.
     */
    ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/sbnd/larcv_sbnd_bnb_cosmics_spine.flat.root");
    analysis.AddLoader("mc", &mc, true);

    ana::SpectrumLoader intime("/pnfs/icarus/persistent/users/mueller/sbnd/larcv_sbnd_intime_spine.flat.root");
    analysis.AddLoader("intime", &intime, true);
    
    /**
     * @brief Add a set of variables for selected interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    #define CUT cuts::muon2024::all_1muNp_cut
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_selected_nu;
    vars_selected_nu.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_nu.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_nu.insert({"cc", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_cc, &CUT, &TCUT)});
    vars_selected_nu.insert({"pdg", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_pdg, &CUT, &TCUT)});
    vars_selected_nu.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_selected_nu.insert({"interaction_type", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_type, &CUT, &TCUT)});
    vars_selected_nu.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::muon2024::category, &CUT, &TCUT)});
    vars_selected_nu.insert({"true_energy", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_energy, &CUT, &TCUT)});
    vars_selected_nu.insert({"true_edep", SpineVar<TTYPE,RTYPE>(&vars::visible_energy, &CUT, &TCUT)});
    vars_selected_nu.insert({"reco_edep", SpineVar<RTYPE,RTYPE>(&vars::visible_energy, &CUT, &TCUT)});
    vars_selected_nu.insert({"true_edep_calosub", SpineVar<TTYPE,RTYPE>(&vars::visible_energy_calosub, &CUT, &TCUT)});
    vars_selected_nu.insert({"reco_edep_calosub", SpineVar<RTYPE,RTYPE>(&vars::visible_energy_calosub, &CUT, &TCUT)});
    vars_selected_nu.insert({"true_muon_x", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::end_x, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"reco_muon_x", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_x, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"true_muon_y", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::end_y, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"reco_muon_y", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_y, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"true_muon_z", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::end_z, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"reco_muon_z", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_z, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"true_proton_x", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::end_x, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"reco_proton_x", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_x, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"true_proton_y", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::end_y, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"reco_proton_y", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_y, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"true_proton_z", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::end_z, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"reco_proton_z", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_z, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"true_tmuon", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::ke, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"reco_tmuon", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::ke, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"true_lmuon", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::length, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"reco_lmuon", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::length, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"true_tproton", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::ke, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"reco_tproton", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::ke, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"true_lproton", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::length, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"reco_lproton", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::length, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"true_ptmuon", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::dpT, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"reco_ptmuon", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::dpT, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"true_ptproton", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::dpT, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"reco_ptproton", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::dpT, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"true_theta_mu", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::polar_angle, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"reco_theta_mu", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::polar_angle, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"true_phi_mu", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::azimuthal_angle, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"reco_phi_mu", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::azimuthal_angle, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"true_opening_angle", SpineVar<TTYPE,RTYPE>(&vars::muon2024::opening_angle, &CUT, &TCUT)});
    vars_selected_nu.insert({"reco_opening_angle", SpineVar<RTYPE,RTYPE>(&vars::muon2024::opening_angle, &CUT, &TCUT)});
    vars_selected_nu.insert({"true_dpT", SpineVar<TTYPE,RTYPE>(&vars::dpT, &CUT, &TCUT)});
    vars_selected_nu.insert({"reco_dpT", SpineVar<RTYPE,RTYPE>(&vars::dpT, &CUT, &TCUT)});
    vars_selected_nu.insert({"true_dpT_lp", SpineVar<TTYPE,RTYPE>(&vars::dpT_lp, &CUT, &TCUT)});
    vars_selected_nu.insert({"reco_dpT_lp", SpineVar<RTYPE,RTYPE>(&vars::dpT_lp, &CUT, &TCUT)});
    vars_selected_nu.insert({"true_dphiT", SpineVar<TTYPE,RTYPE>(&vars::phiT, &CUT, &TCUT)});
    vars_selected_nu.insert({"reco_dphiT", SpineVar<RTYPE,RTYPE>(&vars::phiT, &CUT, &TCUT)});
    vars_selected_nu.insert({"true_dalphaT", SpineVar<TTYPE,RTYPE>(&vars::alphaT, &CUT, &TCUT)});
    vars_selected_nu.insert({"reco_dalphaT", SpineVar<RTYPE,RTYPE>(&vars::alphaT, &CUT, &TCUT)});
    vars_selected_nu.insert({"true_dpL", SpineVar<TTYPE,RTYPE>(&vars::dpL, &CUT, &TCUT)});
    vars_selected_nu.insert({"reco_dpL", SpineVar<RTYPE,RTYPE>(&vars::dpL, &CUT, &TCUT)});
    vars_selected_nu.insert({"true_dpL_lp", SpineVar<TTYPE,RTYPE>(&vars::dpL_lp, &CUT, &TCUT)});
    vars_selected_nu.insert({"reco_dpL_lp", SpineVar<RTYPE,RTYPE>(&vars::dpL_lp, &CUT, &TCUT)});
    vars_selected_nu.insert({"true_pn", SpineVar<TTYPE,RTYPE>(&vars::pn, &CUT, &TCUT)});
    vars_selected_nu.insert({"reco_pn", SpineVar<RTYPE,RTYPE>(&vars::pn, &CUT, &TCUT)});
    vars_selected_nu.insert({"true_pn_lp", SpineVar<TTYPE,RTYPE>(&vars::pn_lp, &CUT, &TCUT)});
    vars_selected_nu.insert({"reco_pn_lp", SpineVar<RTYPE,RTYPE>(&vars::pn_lp, &CUT, &TCUT)});
    vars_selected_nu.insert({"true_vertex_x", SpineVar<TTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_nu.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_nu.insert({"true_vertex_y", SpineVar<TTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_nu.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_nu.insert({"true_vertex_z", SpineVar<TTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    vars_selected_nu.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    vars_selected_nu.insert({"muon_primary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::primary_softmax, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"muon_secondary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::secondary_softmax, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"muon_muon_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::muon_softmax, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"muon_pion_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::pion_softmax, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"muon_proton_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::proton_softmax, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"muon_mip_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::mip_softmax, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"muon_hadron_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::hadron_softmax, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_nu.insert({"proton_primary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::primary_softmax, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"proton_secondary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::secondary_softmax, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"proton_muon_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::muon_softmax, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"proton_pion_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::pion_softmax, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"proton_proton_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::proton_softmax, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"proton_mip_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::mip_softmax, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"proton_hadron_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::hadron_softmax, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_nu.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &TCUT)});
    vars_selected_nu.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &CUT, &TCUT)});
    vars_selected_nu.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &CUT, &TCUT)});

    analysis.AddTree("selectedNu", vars_selected_nu, false);

    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_selected_cos;
    vars_selected_cos.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_selected_cos.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_selected_cos.insert({"cc", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_cc, &CUT, &TCUT)});
    vars_selected_cos.insert({"pdg", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_pdg, &CUT, &TCUT)});
    vars_selected_cos.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_selected_cos.insert({"interaction_type", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_type, &CUT, &TCUT)});
    vars_selected_cos.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::muon2024::category, &CUT, &TCUT)});
    vars_selected_cos.insert({"true_energy", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_energy, &CUT, &TCUT)});
    vars_selected_cos.insert({"true_edep", SpineVar<TTYPE,RTYPE>(&vars::visible_energy, &CUT, &TCUT)});
    vars_selected_cos.insert({"reco_edep", SpineVar<RTYPE,RTYPE>(&vars::visible_energy, &CUT, &TCUT)});
    vars_selected_cos.insert({"true_edep_calosub", SpineVar<TTYPE,RTYPE>(&vars::visible_energy_calosub, &CUT, &TCUT)});
    vars_selected_cos.insert({"reco_edep_calosub", SpineVar<RTYPE,RTYPE>(&vars::visible_energy_calosub, &CUT, &TCUT)});
    vars_selected_cos.insert({"true_muon_x", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::end_x, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"reco_muon_x", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_x, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"true_muon_y", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::end_y, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"reco_muon_y", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_y, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"true_muon_z", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::end_z, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"reco_muon_z", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_z, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"true_proton_x", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::end_x, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"reco_proton_x", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_x, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"true_proton_y", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::end_y, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"reco_proton_y", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_y, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"true_proton_z", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::end_z, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"reco_proton_z", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_z, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"true_tmuon", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::ke, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"reco_tmuon", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::ke, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"true_lmuon", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::length, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"reco_lmuon", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::length, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"true_tproton", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::ke, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"reco_tproton", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::ke, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"true_lproton", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::length, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"reco_lproton", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::length, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"true_ptmuon", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::dpT, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"reco_ptmuon", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::dpT, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"true_ptproton", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::dpT, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"reco_ptproton", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::dpT, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"true_theta_mu", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::polar_angle, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"reco_theta_mu", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::polar_angle, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"true_phi_mu", SpineVar<TTYPEP,RTYPE,TTYPE>(&pvars::azimuthal_angle, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"reco_phi_mu", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::azimuthal_angle, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"true_opening_angle", SpineVar<TTYPE,RTYPE>(&vars::muon2024::opening_angle, &CUT, &TCUT)});
    vars_selected_cos.insert({"reco_opening_angle", SpineVar<RTYPE,RTYPE>(&vars::muon2024::opening_angle, &CUT, &TCUT)});
    vars_selected_cos.insert({"true_dpT", SpineVar<TTYPE,RTYPE>(&vars::dpT, &CUT, &TCUT)});
    vars_selected_cos.insert({"reco_dpT", SpineVar<RTYPE,RTYPE>(&vars::dpT, &CUT, &TCUT)});
    vars_selected_cos.insert({"true_dpT_lp", SpineVar<TTYPE,RTYPE>(&vars::dpT_lp, &CUT, &TCUT)});
    vars_selected_cos.insert({"reco_dpT_lp", SpineVar<RTYPE,RTYPE>(&vars::dpT_lp, &CUT, &TCUT)});
    vars_selected_cos.insert({"true_dphiT", SpineVar<TTYPE,RTYPE>(&vars::phiT, &CUT, &TCUT)});
    vars_selected_cos.insert({"reco_dphiT", SpineVar<RTYPE,RTYPE>(&vars::phiT, &CUT, &TCUT)});
    vars_selected_cos.insert({"true_dalphaT", SpineVar<TTYPE,RTYPE>(&vars::alphaT, &CUT, &TCUT)});
    vars_selected_cos.insert({"reco_dalphaT", SpineVar<RTYPE,RTYPE>(&vars::alphaT, &CUT, &TCUT)});
    vars_selected_cos.insert({"true_dpL", SpineVar<TTYPE,RTYPE>(&vars::dpL, &CUT, &TCUT)});
    vars_selected_cos.insert({"reco_dpL", SpineVar<RTYPE,RTYPE>(&vars::dpL, &CUT, &TCUT)});
    vars_selected_cos.insert({"true_dpL_lp", SpineVar<TTYPE,RTYPE>(&vars::dpL_lp, &CUT, &TCUT)});
    vars_selected_cos.insert({"reco_dpL_lp", SpineVar<RTYPE,RTYPE>(&vars::dpL_lp, &CUT, &TCUT)});
    vars_selected_cos.insert({"true_pn", SpineVar<TTYPE,RTYPE>(&vars::pn, &CUT, &TCUT)});
    vars_selected_cos.insert({"reco_pn", SpineVar<RTYPE,RTYPE>(&vars::pn, &CUT, &TCUT)});
    vars_selected_cos.insert({"true_pn_lp", SpineVar<TTYPE,RTYPE>(&vars::pn_lp, &CUT, &TCUT)});
    vars_selected_cos.insert({"reco_pn_lp", SpineVar<RTYPE,RTYPE>(&vars::pn_lp, &CUT, &TCUT)});
    vars_selected_cos.insert({"true_vertex_x", SpineVar<TTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_cos.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_selected_cos.insert({"true_vertex_y", SpineVar<TTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_cos.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_selected_cos.insert({"true_vertex_z", SpineVar<TTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    vars_selected_cos.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    vars_selected_cos.insert({"muon_primary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::primary_softmax, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"muon_secondary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::secondary_softmax, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"muon_muon_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::muon_softmax, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"muon_pion_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::pion_softmax, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"muon_proton_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::proton_softmax, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"muon_mip_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::mip_softmax, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"muon_hadron_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::hadron_softmax, &CUT, &TCUT, &utilities::leading_muon_index)});
    vars_selected_cos.insert({"proton_primary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::primary_softmax, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"proton_secondary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::secondary_softmax, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"proton_muon_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::muon_softmax, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"proton_pion_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::pion_softmax, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"proton_proton_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::proton_softmax, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"proton_mip_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::mip_softmax, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"proton_hadron_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::hadron_softmax, &CUT, &TCUT, &utilities::leading_proton_index)});
    vars_selected_cos.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &CUT, &TCUT)});
    vars_selected_cos.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &CUT, &TCUT)});
    vars_selected_cos.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &CUT, &TCUT)});

    analysis.AddTree("selectedCos", vars_selected_cos, false);

    #undef TCUT
    #define TCUT cuts::neutrino
    std::map<std::string, ana::SpillMultiVar> vars_purity_nu;
    vars_purity_nu.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::muon2024::category, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"reco_edep", SpineVar<RTYPE,RTYPE>(&vars::visible_energy, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"reco_edep_calosub", SpineVar<RTYPE,RTYPE>(&vars::visible_energy_calosub, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"reco_muon_x", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_x, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_nu.insert({"reco_muon_y", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_y, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_nu.insert({"reco_muon_z", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_z, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_nu.insert({"reco_proton_x", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_x, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_nu.insert({"reco_proton_y", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_y, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_nu.insert({"reco_proton_z", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_z, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_nu.insert({"reco_tmuon", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::ke, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_nu.insert({"reco_lmuon", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::length, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_nu.insert({"reco_tproton", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::ke, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_nu.insert({"reco_lproton", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::length, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_nu.insert({"reco_ptmuon", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::dpT, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_nu.insert({"reco_ptproton", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::dpT, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_nu.insert({"reco_theta_mu", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::polar_angle, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_nu.insert({"reco_phi_mu", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::azimuthal_angle, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_nu.insert({"reco_opening_angle", SpineVar<RTYPE,RTYPE>(&vars::muon2024::opening_angle, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"reco_dpT", SpineVar<RTYPE,RTYPE>(&vars::dpT, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"reco_dpT_lp", SpineVar<RTYPE,RTYPE>(&vars::dpT_lp, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"reco_dphiT", SpineVar<RTYPE,RTYPE>(&vars::phiT, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"reco_dalphaT", SpineVar<RTYPE,RTYPE>(&vars::alphaT, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"reco_dpL", SpineVar<RTYPE,RTYPE>(&vars::dpL, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"reco_dpL_lp", SpineVar<RTYPE,RTYPE>(&vars::dpL_lp, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"reco_pn", SpineVar<RTYPE,RTYPE>(&vars::pn, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"reco_pn_lp", SpineVar<RTYPE,RTYPE>(&vars::pn_lp, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"muon_primary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::primary_softmax, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_nu.insert({"muon_secondary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::secondary_softmax, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_nu.insert({"muon_muon_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::muon_softmax, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_nu.insert({"muon_pion_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::pion_softmax, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_nu.insert({"muon_proton_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::proton_softmax, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_nu.insert({"muon_mip_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::mip_softmax, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_nu.insert({"muon_hadron_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::hadron_softmax, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_nu.insert({"proton_primary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::primary_softmax, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_nu.insert({"proton_secondary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::secondary_softmax, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_nu.insert({"proton_muon_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::muon_softmax, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_nu.insert({"proton_pion_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::pion_softmax, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_nu.insert({"proton_proton_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::proton_softmax, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_nu.insert({"proton_mip_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::mip_softmax, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_nu.insert({"proton_hadron_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::hadron_softmax, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_nu.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"fiducial_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::fiducial_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"containment_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::containment_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"flash_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::flash_cut), &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"has_no_charged_pions", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::no_charged_pions), &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"has_no_showers", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::no_showers), &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"has_single_muon", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::has_single_muon), &cuts::no_cut, &TCUT)});
    vars_purity_nu.insert({"has_multiple_protons", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::has_nonzero_protons), &cuts::no_cut, &TCUT)});

    if constexpr(WRITE_PURITY_TREES)
        analysis.AddTree("purityNu", vars_purity_nu, false);

    #undef TCUT
    #define TCUT cuts::cosmic
    std::map<std::string, ana::SpillMultiVar> vars_purity_cos;
    vars_purity_cos.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::muon2024::category, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"reco_edep", SpineVar<RTYPE,RTYPE>(&vars::visible_energy, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"reco_edep_calosub", SpineVar<RTYPE,RTYPE>(&vars::visible_energy_calosub, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"reco_muon_x", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_x, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_cos.insert({"reco_muon_y", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_y, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_cos.insert({"reco_muon_z", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_z, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_cos.insert({"reco_proton_x", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_x, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_cos.insert({"reco_proton_y", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_y, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_cos.insert({"reco_proton_z", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::end_z, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_cos.insert({"reco_tmuon", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::ke, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_cos.insert({"reco_lmuon", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::length, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_cos.insert({"reco_tproton", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::ke, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_cos.insert({"reco_lproton", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::length, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_cos.insert({"reco_ptmuon", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::dpT, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_cos.insert({"reco_ptproton", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::dpT, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_cos.insert({"reco_theta_mu", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::polar_angle, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_cos.insert({"reco_phi_mu", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::azimuthal_angle, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_cos.insert({"reco_opening_angle", SpineVar<RTYPE,RTYPE>(&vars::muon2024::opening_angle, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"reco_dpT", SpineVar<RTYPE,RTYPE>(&vars::dpT, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"reco_dpT_lp", SpineVar<RTYPE,RTYPE>(&vars::dpT_lp, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"reco_dphiT", SpineVar<RTYPE,RTYPE>(&vars::phiT, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"reco_dalphaT", SpineVar<RTYPE,RTYPE>(&vars::alphaT, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"reco_dpL", SpineVar<RTYPE,RTYPE>(&vars::dpL, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"reco_dpL_lp", SpineVar<RTYPE,RTYPE>(&vars::dpL_lp, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"reco_pn", SpineVar<RTYPE,RTYPE>(&vars::pn, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"reco_pn_lp", SpineVar<RTYPE,RTYPE>(&vars::pn_lp, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"muon_primary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::primary_softmax, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_cos.insert({"muon_secondary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::secondary_softmax, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_cos.insert({"muon_muon_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::muon_softmax, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_cos.insert({"muon_pion_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::pion_softmax, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_cos.insert({"muon_proton_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::proton_softmax, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_cos.insert({"muon_mip_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::mip_softmax, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_cos.insert({"muon_hadron_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::hadron_softmax, &cuts::no_cut, &TCUT, &utilities::leading_muon_index)});
    vars_purity_cos.insert({"proton_primary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::primary_softmax, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_cos.insert({"proton_secondary_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::secondary_softmax, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_cos.insert({"proton_muon_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::muon_softmax, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_cos.insert({"proton_pion_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::pion_softmax, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_cos.insert({"proton_proton_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::proton_softmax, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_cos.insert({"proton_mip_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::mip_softmax, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_cos.insert({"proton_hadron_softmax", SpineVar<RTYPEP,RTYPE,RTYPE>(&pvars::hadron_softmax, &cuts::no_cut, &TCUT, &utilities::leading_proton_index)});
    vars_purity_cos.insert({"flash_time", SpineVar<RTYPE,RTYPE>(&vars::flash_time, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"flash_total", SpineVar<RTYPE,RTYPE>(&vars::flash_total_pe, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"flash_hypothesis", SpineVar<RTYPE,RTYPE>(&vars::flash_hypothesis, &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"fiducial_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::fiducial_cut), &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"containment_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::containment_cut), &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"flash_cut", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::flash_cut), &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"has_no_charged_pions", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::no_charged_pions), &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"has_no_showers", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::no_showers), &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"has_single_muon", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::has_single_muon), &cuts::no_cut, &TCUT)});
    vars_purity_cos.insert({"has_multiple_protons", SpineVar<RTYPE,RTYPE>(WRAP_BOOL(cuts::has_nonzero_protons), &cuts::no_cut, &TCUT)});

    if constexpr(WRITE_PURITY_TREES)
        analysis.AddTree("purityCos", vars_purity_cos, false);

    /**
     * @brief Add a set of variables for signal interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    #define SIGCUT cuts::muon2024::signal_1muNp
    std::map<std::string, ana::SpillMultiVar> vars_signal;
    vars_signal.insert({"nu_id", SpineVar<TTYPE,TTYPE>(&vars::neutrino_id, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"baseline", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_baseline, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"cc", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_cc, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"pdg", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_pdg, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"interaction_mode", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_mode, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"interaction_type", SpineVar<MCTRUTH,TTYPE>(&mctruth::interaction_type, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"category", SpineVar<TTYPE,TTYPE>(&vars::muon2024::category, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_energy", SpineVar<MCTRUTH,TTYPE>(&mctruth::true_neutrino_energy, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_edep", SpineVar<TTYPE,TTYPE>(&vars::visible_energy, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_edep_calosub", SpineVar<TTYPE,TTYPE>(&vars::visible_energy_calosub, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_muon_x", SpineVar<TTYPEP,TTYPE,TTYPE>(&pvars::end_x, &SIGCUT, &SIGCUT, &utilities::leading_muon_index)});
    vars_signal.insert({"true_muon_y", SpineVar<TTYPEP,TTYPE,TTYPE>(&pvars::end_y, &SIGCUT, &SIGCUT, &utilities::leading_muon_index)});
    vars_signal.insert({"true_muon_z", SpineVar<TTYPEP,TTYPE,TTYPE>(&pvars::end_z, &SIGCUT, &SIGCUT, &utilities::leading_muon_index)});
    vars_signal.insert({"true_proton_x", SpineVar<TTYPEP,TTYPE,TTYPE>(&pvars::end_x, &SIGCUT, &SIGCUT, &utilities::leading_proton_index)});
    vars_signal.insert({"true_proton_y", SpineVar<TTYPEP,TTYPE,TTYPE>(&pvars::end_y, &SIGCUT, &SIGCUT, &utilities::leading_proton_index)});
    vars_signal.insert({"true_proton_z", SpineVar<TTYPEP,TTYPE,TTYPE>(&pvars::end_z, &SIGCUT, &SIGCUT, &utilities::leading_proton_index)});
    vars_signal.insert({"true_tmuon", SpineVar<TTYPEP,TTYPE,TTYPE>(&pvars::ke, &SIGCUT, &SIGCUT, &utilities::leading_muon_index)});
    vars_signal.insert({"true_lmuon", SpineVar<TTYPEP,TTYPE,TTYPE>(&pvars::length, &SIGCUT, &SIGCUT, &utilities::leading_muon_index)});
    vars_signal.insert({"true_tproton", SpineVar<TTYPEP,TTYPE,TTYPE>(&pvars::ke, &SIGCUT, &SIGCUT, &utilities::leading_proton_index)});
    vars_signal.insert({"true_lproton", SpineVar<TTYPEP,TTYPE,TTYPE>(&pvars::length, &SIGCUT, &SIGCUT, &utilities::leading_proton_index)});
    vars_signal.insert({"true_ptmuon", SpineVar<TTYPEP,TTYPE,TTYPE>(&pvars::dpT, &SIGCUT, &SIGCUT, &utilities::leading_muon_index)});
    vars_signal.insert({"true_ptproton", SpineVar<TTYPEP,TTYPE,TTYPE>(&pvars::dpT, &SIGCUT, &SIGCUT, &utilities::leading_proton_index)});
    vars_signal.insert({"true_theta_mu", SpineVar<TTYPEP,TTYPE,TTYPE>(&pvars::polar_angle, &SIGCUT, &SIGCUT, &utilities::leading_muon_index)});
    vars_signal.insert({"true_phi_mu", SpineVar<TTYPEP,TTYPE,TTYPE>(&pvars::azimuthal_angle, &SIGCUT, &SIGCUT, &utilities::leading_muon_index)});
    vars_signal.insert({"true_opening_angle", SpineVar<TTYPE,TTYPE>(&vars::muon2024::opening_angle, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_dpT", SpineVar<TTYPE,TTYPE>(&vars::dpT, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_dpT_lp", SpineVar<TTYPE,TTYPE>(&vars::dpT_lp, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_dphiT", SpineVar<TTYPE,TTYPE>(&vars::phiT, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_dalphaT", SpineVar<TTYPE,TTYPE>(&vars::alphaT, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_dpL", SpineVar<TTYPE,TTYPE>(&vars::dpL, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_dpL_lp", SpineVar<TTYPE,TTYPE>(&vars::dpL_lp, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_pn", SpineVar<TTYPE,TTYPE>(&vars::pn, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_pn_lp", SpineVar<TTYPE,TTYPE>(&vars::pn_lp, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_vertex_x", SpineVar<TTYPE,TTYPE>(&vars::vertex_x, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_vertex_y", SpineVar<TTYPE,TTYPE>(&vars::vertex_y, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"true_vertex_z", SpineVar<TTYPE,TTYPE>(&vars::vertex_z, &SIGCUT, &SIGCUT)});
    vars_signal.insert({"muon_primary_softmax", SpineVar<RTYPEP,TTYPE,TTYPE>(&pvars::primary_softmax, &SIGCUT, &SIGCUT, &utilities::leading_muon_index)});
    vars_signal.insert({"muon_secondary_softmax", SpineVar<RTYPEP,TTYPE,TTYPE>(&pvars::secondary_softmax, &SIGCUT, &SIGCUT, &utilities::leading_muon_index)});
    vars_signal.insert({"muon_muon_softmax", SpineVar<RTYPEP,TTYPE,TTYPE>(&pvars::muon_softmax, &SIGCUT, &SIGCUT, &utilities::leading_muon_index)});
    vars_signal.insert({"muon_pion_softmax", SpineVar<RTYPEP,TTYPE,TTYPE>(&pvars::pion_softmax, &SIGCUT, &SIGCUT, &utilities::leading_muon_index)});
    vars_signal.insert({"muon_proton_softmax", SpineVar<RTYPEP,TTYPE,TTYPE>(&pvars::proton_softmax, &SIGCUT, &SIGCUT, &utilities::leading_muon_index)});
    vars_signal.insert({"muon_mip_softmax", SpineVar<RTYPEP,TTYPE,TTYPE>(&pvars::mip_softmax, &SIGCUT, &SIGCUT, &utilities::leading_muon_index)});
    vars_signal.insert({"muon_hadron_softmax", SpineVar<RTYPEP,TTYPE,TTYPE>(&pvars::hadron_softmax, &SIGCUT, &SIGCUT, &utilities::leading_muon_index)});
    vars_signal.insert({"proton_primary_softmax", SpineVar<RTYPEP,TTYPE,TTYPE>(&pvars::primary_softmax, &SIGCUT, &SIGCUT, &utilities::leading_proton_index)});
    vars_signal.insert({"proton_secondary_softmax", SpineVar<RTYPEP,TTYPE,TTYPE>(&pvars::secondary_softmax, &SIGCUT, &SIGCUT, &utilities::leading_proton_index)});
    vars_signal.insert({"proton_muon_softmax", SpineVar<RTYPEP,TTYPE,TTYPE>(&pvars::muon_softmax, &SIGCUT, &SIGCUT, &utilities::leading_proton_index)});
    vars_signal.insert({"proton_pion_softmax", SpineVar<RTYPEP,TTYPE,TTYPE>(&pvars::pion_softmax, &SIGCUT, &SIGCUT, &utilities::leading_proton_index)});
    vars_signal.insert({"proton_proton_softmax", SpineVar<RTYPEP,TTYPE,TTYPE>(&pvars::proton_softmax, &SIGCUT, &SIGCUT, &utilities::leading_proton_index)});
    vars_signal.insert({"proton_mip_softmax", SpineVar<RTYPEP,TTYPE,TTYPE>(&pvars::mip_softmax, &SIGCUT, &SIGCUT, &utilities::leading_proton_index)});
    vars_signal.insert({"proton_hadron_softmax", SpineVar<RTYPEP,TTYPE,TTYPE>(&pvars::hadron_softmax, &SIGCUT, &SIGCUT, &utilities::leading_proton_index)});
    vars_signal.insert({"fiducial_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::fiducial_cut), &SIGCUT, &SIGCUT)});
    vars_signal.insert({"containment_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::containment_cut), &SIGCUT, &SIGCUT)});
    vars_signal.insert({"flash_cut", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::flash_cut), &SIGCUT, &SIGCUT)});
    vars_signal.insert({"has_no_charged_pions", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::no_charged_pions), &SIGCUT, &SIGCUT)});
    vars_signal.insert({"has_no_showers", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::no_showers), &SIGCUT, &SIGCUT)});
    vars_signal.insert({"has_single_muon", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::has_single_muon), &SIGCUT, &SIGCUT)});
    vars_signal.insert({"has_multiple_protons", SpineVar<RTYPE,TTYPE>(WRAP_BOOL(cuts::has_nonzero_protons), &SIGCUT, &SIGCUT)});
    
    analysis.AddTree("signal", vars_signal, true);

    std::map<std::string, ana::SpillMultiVar> vars_pid_muon;
    auto primary_muon = [](const TTYPEP & p) { return PIDFUNC(p) == 2 && p.is_primary; };
    auto primary_muon_reco = [](const RTYPEP & p) -> double { return PIDFUNC(p) == 2 && p.is_primary; };
    auto pid = [](const RTYPEP & p) -> double { return PIDFUNC(p); };
    vars_pid_muon.insert({"muon_primary_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::primary_softmax, primary_muon, &cuts::no_cut)});
    vars_pid_muon.insert({"muon_secondary_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::secondary_softmax, primary_muon, &cuts::no_cut)});
    vars_pid_muon.insert({"muon_muon_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::muon_softmax, primary_muon, &cuts::no_cut)});
    vars_pid_muon.insert({"muon_pion_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::pion_softmax, primary_muon, &cuts::no_cut)});
    vars_pid_muon.insert({"muon_proton_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::proton_softmax, primary_muon, &cuts::no_cut)});
    vars_pid_muon.insert({"muon_mip_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::mip_softmax, primary_muon, &cuts::no_cut)});
    vars_pid_muon.insert({"muon_hadron_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::hadron_softmax, primary_muon, &cuts::no_cut)});
    vars_pid_muon.insert({"true_tmuon", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::ke, primary_muon, &cuts::no_cut)});
    vars_pid_muon.insert({"true_ptmuon", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::dpT, primary_muon, &cuts::no_cut)});
    vars_pid_muon.insert({"true_theta_mu", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::polar_angle, primary_muon, &cuts::no_cut)});
    vars_pid_muon.insert({"true_phi_mu", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::azimuthal_angle, primary_muon, &cuts::no_cut)});
    vars_pid_muon.insert({"pid", SpineVar<RTYPEP,TTYPEP,TTYPE>(pid, primary_muon, &cuts::no_cut)});
    vars_pid_muon.insert({"primary", SpineVar<RTYPEP,TTYPEP,TTYPE>(WRAP_BOOL(pcuts::is_primary), primary_muon, &cuts::no_cut)});
    vars_pid_muon.insert({"is_primary_muon", SpineVar<RTYPEP,TTYPEP,TTYPE>(primary_muon_reco, primary_muon, &cuts::no_cut)});

    analysis.AddTree("pidMuon", vars_pid_muon, true);

    std::map<std::string, ana::SpillMultiVar> vars_pid_proton;
    auto primary_proton = [](const TTYPEP & p) { return PIDFUNC(p) == 4 && p.is_primary; };
    auto primary_proton_reco = [](const RTYPEP & p) -> double { return PIDFUNC(p) == 4 && p.is_primary; };
    vars_pid_proton.insert({"proton_primary_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::primary_softmax, primary_proton, &cuts::no_cut)});
    vars_pid_proton.insert({"proton_secondary_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::secondary_softmax, primary_proton, &cuts::no_cut)});
    vars_pid_proton.insert({"proton_muon_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::muon_softmax, primary_proton, &cuts::no_cut)});
    vars_pid_proton.insert({"proton_pion_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::pion_softmax, primary_proton, &cuts::no_cut)});
    vars_pid_proton.insert({"proton_proton_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::proton_softmax, primary_proton, &cuts::no_cut)});
    vars_pid_proton.insert({"proton_mip_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::mip_softmax, primary_proton, &cuts::no_cut)});
    vars_pid_proton.insert({"proton_hadron_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::hadron_softmax, primary_proton, &cuts::no_cut)});
    vars_pid_proton.insert({"true_tproton", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::ke, primary_proton, &cuts::no_cut)});
    vars_pid_proton.insert({"true_ptproton", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::dpT, primary_proton, &cuts::no_cut)});
    vars_pid_proton.insert({"true_theta_proton", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::polar_angle, primary_proton, &cuts::no_cut)});
    vars_pid_proton.insert({"true_phi_proton", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::azimuthal_angle, primary_proton, &cuts::no_cut)});
    vars_pid_proton.insert({"pid", SpineVar<RTYPEP,TTYPEP,TTYPE>(pid, primary_proton, &cuts::no_cut)});
    vars_pid_proton.insert({"primary", SpineVar<RTYPEP,TTYPEP,TTYPE>(WRAP_BOOL(pcuts::is_primary), primary_proton, &cuts::no_cut)});
    vars_pid_proton.insert({"is_primary_proton", SpineVar<RTYPEP,TTYPEP,TTYPE>(primary_proton_reco, primary_proton, &cuts::no_cut)});

    analysis.AddTree("pidProton", vars_pid_proton, true);

    /**
     * @brief Run the analysis.
     * @details This runs the analysis on the samples specified by the
     * SpectrumLoaders and variables added to the Analysis class. It loops over
     * each sample (here only one), applies the cuts and variables to the data,
     * and stores the results in a TFile.
     */
    analysis.Go();
}