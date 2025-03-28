/**
 * @file spine_benchmark.cc
 * @brief Main analysis script for benchmarking the SPINE reconstruction.
 * @details This script is used to benchmark the particle-level reconstruction
 * of the SPINE reconstruction. This includes things like PID and primary
 * classification efficiency, energy resolution, and variables that may have
 * some correlation with the reconstruction quality such as angle with respect
 * to wire angles.
 * @author mueller@fnal.gov
 */
#define PLACEHOLDERVALUE std::numeric_limits<double>::quiet_NaN()
#define PRIMARYFUNC pvars::primary_classification
#define PIDFUNC pvars::pid
#define PROTON_BINDING_ENERGY 30.9 // MeV
#define BEAM_IS_NUMI false
#define WRITE_PURITY_TREES false

#include "include/mctruth.h"
#include "include/variables.h"
#include "include/cuts.h"
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
    ana::Analysis analysis("spine_benchmark");

    /**
     * @brief Add a sample to the analysis.
     * @details This adds a sample to the analysis by creating a SpectrumLoader
     * object and adding it to the Analysis class. The SpectrumLoader object
     * represents the sample in the analysis, and is used to load the data from
     * the ROOT file and apply the cuts and variables. The name passed to the
     * AddLoader function is used to create a directory in the output ROOT file
     * to store the results of the analysis.
     */
    ana::SpectrumLoader icarus("/pnfs/icarus/persistent/users/mueller/production/simulation/nominal/flat/input*.flat.root");
    analysis.AddLoader("icarus", &icarus, true);

    auto pid_category = [](const TTYPEP & p) -> double { return 5 * int(pcuts::is_primary(p)) + int(PIDFUNC(p)); };
    auto pid_category_reco = [](const RTYPEP & p) -> double { return 5 * int(pcuts::is_primary(p)) + int(PIDFUNC(p)); };
    auto particle_contained = [](const TTYPEP & p) -> double { return p.is_contained; };
    auto particle_time = [](const TTYPEP & p) -> double { return p.t; };

    std::map<std::string, ana::SpillMultiVar> vars_pid_efficiency;
    vars_pid_efficiency.insert({"pid_category", SpineVar<TTYPEP,TTYPEP,TTYPE>(pid_category, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"semantic_type", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::semantic_type, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"reco_semantic_type", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::semantic_type, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"particle_time", SpineVar<TTYPEP,TTYPEP,TTYPE>(particle_time, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"primary_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::primary_softmax, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"secondary_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::secondary_softmax, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"photon_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::photon_softmax, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"electron_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::electron_softmax, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"muon_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::muon_softmax, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"pion_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::pion_softmax, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"proton_softmax", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::proton_softmax, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"contained", SpineVar<TTYPEP,TTYPEP,TTYPE>(particle_contained, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"ke", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::ke, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"dpT", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::dpT, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"x", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::start_x, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"y", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::start_y, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"z", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::start_z, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"theta", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::polar_angle, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"phi", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::azimuthal_angle, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"theta_xw_horizontal", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::theta_xw_horizontal, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"theta_xw_vertical", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::theta_xw_vertical, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"theta_xw_p60", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::theta_xw_p60, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"theta_xw_p30", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::theta_xw_p30, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"theta_xw_m60", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::theta_xw_m60, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"theta_xw_m30", SpineVar<TTYPEP,TTYPEP,TTYPE>(&pvars::theta_xw_m30, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"reco_pid_category", SpineVar<RTYPEP,TTYPEP,TTYPE>(pid_category_reco, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"reco_ke", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::ke, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"calo_ke", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::calo_ke, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"mcs_ke", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::mcs_ke, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"csda_ke", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::csda_ke, &cuts::no_cut, &cuts::no_cut)});
    vars_pid_efficiency.insert({"iou", SpineVar<RTYPEP,TTYPEP,TTYPE>(&pvars::iou, &cuts::no_cut, &cuts::no_cut)});

    analysis.AddTree("pid_efficiency", vars_pid_efficiency, true);

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
