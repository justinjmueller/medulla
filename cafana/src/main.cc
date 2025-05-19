/**
 * @file main.cc
 * @brief Main file for the SPINE analysis framework.
 * @details This file contains the main function for the SPINE analysis
 * framework. The main function is responsible for loading the configuration
 * file, initializing the analysis framework, and running the analysis.
 * @author mueller@fnal.gov
 */
#define PLACEHOLDERVALUE std::numeric_limits<double>::quiet_NaN()
#define PRIMARYFUNC pvars::custom_primary_classification
#define PIDFUNC pvars::pid
#define PROTON_BINDING_ENERGY 30.9 // MeV
#define BEAM_IS_NUMI false

#include <iostream>
#include <string>

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "configuration.h"
#include "framework.h"
#include "cuts.h"
#include "variables.h"
#include "include/analysis.h"

REGISTER_CUT_SCOPE(RegistrationScope::Both, flash_cut,      cuts::flash_cut);
REGISTER_CUT_SCOPE(RegistrationScope::Both, fiducial_cut,   cuts::fiducial_cut);

REGISTER_VAR_SCOPE(RegistrationScope::True, neutrino_id,    vars::neutrino_id);
REGISTER_VAR_SCOPE(RegistrationScope::Both, flash_time,     vars::flash_time);
REGISTER_VAR_SCOPE(RegistrationScope::Both, visible_energy, vars::visible_energy);

int main(int argc, char * argv[])
{
    // Check if the configuration file is provided as a command line argument
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <configuration_file>" << std::endl;
        return 1;
    }

    // SpectrumLoader
    ana::Analysis analysis("tomlexample");
    ana::SpectrumLoader mc("/pnfs/sbnd/persistent/users/mueller/v10_04_07/nue_v10_04_07.flat.root");
    analysis.AddLoader("mc", &mc, true);

    // Load the configuration file
    sys::cfg::ConfigurationTable config;
    try
    {
        // Load the configuration file
        config.set_config(argv[1]);
        // Check if the required fields are present in the configuration file
        std::vector<sys::cfg::ConfigurationTable> trees(config.get_subtables("tree"));
        for(const auto & tree : trees)
        {
            std::vector<sys::cfg::ConfigurationTable> cuts = tree.get_subtables("cut");
            std::vector<sys::cfg::ConfigurationTable> vars = tree.get_subtables("branch");
            
            std::map<std::string, ana::SpillMultiVar> vars_map;
            for(const auto & var : vars)
            {
                ana::SpillMultiVar thisvar = construct(cuts, var);
                vars_map.try_emplace(var.get_string_field("name"), thisvar);
            }
            analysis.AddTree(tree.get_string_field("name"), vars_map, tree.get_bool_field("sim_only"));
        }

        analysis.Go();
    }
    catch(const sys::cfg::ConfigurationError &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}