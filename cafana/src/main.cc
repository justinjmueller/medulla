/**
 * @file main.cc
 * @brief Main file for the SPINE analysis framework.
 * @details This file contains the main function for the SPINE analysis
 * framework. The main function is responsible for loading the configuration
 * file, initializing the analysis framework, and running the analysis.
 * @author mueller@fnal.gov
 */
#define PLACEHOLDERVALUE std::numeric_limits<double>::quiet_NaN()
#define PRIMARYFUNC pvars::lax_primary_classification
#define PIDFUNC pvars::pid
#define PROTON_BINDING_ENERGY 30.9 // MeV
#define BEAM_IS_NUMI false

#include <iostream>
#include <string>
#include <memory>

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "configuration.h"
#include "framework.h"
#include "cuts.h"
#include "muon2024/cuts_muon2024.h"
#include "variables.h"
#include "muon2024/variables_muon2024.h"
#include "mctruth.h"
#include "include/analysis.h"

int main(int argc, char * argv[])
{
    // Check if the configuration file is provided as a command line argument
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <configuration_file>" << std::endl;
        return 1;
    }

    // Load the configuration file
    sys::cfg::ConfigurationTable config;
    try
    {
        // Load the configuration file
        config.set_config(argv[1]);

        // SpectrumLoader
        ana::Analysis analysis(config.get_string_field("general.output"));

        // Configure the samples in the analysis
        std::vector<sys::cfg::ConfigurationTable> samples = config.get_subtables("sample");
        std::vector<std::unique_ptr<ana::SpectrumLoader>> loaders;
        loaders.reserve(samples.size());
        for(const auto & sample : samples)
        {
            // Create a SpectrumLoader for each sample
            std::unique_ptr<ana::SpectrumLoader> loader = std::make_unique<ana::SpectrumLoader>(sample.get_string_field("path"));
            analysis.AddLoader(sample.get_string_field("name"), loader.get(), sample.get_bool_field("ismc"));
            loaders.push_back(std::move(loader));
        }

        // Main loop over the trees defined in the configuration
        std::vector<sys::cfg::ConfigurationTable> trees(config.get_subtables("tree"));
        for(const auto & tree : trees)
        {
            std::vector<sys::cfg::ConfigurationTable> cuts = tree.get_subtables("cut");
            std::vector<sys::cfg::ConfigurationTable> vars = tree.get_subtables("branch");
            std::string mode = tree.get_string_field("mode");
            
            std::map<std::string, ana::SpillMultiVar> vars_map;
            for(const auto & var : vars)
            {
                // If the variable type is "both", we need to construct two
                // variables: one for "true" and one for "reco".
                if(var.get_string_field("type") == "both")
                {
                    NamedSpillMultiVar thisvar_true = construct(cuts, var, mode, "true");
                    NamedSpillMultiVar thisvar_reco = construct(cuts, var, mode, "reco");
                    vars_map.try_emplace(thisvar_true.first, thisvar_true.second);
                    vars_map.try_emplace(thisvar_reco.first, thisvar_reco.second);
                }
                else if(var.get_string_field("type") == "both_particle")
                {
                    NamedSpillMultiVar thisvar_true = construct(cuts, var, mode, "true_particle");
                    NamedSpillMultiVar thisvar_reco = construct(cuts, var, mode, "reco_particle");
                    vars_map.try_emplace(thisvar_true.first, thisvar_true.second);
                    vars_map.try_emplace(thisvar_reco.first, thisvar_reco.second);
                }
                else if(var.get_string_field("type") == "true"
                        || var.get_string_field("type") == "reco"
                        || var.get_string_field("type") == "mctruth"
                        || var.get_string_field("type") == "true_particle"
                        || var.get_string_field("type") == "reco_particle")
                {
                    NamedSpillMultiVar thisvar = construct(cuts, var, mode);
                    vars_map.try_emplace(thisvar.first, thisvar.second);
                }
                else
                {
                    throw std::runtime_error("Illegal variable type '" + var.get_string_field("type") + "' for branch " + tree.get_string_field("name") +  ":" + var.get_string_field("name"));
                }
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