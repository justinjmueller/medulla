/**
 * @file main.cc
 * @brief Main file for the SPINE analysis framework.
 * @details This file contains the main function for the SPINE analysis
 * framework. The main function is responsible for loading the configuration
 * file, initializing the analysis framework, and running the analysis.
 * @author mueller@fnal.gov
 * @author rvizarr@fnal.gov
 */
#define PLACEHOLDERVALUE std::numeric_limits<double>::quiet_NaN()
#define PROTON_BINDING_ENERGY 30.9 // MeV
#define BEAM_IS_NUMI true

#include <iostream>
#include <string>
#include <memory>

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "TError.h"

#include "configuration.h"
#include "framework.h"
#include "scorers.h"
#include "cuts.h"
#include "variables.h"
#include "mctruth_variables.h"
#include "mctruth_cuts.h"
#include "bivariables.h"
#include "event_cuts.h"
#include "event_variables.h"
#include "spill_cuts.h"
#include "selectors.h"
#include "biselectors.h"
#include "analysis.h"

std::shared_ptr<VarFn<RParticleType>> pvars::primfn = std::make_shared<VarFn<RParticleType>>(pvars::default_primary_classification<RParticleType>);
std::shared_ptr<VarFn<RParticleType>> pvars::pidfn = std::make_shared<VarFn<RParticleType>>(pvars::default_pid<RParticleType>);

/**
 * @brief Set a function pointer for a variable function.
 * @details This function sets a function pointer for a variable function of
 * type T. It is intended to be used to set scoring function to the user-
 * defined function registered in the framework.
 * @tparam T The type of the variable function.
 * @param fcn A shared pointer to the variable function to be set.
 * @param name The name of the variable function to be set.
 */
template<typename T>
void set_fcn(std::shared_ptr<VarFn<T>> & fcn, const std::string & name)
{
    std::string var_name;
    if constexpr(std::is_same_v<T, RParticleType>)
        var_name = "reco_particle_" + name;
    auto factory = VarFactoryRegistry<T>::instance().get(var_name);
    auto var_fn = factory({});
    fcn = std::make_shared<VarFn<T>>(var_fn);
}

/**
 * @brief An error handler for ROOT errors related to XRootD authentication.
 * @details This function is a custom error handler for ROOT errors. It checks
 * for specific error messages related to XRootD authentication and throws
 * a runtime error with a more user-friendly message if such errors are
 * detected. If the error level is greater than kWarning, it will also
 * call the default error handler to handle other errors.
 * @param level The error level (e.g., kError, kWarning).
 * @param abort Whether to abort the program on error.
 * @param location The location of the error (file and line number).
 * @param message The error message.
 */
void error_handler(int level, bool abort, const char * location, const char * message)
{
    if(level > kWarning)
    {
        // Check for XRootD authentication errors
        if(std::string(message).find("Auth failed: No protocols left to try") != std::string::npos ||
           std::string(message).find("Server responded with an error") != std::string::npos)
        {
            std::string error_message = "Authentication error: No valid token found for XRootD access.";
            error_message += "\n\tPlease ensure you have a valid token with:";
            error_message += "\n\thtgettoken -a htvaultprod.fnal.gov -i <experiment>";
            throw std::runtime_error(error_message);
        }
    }
    ::DefaultErrorHandler(level, abort, location, message);
}

int main(int argc, char * argv[])
{
    // Set the ROOT error handler to our custom error handler. This allows us
    // to catch errors related to XRootD authentication.
    SetErrorHandler(error_handler);

    // Check if the configuration file is provided as a command line argument
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <configuration_file>" << std::endl;
        return 1;
    }

    // Load the configuration file
    cfg::ConfigurationTable config;
    try
    {
        // Load the configuration file
        config.set_config(argv[1]);

        // Construct the "final_state_signal" particle-level cut function.
        if(config.has_field("general.fsthresh"))
        {
            // Retrieve the threshold for final state signal particles.
            std::vector<double> fsthresh = config.get_double_vector("general.fsthresh");

            // Set the global vector for final state signal thresholds.
            pcuts::final_state_signal_thresholds = fsthresh;
        }

        // Each category is a (TType cut, MCTruth cut) pair; the vector index
        // is the value emitted for interactions that match that category.
        CategoryFns category_cut_functions;

        // Construct the category function.
        if(config.has_field("category"))
        {
            // Iterate over the categories and construct the cut functions.
            std::vector<cfg::ConfigurationTable> categories(config.get_subtables("category"));
            for(const auto & category : categories)
            {
                std::vector<CutFn<TType>> true_cut_functions;
                std::vector<CutFn<MCTruth>> mctruth_cut_functions;
                std::vector<cfg::ConfigurationTable> cuts = category.get_subtables("cuts");
                for(const auto & cut : cuts)
                {
                    // Cuts with type = "mctruth" are prefixed with "mctruth_"
                    // and looked up in CutFactoryRegistry<MCTruth>
                    std::string name = cut.get_string_field("name");
                    bool invert = false;
                    if(name.at(0) == '!')
                    {
                        invert = true;
                        name = name.substr(1);
                    }

                    std::vector<double> params;
                    if(cut.has_field("parameters"))
                        params = cut.get_double_vector("parameters");

                    std::string type = cut.get_string_field("type", "true");

                    if(type == "true")
                    {
                        std::string cut_name = "true_" + name;
                        auto factory = CutFactoryRegistry<TType>::instance().get(cut_name);
                        if(invert)
                        {
                            auto fn = factory(params);
                            true_cut_functions.push_back([fn](const TType & e) { return !fn(e); });
                        }
                        else
                            true_cut_functions.push_back(factory(params));
                    }
                    else
                    {
                        std::string cut_name = "mctruth_" + name;
                        auto factory = CutFactoryRegistry<MCTruth>::instance().get(cut_name);
                        if(invert)
                        {
                            auto fn = factory(params);
                            mctruth_cut_functions.push_back([fn](const MCTruth & m) { return !fn(m); });
                        }
                        else
                            mctruth_cut_functions.push_back(factory(params));
                    }
                }

                // Compose the category cut functions. The true_cut function is
                // applied to the true interaction to determine if it passes
                // the category definition, while the mctruth_cut function is
                // applied to the corresponding MCTruth instance.
                auto true_cut = [true_cut_functions](const TType & e) -> bool {
                    return std::all_of(true_cut_functions.begin(), true_cut_functions.end(), [&e](auto & f) { return f(e); });
                };
                auto mctruth_cut = [mctruth_cut_functions](const MCTruth & m) -> bool {
                    return std::all_of(mctruth_cut_functions.begin(), mctruth_cut_functions.end(), [&m](auto & f) { return f(m); });
                };
                category_cut_functions.push_back({true_cut, mctruth_cut});

            }
        }

        // Initialize the analysis framework with the output file name 
        // specified in the configuration.
        ana::Analysis analysis(config.get_string_field("general.output"));

        // Set the PID functions.
        set_fcn(pvars::primfn, config.get_string_field("general.primfn", "default_primary_classification"));
        set_fcn(pvars::pidfn, config.get_string_field("general.pidfn", "default_pid"));

        // Configure the samples in the analysis
        std::vector<cfg::ConfigurationTable> samples = config.get_subtables("sample");
        std::vector<std::unique_ptr<ana::SpectrumLoader>> loaders;
        loaders.reserve(samples.size());
        for(const auto & sample : samples)
        {
            if(sample.get_bool_field("disable", false))
            {
                std::cout << "Sample '" << sample.get_string_field("name") << "' is disabled, skipping." << std::endl;
                continue;
            }

            bool ismc = sample.get_bool_field("ismc", false);
            std::unique_ptr<ana::SpectrumLoader> loader;
            try
            {
                sample.get_string_field("path");
                loader = std::make_unique<ana::SpectrumLoader>(sample.get_string_field("path"));
            }
            catch(const cfg::ConfigurationError &)
            {
                loader = std::make_unique<ana::SpectrumLoader>(sample.get_string_vector("path"));
            }
            analysis.AddLoader(sample.get_string_field("name"), loader.get(), ismc);
            loaders.push_back(std::move(loader));

            // Main loop over the trees defined in the configuration
            std::vector<cfg::ConfigurationTable> trees(config.get_subtables("tree"));
            for(const auto & tree : trees)
            {
                std::vector<cfg::ConfigurationTable> cuts = tree.get_subtables("cut");
                std::vector<cfg::ConfigurationTable> vars = tree.get_subtables("branch");
                std::string mode = tree.get_string_field("mode");

                std::map<std::string, ana::SpillMultiVar> vars_map;
                for(const auto & var : vars)
                {
                    const std::string var_type = var.get_string_field("type");
                    
                    // If the variable type is "both", we need to construct two
                    // variables: one for "true" and one for "reco".
                    if(var_type == "both")
                    {
                        NamedSpillMultiVar thisvar_true = construct(cuts, var, mode, "true", ismc);
                        NamedSpillMultiVar thisvar_reco = construct(cuts, var, mode, "reco", ismc);
                        vars_map.try_emplace(thisvar_true.first, thisvar_true.second);
                        vars_map.try_emplace(thisvar_reco.first, thisvar_reco.second);
                    }
                    else if(var_type == "both_particle")
                    {
                        NamedSpillMultiVar thisvar_true = construct(cuts, var, mode, "true_particle", ismc);
                        NamedSpillMultiVar thisvar_reco = construct(cuts, var, mode, "reco_particle", ismc);
                        vars_map.try_emplace(thisvar_true.first, thisvar_true.second);
                        vars_map.try_emplace(thisvar_reco.first, thisvar_reco.second);
                    }
                    else if(var_type == "both_bivar")
                    {
                        NamedSpillMultiVar thisvar_true = construct(cuts, var, mode, "true_bivar", ismc);
                        NamedSpillMultiVar thisvar_reco = construct(cuts, var, mode, "reco_bivar", ismc);
                        vars_map.try_emplace(thisvar_true.first, thisvar_true.second);
                        vars_map.try_emplace(thisvar_reco.first, thisvar_reco.second);
                    }
                    else if(var_type == "true"
                            || var_type == "reco"
                            || var_type == "mctruth"
                            || var_type == "true_particle"
                            || var_type == "reco_particle"
                            || var_type == "true_bivar"
                            || var_type == "reco_bivar"
                            || var_type == "event")
                    {
                        if(var.get_string_field("name") == "category" && var_type == "true")
                        {
                            NamedSpillMultiVar thisvar = construct_category(cuts, category_cut_functions, mode, ismc);
                            vars_map.try_emplace(thisvar.first, thisvar.second);
                        }
                        else
                        {
                            NamedSpillMultiVar thisvar = construct(cuts, var, mode, var_type, ismc);
                            vars_map.try_emplace(thisvar.first, thisvar.second);
                        }
                    }

                    else
                    {
                        throw std::runtime_error("Illegal variable type '" + var_type + "' for branch " + tree.get_string_field("name") + ":" + var.get_string_field("name"));
                    }
                }
                analysis.AddTreeForSample(sample.get_string_field("name"), tree.get_string_field("name"), vars_map, tree.get_bool_field("sim_only"));

                // Add the exposure tree.
                if(tree.get_bool_field("add_exposure", false))
                {
                    // Construct the exposure variables.
                    std::map<std::string, ana::SpillMultiVar> exposure_vars_map;
                    std::vector<NamedSpillMultiVar> exposure_vars = construct_exposure_vars(cuts);

                    // Add the exposure variables to the map.
                    for(const auto & exposure_var : exposure_vars)
                        exposure_vars_map.try_emplace(exposure_var.first, exposure_var.second);

                    // Add the exposure tree for the sample.
                    analysis.AddTreeForSample(sample.get_string_field("name"), tree.get_string_field("name")+"_exposure", exposure_vars_map, tree.get_bool_field("sim_only"));
                }
            }
        }

        analysis.Go();
    }
    catch(const cfg::ConfigurationError &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
