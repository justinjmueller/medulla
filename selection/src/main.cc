/**
 * @file main.cc
 * @brief Main file for the SPINE analysis framework.
 * @details This file contains the main function for the SPINE analysis
 * framework. The main function is responsible for loading the configuration
 * file, initializing the analysis framework, and running the analysis.
 * @author mueller@fnal.gov
 */
#define PLACEHOLDERVALUE std::numeric_limits<double>::quiet_NaN()
#define PROTON_BINDING_ENERGY 30.9 // MeV
#define BEAM_IS_NUMI false

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
#include "mctruth.h"
#include "event_cuts.h"
#include "event_variables.h"
#include "spill_cuts.h"
#include "selectors.h"
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

        // Construct the category function.
        if(config.has_field("category"))
        {
            // Map of category enumeration to cut functions.
            std::map<double, CutFn<TType>> category_cut_functions;

            // Iterate over the categories and construct the cut functions.
            std::vector<cfg::ConfigurationTable> categories(config.get_subtables("category"));
            for(const auto & category : categories)
            {
                std::vector<CutFn<TType>> true_cut_functions;
                std::vector<cfg::ConfigurationTable> cuts = category.get_subtables("cuts");
                for(const auto & cut : cuts)
                {
                    // Retrieve the cut name and check for negation.
                    std::string name = cut.get_string_field("name");
                    bool invert = false;
                    if(name.at(0) == '!')
                    {
                        invert = true;
                        name = name.substr(1); // Remove the negation character.
                    }
                    name = "true_" + name;

                    // Load parameters (if any) for the cut.
                    std::vector<double> params;
                    if(cut.has_field("parameters"))
                        params = cut.get_double_vector("parameters");

                    auto factory = CutFactoryRegistry<TType>::instance().get(name);
                    if(invert)
                    {
                        // If the cut is inverted, we need to negate the function.
                        auto fn = factory(params);
                        true_cut_functions.push_back([fn](const TType & e) { return !fn(e); });
                    }
                    else
                        // Otherwise, we just add the function as is.
                        true_cut_functions.push_back(factory(params));
                }
                // Compose a common cut function for the category.
                auto category_cut = [true_cut_functions](const TType & e) -> bool {
                    return std::all_of(true_cut_functions.begin(), true_cut_functions.end(), [&e](auto & f) { return f(e); });
                };
                category_cut_functions.try_emplace(
                    category_cut_functions.size(),
                    category_cut
                );
            }

            // Create the category function.
            auto category_fn = [category_cut_functions](const TType & e) -> double
            {
                // Iterate over the category cut functions and return the first
                // one that returns true.
                for(const auto & [category, cut_fn] : category_cut_functions)
                {
                    if(cut_fn(e))
                        return category; // Return the category number.
                }
                return PLACEHOLDERVALUE; // No category matched.
            };
            // Register the category function.
            VarFactoryRegistry<TType>::instance().register_fn(
                "true_category",
                [category_fn](const std::vector<double>&) -> VarFn<TType> { return category_fn; }
            );
        }

        // SpectrumLoader
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
            // Check if the sample has the "disable" flag set to true
            if(sample.get_bool_field("disable", false))
            {
                std::cout << "Sample '" << sample.get_string_field("name") << "' is disabled, skipping." << std::endl;
                continue;
            }

            // Create a SpectrumLoader for each sample
            std::unique_ptr<ana::SpectrumLoader> loader = std::make_unique<ana::SpectrumLoader>(sample.get_string_field("path"));
            analysis.AddLoader(sample.get_string_field("name"), loader.get(), sample.get_bool_field("ismc"));
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
                    // If the variable type is "both", we need to construct two
                    // variables: one for "true" and one for "reco".
                    if(var.get_string_field("type") == "both")
                    {
                        NamedSpillMultiVar thisvar_true = construct(cuts, var, mode, "true", sample.get_bool_field("ismc"));
                        NamedSpillMultiVar thisvar_reco = construct(cuts, var, mode, "reco", sample.get_bool_field("ismc"));
                        vars_map.try_emplace(thisvar_true.first, thisvar_true.second);
                        vars_map.try_emplace(thisvar_reco.first, thisvar_reco.second);
                    }
                    else if(var.get_string_field("type") == "both_particle")
                    {
                        NamedSpillMultiVar thisvar_true = construct(cuts, var, mode, "true_particle", sample.get_bool_field("ismc"));
                        NamedSpillMultiVar thisvar_reco = construct(cuts, var, mode, "reco_particle", sample.get_bool_field("ismc"));
                        vars_map.try_emplace(thisvar_true.first, thisvar_true.second);
                        vars_map.try_emplace(thisvar_reco.first, thisvar_reco.second);
                    }
                    else if(var.get_string_field("type") == "true"
                            || var.get_string_field("type") == "reco"
                            || var.get_string_field("type") == "mctruth"
                            || var.get_string_field("type") == "true_particle"
                            || var.get_string_field("type") == "reco_particle"
                            || var.get_string_field("type") == "event")
                    {
                        NamedSpillMultiVar thisvar = construct(cuts, var, mode, var.get_string_field("type"), sample.get_bool_field("ismc"));
                        vars_map.try_emplace(thisvar.first, thisvar.second);
                    }
                    else
                    {
                        throw std::runtime_error("Illegal variable type '" + var.get_string_field("type") + "' for branch " + tree.get_string_field("name") +  ":" + var.get_string_field("name"));
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