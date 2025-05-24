/**
 * @file framework.cc
 * @brief Implementation of the components of the SPINE analysis framework.
 * @details This file contains the implementation of the SPINE analysis
 * framework. The framework is designed to be modular and extensible, allowing
 * for easy integration and application of cuts and variables.
 * @author mueller@fnal.gov
 */
#include <map>
#include <string>
#include <functional>
#include <stdexcept>

#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "framework.h"
#include "configuration.h"

// Set a sensible default for a no-match scenario.
constexpr size_t kNoMatch = std::numeric_limits<size_t>::max();

// Get the singleton instance of the Registry.
template<typename ValueT>
Registry<ValueT> & Registry<ValueT>::instance()
{
    static Registry instance;
    return instance;
}

// Check if a function is registered under the specified name.
template<typename ValueT>
bool Registry<ValueT>::is_registered(const std::string & name) const
{
    // Check if the function is registered
    return registry_.find(name) != registry_.end();
}

// Register under the specified name.
template<typename ValueT>
void Registry<ValueT>::register_fn(const std::string & name, ValueT fn)
{
    // Check if the function is already registered
    if(is_registered(name))
    {
        throw std::runtime_error("Function " + name + " is already registered.");
    }
    // Register the function
    registry_[name] = std::move(fn);
}

// Retrieve a previously registered function by name.
template<typename ValueT>
ValueT Registry<ValueT>::get(const std::string & name)
{
    // Check if the function is registered
    if(!is_registered(name))
    {
        throw std::runtime_error("Function " + name + " is not registered.");
    }
    // Retrieve the function
    return registry_[name];
}

// Build a single SpillMultiVar for a single branch variable.
NamedSpillMultiVar construct(const std::vector<sys::cfg::ConfigurationTable> & cuts,
                             const sys::cfg::ConfigurationTable & var,
                             const std::string & mode,
                             const std::string & override_type)
{
    /**
     * @brief Determine the type of the cuts.
     * @details The type of the first cut is used to determine the type of the
     * cuts, and therefore whether the selection is applied in a loop over true
     * or reco events. This type is used to branch the code into the appropriate
     * path. First, we check if the mode is a valid option. If not, we throw an
     * exception.
     */
    Mode exec_mode;
    if(mode == "true") exec_mode = Mode::True;
    else if(mode == "reco") exec_mode = Mode::Reco;
    else throw std::runtime_error("Illegal mode '" + mode + "' for variable " + var.get_string_field("name"));

    std::vector<CutFn<TType>> true_cut_functions;
    std::vector<CutFn<RType>> reco_cut_functions;
    for(const auto & cut : cuts)
    {
        if(!cut.has_field("type"))
            throw std::runtime_error("Cut " + cut.get_string_field("name") + " does not have a type field.");
        if(cut.get_string_field("type") == "true")
        {
            std::string cut_name = "true_" + cut.get_string_field("name");
            std::vector<double> params;
            if(cut.has_field("parameters"))
                params = cut.get_double_vector("parameters");
            auto factory = CutFactoryRegistry<TType>::instance().get(cut_name);
            true_cut_functions.push_back(factory(params));
        }
        else if(cut.get_string_field("type") == "reco")
        {
            std::string cut_name = "reco_" + cut.get_string_field("name");
            std::vector<double> params;
            if(cut.has_field("parameters"))
                params = cut.get_double_vector("parameters");
            auto factory = CutFactoryRegistry<RType>::instance().get(cut_name);
            reco_cut_functions.push_back(factory(params));
        }
        else
        {
            throw std::runtime_error("Illegal cut type '" + cut.get_string_field("type") + "' for cut " + cut.get_string_field("name"));
        }
    }

    /**
     * @brief Compose a common cut function.
     * @details This function composes a common cut function from the
     * subtables of the cuts vector. This is a logical "and" of all
     * configured cuts constructed using std::all_of.
     */
    auto true_cut = [true_cut_functions](const TType & e) -> bool {
        return std::all_of(true_cut_functions.begin(), true_cut_functions.end(), [&e](auto & f) { return f(e); });
    };
    auto reco_cut = [reco_cut_functions](const RType & e) -> bool {
        return std::all_of(reco_cut_functions.begin(), reco_cut_functions.end(), [&e](auto & f) { return f(e); });
    };

    if(exec_mode == Mode::True)
    {
        /**
         * @brief Read the branch variable configuration.
         * @details This function constructs the branch variable from the TOML
         * configuration of the variable. The variable name is used to retrieve the
         * function from the registry.
         */
        std::string var_name = var.get_string_field("name");
        std::string var_type = (override_type.empty() ? var.get_string_field("type") : override_type);
        std::vector<double> varPars;
        if(var.has_field("parameters"))
            varPars = var.get_double_vector("parameters");

        if(var_type == "true")
        {
            var_name = "true_" + var_name;
            auto factory = VarFactoryRegistry<TType>::instance().get(var_name);
            auto varFn = factory(varPars);
            return std::make_pair(var_name, spill_multivar_helper<TType, RType, TType>(true_cut, reco_cut, varFn));
        }
        else if(var_type == "reco")
        {
            var_name = "reco_" + var_name;
            auto factory = VarFactoryRegistry<RType>::instance().get(var_name);
            auto varFn = factory(varPars);
            return std::make_pair(var_name, spill_multivar_helper<TType, RType, RType>(true_cut, reco_cut, varFn));
        }
        else
        {
            throw std::runtime_error("Illegal variable type '" + var_type + "' for variable " + var_name);
        }
    }
    else
    {
        /**
         * @brief Read the branch variable configuration.
         * @details This function constructs the branch variable from the TOML
         * configuration of the variable. The variable name is used to retrieve the
         * function from the registry.
         */
        std::string var_name = var.get_string_field("name");
        std::string var_type = (override_type.empty() ? var.get_string_field("type") : override_type);
        std::vector<double> varPars;
        if(var.has_field("parameters"))
            varPars = var.get_double_vector("parameters");

        if(var_type == "true")
        {
            var_name = "true_" + var_name;
            auto factory = VarFactoryRegistry<TType>::instance().get(var_name);
            auto varFn = factory(varPars);
            return std::make_pair(var_name, spill_multivar_helper<RType, TType, TType>(reco_cut, true_cut, varFn));
        }
        else if(var_type == "reco")
        {
            var_name = "reco_" + var_name;
            auto factory = VarFactoryRegistry<RType>::instance().get(var_name);
            auto varFn = factory(varPars);
            return std::make_pair(var_name, spill_multivar_helper<RType, TType, RType>(reco_cut, true_cut, varFn));
        }
        else
        {
            throw std::runtime_error("Illegal variable type '" + var_type + "' for variable " + var_name);
        }
    }
}

// Helper method for constructing a SpillMultiVar object.
template<typename CutsOn, typename CompsOn, typename VarOn>
ana::SpillMultiVar spill_multivar_helper(
  CutFn<CutsOn> cuts,
  CutFn<CompsOn> comps,
  VarFn<VarOn> var)
{
    return ana::SpillMultiVar([comps, cuts, var](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>
    {
        std::vector<double> values;
        if constexpr (std::is_same_v<CutsOn, TType>)
        {
            for(auto const& i : sr->dlp_true)
            {
                // Check for match
                size_t match_id = (i.match_ids.size() > 0) ? (size_t)i.match_ids[0] : kNoMatch;

                if constexpr(std::is_same_v<VarOn, RType>)
                {
                    if(cuts(i) && match_id != kNoMatch && comps(sr->dlp[match_id]))                    
                    {
                        values.push_back(var(sr->dlp[match_id]));
                    }
                }
                else if constexpr(std::is_same_v<VarOn, TType>)
                {
                    if(cuts(i) && match_id != kNoMatch && comps(sr->dlp[match_id]))
                    {
                        values.push_back(var(i));
                    }
                }
            }
        }
        else if constexpr(std::is_same_v<CutsOn, RType>)
        {
            for(auto const& i : sr->dlp)
            {
                // Check for match
                size_t match_id = (i.match_ids.size() > 0) ? (size_t)i.match_ids[0] : kNoMatch;

                if constexpr(std::is_same_v<VarOn, TType>)
                {
                    if(cuts(i) && match_id != kNoMatch && comps(sr->dlp_true[match_id]))
                    {
                        values.push_back(var(sr->dlp_true[match_id]));
                    }
                }
                else if constexpr(std::is_same_v<VarOn, RType>)
                {
                    if(cuts(i) && match_id != kNoMatch && comps(sr->dlp_true[match_id]))
                    {
                        values.push_back(var(i));
                    }
                }
            }
        }
        return values;
    });
}

// Explicitly instantiate Registry for the factory types we use:
template class Registry<CutFactory<TType>>;
template class Registry<CutFactory<RType>>;
template class Registry<VarFactory<TType>>;
template class Registry<VarFactory<RType>>;