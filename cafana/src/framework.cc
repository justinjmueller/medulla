/**
 * @file framework.cc
 * @brief Implementation of the SPINE analysis framework.
 * @details This file contains the implementation of the SPINE analysis
 * framework. The framework is designed to be modular and extensible, allowing for
 * easy integration and application of cuts and variables.
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

constexpr size_t kNoMatch = std::numeric_limits<size_t>::max();

// Get the singleton instance of the Registry.
template<typename EventT, typename RegistryT>
Registry<EventT, RegistryT> & Registry<EventT, RegistryT>::instance()
{
    static Registry<EventT, RegistryT> instance;
    return instance;
}

// Check if a function is registered under the specified name.
template<typename EventT, typename RegistryT>
bool Registry<EventT, RegistryT>::is_registered(const std::string & name) const
{
    // Check if the function is registered
    return registry_.find(name) != registry_.end();
}

// Register under the specified name.
template<typename EventT, typename RegistryT>
void Registry<EventT, RegistryT>::register_fn(const std::string & name, RegistryT (*fn)(const EventT &))
{
    // Check if the function is already registered
    if(is_registered(name))
    {
        throw std::runtime_error("Function " + name + " is already registered.");
    }
    // Register the function
    registry_[name] = fn;
}

// Retrieve a previously registered function by name.
template<typename EventT, typename RegistryT>
typename Registry<EventT, RegistryT>::Fn Registry<EventT, RegistryT>::get(const std::string & name)
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
NamedSpillMultiVar construct(const std::vector<sys::cfg::ConfigurationTable> & cuts, const sys::cfg::ConfigurationTable & var, const std::string & override_type)
{
    /**
     * @brief Determine the type of the cuts.
     * @details The type of the first cut is used to determine the type of the
     * cuts, and therefore whether the selection is applied in a loop over true
     * or reco events. This type is used to branch the code into the appropriate
     * path.
     */
    bool is_mc(cuts[0].get_string_field("type") == "true");

    if(is_mc)
    {
        /**
         * @brief Compose a common cut function.
         * @details This function composes a common cut function from the
         * subtables of the cuts vector. This is a logical "and" of all
         * configured cuts constructed using std::all_of.
         */
        std::vector<std::function<bool(const TType &)>> cut_functions;
        cut_functions.reserve(cuts.size());
        for(const auto & cut : cuts)
        {
            std::string cut_name = "true_" + cut.get_string_field("name");
            std::vector<double> params;
            if(cut.has_field("parameters"))
              params = cut.get_double_vector("parameters");
            if(CutFactoryRegistry<TType>::Instance().Create(cut_name, params).target_type() != typeid(void))
            {
                cut_functions.push_back(CutFactoryRegistry<TType>::Instance().Create(cut_name, params));
            }
            else if(CutRegistry<TType>::instance().is_registered(cut_name))
            {
                cut_functions.push_back(CutRegistry<TType>::instance().get(cut_name));
            }
            else
            {
                throw std::runtime_error("Cut " + cut_name + " is not registered.");
            }
        }
        auto cut = [cut_functions](const TType & e) -> bool {
            return std::all_of(cut_functions.begin(), cut_functions.end(), [&e](auto & f) { return f(e); });
        };

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

        if(var_type == "true") {
            var_name = "true_" + var_name;
            auto varFn = VarFactoryRegistry<TType>::Instance().Create(var_name, varPars);
            return std::make_pair(var_name, spill_multivar_helper<TType, TType>(cut, varFn));
        } else if(var_type == "reco") {
            var_name = "reco_" + var_name;
            auto varFn = VarFactoryRegistry<RType>::Instance().Create(var_name, varPars);
            return std::make_pair(var_name, spill_multivar_helper<TType, RType>(cut, varFn));
        } else {
            throw std::runtime_error("Illegal variable type '" + var_type + "' for variable " + var_name);
        }
    }
    else
    {
        /**
         * @brief Compose a common cut function.
         * @details This function composes a common cut function from the
         * subtables of the cuts vector. This is a logical "and" of all
         * configured cuts constructed using std::all_of.
         */
        std::vector<std::function<bool(const RType &)>> cut_functions;
        cut_functions.reserve(cuts.size());
        for(const auto & cut : cuts)
        {
            std::string cut_name = "reco_" + cut.get_string_field("name");
            std::vector<double> params;
            if(cut.has_field("parameters"))
              params = cut.get_double_vector("parameters");
            if(CutFactoryRegistry<RType>::Instance().Create(cut_name, params).target_type() != typeid(void))
            {
                cut_functions.push_back(CutFactoryRegistry<RType>::Instance().Create(cut_name, params));
            }
            else if(CutRegistry<RType>::instance().is_registered(cut_name))
            {
                cut_functions.push_back(CutRegistry<RType>::instance().get(cut_name));
            }
            else
            {
                throw std::runtime_error("Cut " + cut_name + " is not registered.");
            }
        }
        auto cut = [cut_functions](const RType & e) -> bool {
            return std::all_of(cut_functions.begin(), cut_functions.end(), [&e](auto & f) { return f(e); });
        };

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

        if(var_type == "true") {
            var_name = "true_" + var_name;
            auto varFn = VarFactoryRegistry<TType>::Instance().Create(var_name, varPars);
            return std::make_pair(var_name, spill_multivar_helper<RType, TType>(cut, varFn));
        } else if(var_type == "reco") {
            var_name = "reco_" + var_name;
            auto varFn = VarFactoryRegistry<RType>::Instance().Create(var_name, varPars);
            return std::make_pair(var_name, spill_multivar_helper<RType, RType>(cut, varFn));
        } else {
            throw std::runtime_error("Illegal variable type '" + var_type + "' for variable " + var_name);
        }
    }
}

// Helper method for constructing a SpillMultiVar object.
template<typename CutsOn, typename VarOn>
ana::SpillMultiVar spill_multivar_helper(std::function<bool(const CutsOn &)> cuts, std::function<double(const VarOn &)> var)
{
    return ana::SpillMultiVar([cuts, var](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>
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
                    if(cuts(i) && match_id != kNoMatch)
                    {
                        values.push_back(var(sr->dlp[match_id]));
                    }
                }
                else if constexpr(std::is_same_v<VarOn, TType>)
                {
                    if(cuts(i) && match_id != kNoMatch)
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
                    if(cuts(i) && match_id != kNoMatch)
                    {
                        values.push_back(var(sr->dlp_true[match_id]));
                    }
                }
                else if constexpr(std::is_same_v<VarOn, RType>)
                {
                    if(cuts(i) && match_id != kNoMatch)
                    {
                        values.push_back(var(i));
                    }
                }
            }
        }
        return values;
    });
}

using Truth = caf::Proxy<caf::SRInteractionTruthDLP>;
using Reco  = caf::Proxy<caf::SRInteractionDLP>;
template class Registry<Truth,bool>;
template class Registry<Reco, bool>;
template class Registry<Truth,double>;
template class Registry<Reco, double>;