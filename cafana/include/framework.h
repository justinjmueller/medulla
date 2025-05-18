/**
 * @file framework.h
 * @brief Header file for the SPINE analysis framework.
 * @details This file contains the header for the SPINE analysis framework. The
 * framework is designed to be modular and extensible, allowing for easy
 * integration and application of cuts and variables.
 * @author mueller@fnal.gov
 */
#ifndef FRAMEWORK_H
#define FRAMEWORK_H
#include <map>
#include <string>
#include <functional>
#include <stdexcept>

#include "sbnana/CAFAna/Core/MultiVar.h"
#include "configuration.h"

typedef caf::SRInteractionTruthDLPProxy TType;
typedef caf::SRInteractionDLPProxy RType;

template<typename EventT, typename RegistryT>
class Registry
{
    public:
        /**
         * @brief The function signature on EventT.
         */
        using Fn = std::function<RegistryT(const EventT &)>;

        /**
         * @brief Get the singleton instance of the Registry.
         * @details This function returns a reference to the singleton instance
         * of the Registry. The instance is created on the first call to this
         * function.
         * @return A reference to the singleton instance of the Registry.
         */
        static Registry & instance();

        /**
         * @brief Check if a function is registered under the specified name.
         * @details This function checks if a function is registered under the
         * specified name. This is intended to be used both to make sure that
         * lookups are not done on unregistered functions and to ensure that a
         * name is not reused.
         * @param name The name to check.
         * @return A boolean value indicating whether the function is registered.
         */
        bool is_registered(const std::string & name) const;

        /**
         * @brief Register under the specified name.
         * @details This function registers under the specified name. The name
         * is used to identify the function in the TOML-based configuration
         * file.
         * @param name The name to register the function under.
         * @param fn The function to register.
         * @throw std::runtime_error if the function is already registered.
         */
        void register_fn(const std::string & name, RegistryT (*fn)(const EventT &));

        /**
         * @brief Retrieve a previously registered function by name.
         * @details This function retrieves a previously registered function by
         * name. The name is used to identify the function in the TOML-based
         * configuration file.
         * @param name The name of the function to retrieve.
         * @return A reference to the function 
         * @throw std::runtime_error if the function is not registered.
         */
        Fn get(const std::string & name);

        private:
        /**
         * @brief The registry of functions.
         * @details This is a map of function names to function pointers. The
         * function names are used to identify the functions in the TOML-based
         * configuration file.
         */
        std::map<std::string, Fn> registry_;
};

/**
 * @brief Alias for the registry of Cuts.
 * @details There are two specific registeries employed in this
 * framework: one for Cuts (return bool) and one for Variables (return
 * double). This one is for Cuts.
 * @tparam EventT The type of event.
 */
template<typename EventT>
using CutRegistry = Registry<EventT, bool>;

/**
 * @brief Alias for the registry of Variables.
 * @details There are two specific registeries employed in this
 * framework: one for Cuts (return bool) and one for Variables (return
 * double). This one is for Variables.
 * @tparam EventT The type of event.
 */
template<typename EventT>
using VariableRegistry = Registry<EventT, double>;

/**
 * @brief Macro to register a cut on both true and reco types.
 * @param fn The templated function name (no angle brackets).
 *
 * Expands to:
 *   CutRegistry<TType>::instance().register_fn("true_fn", &fn<TType>);
 *   CutRegistry<RType>::instance().register_fn("reco_fn", &fn<RType>);
 */
#define REGISTER_CUT(fn)                                                           \
namespace                                                                          \
{                                                                                  \
    struct RegistrarFor##fn {                                                      \
    RegistrarFor##fn()                                                             \
    {                                                                              \
        CutRegistry<TType>::instance().register_fn("true_" #fn, fn<TType>);        \
        CutRegistry<RType>::instance().register_fn("reco_" #fn, fn<RType>);        \
    }                                                                              \
    } registrarFor##fn;                                                            \
}

/**
 * @brief Macro to register a variable on both true and reco types.
 * @param fn The templated function name (no angle brackets).
 *
 * Expands to:
 *   VariableRegistry<TType>::instance().register_fn("true_fn", &fn<TType>);
 *   VariableRegistry<RType>::instance().register_fn("reco_fn", &fn<RType>);
 */
#define REGISTER_VARIABLE(fn)                                                      \
namespace                                                                          \
{                                                                                  \
    struct RegistrarFor##fn {                                                      \
    RegistrarFor##fn()                                                             \
    {                                                                              \
        VariableRegistry<TType>::instance().register_fn("true_" #fn, fn<TType>);   \
        VariableRegistry<RType>::instance().register_fn("reco_" #fn, fn<RType>);   \
    }                                                                              \
    } registrarFor##fn;                                                            \
}

/**
 * @brief Build a single SpillMultiVar for a single branch variable.
 * @details Applies the sequence of Cuts from @p cuts to select events, then
 * computes the variable as defined by @p var. Handles "true" and "reco" types
 * by prefixing the branch name and selecting the appropriate event types.
 * @param cuts Vector of [[tree.cut]] subtables with fields:
 *        - name:       string (base cut name)
 *        - type:       string ("true" or "reco")
 *        - parameters: array of floats (parameters for the cut)
 * @param var [[tree.variable]] subtable with fields:
 *        - name:       string (base variable name)
 *        - type:       string ("true" or "reco")
 *        - parameters: array of floats (parameters for the variable)
 * @return A SpillMultiVar object that applies the cuts and computes the variable.
 * @throw std::runtime_error if a function is not registered.
 */
ana::SpillMultiVar construct(const std::vector<sys::cfg::ConfigurationTable> & cuts, const sys::cfg::ConfigurationTable & var);

/**
 * @brief Helper method for constructing a SpillMultiVar object.
 * @details This function is used to construct a SpillMultiVar object from
 * the parameter Cut and Variable objects. It is intended to be called by the
 * @ref construct function.
 * @tparam CutsOn The type (TType or RType) that the cut is applied to.
 * @tparam VarOn The type (TType or RType) that the variable is applied to.
 * @param cuts The callable that implements the cut.
 * @param var The callable that implements the variable.
 * @return A SpillMultiVar object that applies the cuts and computes the variable.
 */
template<typename CutsOn, typename VarOn>
ana::SpillMultiVar spill_multivar_helper(std::function<bool(const CutsOn &)> cuts, std::function<double(const VarOn &)> var);

#endif // FRAMEWORK_H