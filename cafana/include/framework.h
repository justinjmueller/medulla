/**
 * @file framework.h
 * @brief Header file declaring the components of the SPINE analysis framework.
 * @details This file contains the header declarations of the SPINE analysis
 * framework. The framework is designed to be modular and extensible, allowing
 * for easy integration and applications of cuts and variables.
 * @author mueller@fnal.gov
 */
#ifndef FRAMEWORK_H
#define FRAMEWORK_H
#include <map>
#include <vector>
#include <string>
#include <functional>
#include <stdexcept>
#include <type_traits>

#include "sbnana/CAFAna/Core/MultiVar.h"
#include "configuration.h"

/**
 * @brief Type aliases for the event types used in the framework.
 * @details These type aliases are used to simplify the code and make it
 * easier to read. The TType and RType are used to represent the "true" and
 * "reco" event types, respectively. Both refer to a "Proxy" wrapped object.
 */
using TType = caf::SRInteractionTruthDLPProxy;
using RType = caf::SRInteractionDLPProxy;
using MCTruth = caf::Proxy<caf::SRTrueInteraction>;

using NamedSpillMultiVar = std::pair<std::string, ana::SpillMultiVar>;

//-----------------------------------------------------------------------------
// 1) Generic registry template
//-----------------------------------------------------------------------------
/**
 * @brief A singleton registory mapping string names to callables.
 * @details This class is a singleton that maps string names to callable
 * functions. Its primary use is to allow for TOML-based lookup of functions
 * to be called in the analysis framework.
 * @tparam ValueT The type of the value to be registered. This can be a function
 * pointer, a lambda, or any other callable type.
 */
template<typename ValueT>
class Registry
{
    public:
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
        void register_fn(const std::string & name, ValueT fn);

        /**
         * @brief Retrieve a previously registered function by name.
         * @details This function retrieves a previously registered function by
         * name. The name is used to identify the function in the TOML-based
         * configuration file.
         * @param name The name of the function to retrieve.
         * @return A copy of the registered Fn. 
         * @throw std::runtime_error if the function is not registered.
         */
        ValueT get(const std::string & name);

        private:
        /**
         * @brief The registry of functions.
         * @details This is a map of function names to function pointers. The
         * function names are used to identify the functions in the TOML-based
         * configuration file.
         */
        std::map<std::string, ValueT> registry_;
};

//-----------------------------------------------------------------------------
// 2) Raw function registries
//-----------------------------------------------------------------------------
/**
 * @brief Alias for raw Cut functions with signature bool(const EventT&).
 */
template<typename EventT>
using CutFn = std::function<bool(const EventT&)>;

/**
 * @brief Alias for raw Variable functions with signature double(const EventT&).
 */
template<typename EventT>
using VarFn = std::function<double(const EventT&)>;

//-----------------------------------------------------------------------------
// 3) Factory function registries
//-----------------------------------------------------------------------------
/**
 * @brief A factory function: given params, returns a CutFn<EventT>
 */
template<typename EventT>
using CutFactory = std::function<CutFn<EventT>(const std::vector<double>&)>;

template<typename EventT>
using CutFactoryRegistry = Registry<CutFactory<EventT>>;

/**
 * @brief A factory function: given params, returns a VarFn<EventT>
 */
template<typename EventT>
using VarFactory = std::function<VarFn<EventT>(const std::vector<double>&)>;

template<typename EventT>
using VarFactoryRegistry = Registry<VarFactory<EventT>>;

/**
 * @brief Bind a function to a specific parameter set.
 * @details This function binds a function to a specific parameter set. It uses
 * std::is_invocable_v to check if the function can accept a vector of
 * parameters; if so, it binds the function with the parameters. In either
 * case, it returns a std::function<ValueT(const EventT&)> that can be used to
 * apply the function to an event.
 * @tparam F The function to bind.
 * @tparam EventT The type of event: @ref TType or @ref RType.
 * @tparam ValueT The return type of the function.
 * @param pars The parameters to bind to the function.
 * @return A std::function<ValueT(const EventT&)> that applies the cut to an
 * event.
 */
template<auto F, typename EventT, typename ValueT>
inline std::function<ValueT(const EventT&)> bind(const std::vector<double>& pars)
{
    if constexpr(std::is_invocable_v<decltype(F), const EventT&, const std::vector<double>&>)
        return [pars](const EventT& e){ return F(e, pars); };
    else
        return [=](const EventT& e){ return F(e); };
}

/**
 * @brief Scope for registration macros
 * @details This enum class defines the scope of registration for cuts and
 * variables. It can be either "true", "reco," "both," or "mctruth". This is
 * used in the registration macros to determine which type of event the cut or
 * variable can be reasonably applied to.
 */
enum class RegistrationScope { True, Reco, Both, MCTruth };

// Register a cut with scope, auto‐detecting its signature
#define REGISTER_CUT_SCOPE(scope, name, fn)                                                \
namespace                                                                                  \
{                                                                                          \
    const bool _reg_cut_##name = []{                                                       \
        if constexpr((scope)==RegistrationScope::True || (scope)==RegistrationScope::Both) \
            CutFactoryRegistry<TType>::instance().register_fn(                             \
                "true_" #name, bind<+fn<TType>, TType, bool>                               \
            );                                                                             \
        if constexpr((scope)==RegistrationScope::Reco || (scope)==RegistrationScope::Both) \
            CutFactoryRegistry<RType>::instance().register_fn(                             \
                "reco_" #name, bind<+fn<RType>, RType, bool>                               \
            );                                                                             \
        return true;                                                                       \
    }();                                                                                   \
}

// Register a variable with scope, auto‐detecting its signature
#define REGISTER_VAR_SCOPE(scope, name, fn)                                                \
namespace                                                                                  \
{                                                                                          \
    const bool _reg_var_##name = []{                                                       \
        if constexpr((scope)==RegistrationScope::True || (scope)==RegistrationScope::Both) \
            VarFactoryRegistry<TType>::instance().register_fn(                             \
                "true_" #name, bind<fn<TType>, TType, double>                              \
            );                                                                             \
        if constexpr((scope)==RegistrationScope::Reco || (scope)==RegistrationScope::Both) \
            VarFactoryRegistry<RType>::instance().register_fn(                             \
                "reco_" #name, bind<fn<RType>, RType, double>                              \
            );                                                                             \
        if constexpr((scope)==RegistrationScope::MCTruth)                                  \
            VarFactoryRegistry<MCTruth>::instance().register_fn(                           \
                "true_" #name, bind<fn<MCTruth>, MCTruth, double>                          \
            );                                                                             \
        return true;                                                                       \
    }();                                                                                   \
}

/**
 * @brief Operation mode for iteration over data products.
 * @details This enum class defines the operation mode for iteration over data
 * products. It can be either "true" or "reco". This is used in the
 * @ref construct function to determine the object to broadcast the selection
 * over.
 */
enum class Mode { True = 0, Reco = 1 };

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
 * @param mode The mode to use for the main loop ("true" or "reco").
 * @param override_type The type to use for the variable ("true" or "reco").
 * @return A NamedSpillMultiVar object that applies the cuts and computes the variable.
 * @throw std::runtime_error if a function is not registered.
 */
NamedSpillMultiVar construct(const std::vector<sys::cfg::ConfigurationTable> & cuts,
                             const sys::cfg::ConfigurationTable & var,
                             const std::string & mode,
                             const std::string & override_type = "");

/**
 * @brief Helper method for constructing a SpillMultiVar object.
 * @details This function is used to construct a SpillMultiVar object from
 * the parameter Cut and Variable objects. It is intended to be called by the
 * @ref construct function.
 * @tparam CutsOn The type (TType or RType) that the cut is applied to. This
 * also determines which type of object is iterated over in the loop. E.g.,
 * if CutsOn is TType, the loop iterates over the true events and applies
 * the cuts on truth information.
 * @tparam CompsOn The type (TType or RType) that is complementary to CutsOn
 * and is used to partition the events. This can be used, for example, to set
 * a truth cut on reco events passing the CutsOn cuts.
 * @tparam VarOn The type (TType or RType) that the variable is applied to.
 * @param cuts The callable that implements the cuts on the broadcast branch.
 * @param comps The callable that implements the cuts on the selected branch.
 * @param var The callable that implements the variable on the selected branch.
 * @return A SpillMultiVar object that applies the cuts and computes the variable.
 */
template<typename CutsOn, typename CompsOn, typename VarOn>
ana::SpillMultiVar spill_multivar_helper(
    const CutFn<CutsOn> & cuts,
    const CutFn<CompsOn> & comps,
    const VarFn<VarOn> & var
);

#endif // FRAMEWORK_H