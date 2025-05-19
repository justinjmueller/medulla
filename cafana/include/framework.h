#ifndef FRAMEWORK_H
#define FRAMEWORK_H
#include <map>
#include <string>
#include <functional>
#include <stdexcept>
#include <type_traits>

#include "sbnana/CAFAna/Core/MultiVar.h"
#include "configuration.h"

using TType = caf::SRInteractionTruthDLPProxy;
using RType = caf::SRInteractionDLPProxy;
using NamedSpillMultiVar = std::pair<std::string, ana::SpillMultiVar>;

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
         * @return A copy of the registered Fn. 
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

//===----------------------------------------------------------------------===//
// 2) Factory‐based registries: bind params → callable<EventT>
//===----------------------------------------------------------------------===//

/// @brief A factory function: given params, returns a bool‑cut on EventT.
template<typename EventT>
using CutFactory = std::function<std::function<bool(const EventT&)>(const std::vector<double>&)>;

/// @brief Registry of CutFactory<EventT> by name.
template<typename EventT>
class CutFactoryRegistry {
public:
  static CutFactoryRegistry& Instance() {
    static CutFactoryRegistry inst;
    return inst;
  }

  /// Register a factory that binds parameters
  void Register(std::string name, CutFactory<EventT> f) {
    registry_[std::move(name)] = std::move(f);
  }

  /// Create a bound cut by name + params
  std::function<bool(const EventT&)> Create(const std::string& name,
                                             const std::vector<double>& pars) const
  {
    auto it = registry_.find(name);
    if(it == registry_.end())
      throw std::runtime_error("Unknown cut factory: " + name);
    return it->second(pars);
  }

private:
  std::map<std::string, CutFactory<EventT>> registry_;
};

/// @brief Registry of variable factory functions (similar to CutFactoryRegistry).
template<typename EventT>
using VarFactory = std::function<std::function<double(const EventT&)>(const std::vector<double>&)>;

template<typename EventT>
class VarFactoryRegistry {
public:
  static VarFactoryRegistry& Instance() {
    static VarFactoryRegistry inst;
    return inst;
  }

  void Register(std::string name, VarFactory<EventT> f) {
    registry_[std::move(name)] = std::move(f);
  }

  std::function<double(const EventT&)> Create(const std::string& name,
                                             const std::vector<double>& pars) const
  {
    auto it = registry_.find(name);
    if(it == registry_.end())
      throw std::runtime_error("Unknown variable factory: " + name);
    return it->second(pars);
  }

private:
  std::map<std::string, VarFactory<EventT>> registry_;
};

#include <type_traits>
#include <vector>

namespace sys { namespace fw {

// Overload‐detecting binder for cuts:
template<auto F, typename EventT>
std::function<bool(const EventT&)>
BindCut(const std::vector<double>& pars) {
  if constexpr (std::is_invocable_v<decltype(F), const EventT&, const std::vector<double>&>) {
    return [pars](const EventT& e){ return F(e, pars); };
  } else {
    return [=](const EventT& e){ return F(e); };
  }
}

// Overload‐detecting binder for variables:
template<auto F, typename EventT>
std::function<double(const EventT&)>
BindVar(const std::vector<double>& pars) {
  if constexpr (std::is_invocable_v<decltype(F), const EventT&, const std::vector<double>&>) {
    return [pars](const EventT& e){ return F(e, pars); };
  } else {
    return [=](const EventT& e){ return F(e); };
  }
}

}} // namespace sys::fw

/// @brief Scope for registration macros
enum class RegistrationScope { True, Reco, Both };

// Register a cut with scope, auto‐detecting its signature
#define REGISTER_CUT_SCOPE(scope, name, fn)                                            \
namespace {                                                                            \
  const bool _reg_cut_##name = []{                                                     \
    if constexpr((scope)==RegistrationScope::True || (scope)==RegistrationScope::Both) \
      CutFactoryRegistry<TType>::Instance().Register(                                  \
        "true_" #name, sys::fw::BindCut<+fn<TType>, TType>                             \
      );                                                                               \
    if constexpr((scope)==RegistrationScope::Reco || (scope)==RegistrationScope::Both) \
      CutFactoryRegistry<RType>::Instance().Register(                                  \
        "reco_" #name, sys::fw::BindCut<+fn<RType>, RType>                             \
      );                                                                               \
    return true;                                                                       \
  }();                                                                                 \
}

// Register a variable with scope, auto‐detecting its signature
#define REGISTER_VAR_SCOPE(scope, name, fn)                                            \
namespace {                                                                            \
  const bool _reg_var_##name = []{                                                     \
    if constexpr((scope)==RegistrationScope::True || (scope)==RegistrationScope::Both) \
      VarFactoryRegistry<TType>::Instance().Register(                                  \
        "true_" #name, sys::fw::BindVar<+fn<TType>, TType>                             \
      );                                                                               \
    if constexpr((scope)==RegistrationScope::Reco || (scope)==RegistrationScope::Both) \
      VarFactoryRegistry<RType>::Instance().Register(                                  \
        "reco_" #name, sys::fw::BindVar<+fn<RType>, RType>                             \
      );                                                                               \
    return true;                                                                       \
  }();                                                                                 \
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
 * @return A NamedSpillMultiVar object that applies the cuts and computes the variable.
 * @throw std::runtime_error if a function is not registered.
 */
NamedSpillMultiVar construct(const std::vector<sys::cfg::ConfigurationTable> & cuts, const sys::cfg::ConfigurationTable & var, const std::string & override_type = "");

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