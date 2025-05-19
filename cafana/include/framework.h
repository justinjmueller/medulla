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

/// @brief Register a cut that takes (EventT const&, std::vector<double> const&)
#define REGISTER_CUT_PARAMS(name, fn)                                           \
namespace {                                                                     \
  struct Registrar_CutParams_##name {                                           \
    Registrar_CutParams_##name() {                                              \
      CutFactoryRegistry<TType>::Instance().Register(                           \
        "true_" #name,                                                          \
        [](const std::vector<double>& pars) -> std::function<bool(const TType&)> { \
          return [pars](const TType& ev){ return fn<TType>(ev, pars); };        \
        }                                                                       \
      );                                                                        \
      CutFactoryRegistry<RType>::Instance().Register(                           \
        "reco_" #name,                                                          \
        [](const std::vector<double>& pars) -> std::function<bool(const RType&)> { \
          return [pars](const RType& ev){ return fn<RType>(ev, pars); };        \
        }                                                                       \
      );                                                                        \
    }                                                                           \
  } registrar_CutParams_##name;                                                 \
}

/// @brief Register a cut that takes just (EventT const&)
#define REGISTER_CUT_NO_PARAMS(name, fn)                                        \
namespace {                                                                     \
  struct Registrar_CutNP_##name {                                               \
    Registrar_CutNP_##name() {                                                  \
      CutFactoryRegistry<TType>::Instance().Register(                           \
        "true_" #name,                                                          \
        [](const std::vector<double>&) -> std::function<bool(const TType&)> {   \
          return fn<TType>;                                                     \
        }                                                                       \
      );                                                                        \
      CutFactoryRegistry<RType>::Instance().Register(                           \
        "reco_" #name,                                                          \
        [](const std::vector<double>&) -> std::function<bool(const RType&)> {   \
          return fn<RType>;                                                     \
        }                                                                       \
      );                                                                        \
    }                                                                           \
  } registrar_CutNP_##name;                                                     \
}


/// @brief Register a variable that takes (EventT const&, std::vector<double> const&)
#define REGISTER_VAR_PARAMS(name, fn)                                                \
namespace {                                                                          \
  struct Registrar_VarParams_##name {                                                \
    Registrar_VarParams_##name() {                                                   \
      VarFactoryRegistry<TType>::Instance().Register(                                \
        "true_" #name,                                                               \
        [](const std::vector<double>& pars) -> std::function<double(const TType&)> { \
          return [pars](const TType& ev){ return fn<TType>(ev, pars); };             \
        }                                                                            \
      );                                                                             \
      VarFactoryRegistry<RType>::Instance().Register(                                \
        "reco_" #name,                                                               \
        [](const std::vector<double>& pars) -> std::function<double(const RType&)> { \
          return [pars](const RType& ev){ return fn<RType>(ev, pars); };             \
        }                                                                            \
      );                                                                             \
    }                                                                                \
  } registrar_VarParams_##name;                                                      \
}

/// @brief Register a variable that takes just (EventT const&)
#define REGISTER_VAR_NO_PARAMS(name, fn)                                             \
namespace {                                                                          \
  struct Registrar_VarNP_##name {                                                    \
    Registrar_VarNP_##name() {                                                       \
      VarFactoryRegistry<TType>::Instance().Register(                                \
        "true_" #name,                                                               \
        [](const std::vector<double>&) -> std::function<double(const TType&)> {      \
          return fn<TType>;                                                          \
        }                                                                            \
      );                                                                             \
      VarFactoryRegistry<RType>::Instance().Register(                                \
        "reco_" #name,                                                               \
        [](const std::vector<double>&) -> std::function<double(const RType&)> {      \
          return fn<RType>;                                                          \
        }                                                                            \
      );                                                                             \
    }                                                                                \
  } registrar_VarNP_##name;                                                          \
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