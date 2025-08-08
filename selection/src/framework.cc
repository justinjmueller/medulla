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
NamedSpillMultiVar construct(const std::vector<cfg::ConfigurationTable> & cuts,
                             const cfg::ConfigurationTable & var,
                             const std::string & mode,
                             const std::string & override_type,
                             const bool ismc)
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
    else if(mode == "event") exec_mode = Mode::Event;
    else throw std::runtime_error("Illegal mode '" + mode + "' for variable " + var.get_string_field("name"));

    std::vector<CutFn<TType>> true_cut_functions;
    std::vector<CutFn<RType>> reco_cut_functions;
    std::vector<CutFn<TParticleType>> true_particle_cut_functions;
    std::vector<CutFn<RParticleType>> reco_particle_cut_functions;
    std::vector<CutFn<EventType>> event_cut_functions;
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
         
        if(!cut.has_field("type"))
            throw std::runtime_error("Cut " + name + " does not have a type field.");
        if(cut.get_string_field("type") == "true")
        {
            std::string cut_name = "true_" + name;
            std::vector<double> params;
            if(cut.has_field("parameters"))
                params = cut.get_double_vector("parameters");
            auto factory = CutFactoryRegistry<TType>::instance().get(cut_name);
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
        else if(cut.get_string_field("type") == "reco")
        {
            std::string cut_name = "reco_" + name;
            std::vector<double> params;
            if(cut.has_field("parameters"))
                params = cut.get_double_vector("parameters");
            auto factory = CutFactoryRegistry<RType>::instance().get(cut_name);
            if(invert)
            {
                // If the cut is inverted, we need to negate the function.
                auto fn = factory(params);
                reco_cut_functions.push_back([fn](const RType & e) { return !fn(e); });
            }
            else
                // Otherwise, we just add the function as is.
                reco_cut_functions.push_back(factory(params));
        }
        else if(cut.get_string_field("type") == "true_particle")
        {
            std::string cut_name = "true_particle_" + name;
            std::vector<double> params;
            if(cut.has_field("parameters"))
                params = cut.get_double_vector("parameters");
            auto factory = CutFactoryRegistry<TParticleType>::instance().get(cut_name);
            if(invert)
            {
                // If the cut is inverted, we need to negate the function.
                auto fn = factory(params);
                true_particle_cut_functions.push_back([fn](const TParticleType & e) { return !fn(e); });
            }
            else
                // Otherwise, we just add the function as is.
                true_particle_cut_functions.push_back(factory(params));
        }
        else if(cut.get_string_field("type") == "reco_particle")
        {
            std::string cut_name = "reco_particle_" + name;
            std::vector<double> params;
            if(cut.has_field("parameters"))
                params = cut.get_double_vector("parameters");
            auto factory = CutFactoryRegistry<RParticleType>::instance().get(cut_name);
            if(invert)
            {
                // If the cut is inverted, we need to negate the function.
                auto fn = factory(params);
                reco_particle_cut_functions.push_back([fn](const RParticleType & e) { return !fn(e); });
            }
            else
                // Otherwise, we just add the function as is.
                reco_particle_cut_functions.push_back(factory(params));
        }
        else if(cut.get_string_field("type") == "event")
        {
            std::string cut_name = "event_" + name;
            std::vector<double> params;
            if(cut.has_field("parameters"))
                params = cut.get_double_vector("parameters");
            auto factory = CutFactoryRegistry<EventType>::instance().get(cut_name);
            if(invert)
            {
                // If the cut is inverted, we need to negate the function.
                auto fn = factory(params);
                event_cut_functions.push_back([fn](const EventType & e) { return !fn(e); });
            }
            else
                // Otherwise, we just add the function as is.
                event_cut_functions.push_back(factory(params));
        }
        else if(cut.get_string_field("type") == "spill")
        {
            std::string cut_name = "spill_" + name;
            std::vector<double> params;
            if(cut.has_field("parameters"))
                params = cut.get_double_vector("parameters");
            auto factory = CutFactoryRegistry<SpillType>::instance().get(cut_name);

            // Transform this to a simple event-level cut.
            auto fn = [factory, params](const EventType & e) {
                if(!e.hdr.ismc)
                    return factory(params)(e.hdr.spillbnbinfo);
                else
                    return true; // If it's MC, we don't apply the spill cut.
            };
            if(invert)
            {
                // If the cut is inverted, we need to negate the function.
                event_cut_functions.push_back([fn](const EventType & e) {
                    if(!e.hdr.ismc)
                        return !fn(e);
                    else
                        return true; // If it's MC, we don't invert.
                });
            }
            else
                // Otherwise, we just add the function as is.
                event_cut_functions.push_back(fn);
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
    auto true_particle_cut = [true_particle_cut_functions](const TParticleType & e) -> bool {
        return std::all_of(true_particle_cut_functions.begin(), true_particle_cut_functions.end(), [&e](auto & f) { return f(e); });
    };
    auto reco_particle_cut = [reco_particle_cut_functions](const RParticleType & e) -> bool {
        return std::all_of(reco_particle_cut_functions.begin(), reco_particle_cut_functions.end(), [&e](auto & f) { return f(e); });
    };
    auto event_cut = [event_cut_functions](const EventType & e) -> bool {
        return std::all_of(event_cut_functions.begin(), event_cut_functions.end(), [&e](auto & f) { return f(e); });
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

        if(var_type == "true" || (var.has_field("selector") && var_type == "true_particle"))
        {
            if(var.has_field("selector"))
            {
                // Full name for the variable.
                std::string full_name = "true_" + var.get_string_field("selector") + "_" + var_name;

                // Retrieve the selector function.
                std::string selector_name = "true_" + var.get_string_field("selector");
                auto selector_factory = SelectorFactoryRegistry<TType>::instance().get(selector_name);
                auto selector = selector_factory(std::vector<double>{});
                                
                // Retrieve the particle-level variable function.
                var_name = "true_particle_" + var_name;
                auto factory = VarFactoryRegistry<TParticleType>::instance().get(var_name);
                auto var_fn = factory(varPars);
                
                VarFn<TType> var_fn_with_selector = [var_fn, selector](const TType & e) -> double
                {
                    // Apply the selector to the event.
                    size_t idx = selector(e);
                    if(idx == kNoMatch) return kNoMatchValue; // No match found.
                    // Apply the variable function to the selected particle.
                    return var_fn(e.particles[idx]);
                };
                return std::make_pair(full_name, spill_multivar_helper<TType, RType, TParticleType, TType>(
                    true_cut,
                    reco_cut_functions.empty() ? std::nullopt : std::optional<CutFn<RType>>(reco_cut),
                    true_particle_cut,
                    var_fn_with_selector,
                    event_cut,
                    ismc));
            }
            else
            {
                var_name = "true_" + var_name;
                auto factory = VarFactoryRegistry<TType>::instance().get(var_name);
                auto var_fn = factory(varPars);
                return std::make_pair(var_name, spill_multivar_helper<TType, RType, TParticleType, TType>(
                    true_cut,
                    reco_cut_functions.empty() ? std::nullopt : std::optional<CutFn<RType>>(reco_cut),
                    true_particle_cut,
                    var_fn,
                    event_cut,
                    ismc));
            }
        }
        
        else if(var_type == "reco" || (var.has_field("selector") && var_type == "reco_particle"))
        {
            if(var.has_field("selector"))
            {
                // Full name for the variable.
                std::string full_name = "reco_" + var.get_string_field("selector") + "_" + var_name;

                // Retrieve the selector function.
                std::string selector_name = "reco_" + var.get_string_field("selector");
                auto selector_factory = SelectorFactoryRegistry<RType>::instance().get(selector_name);
                auto selector = selector_factory(std::vector<double>{});
                                
                // Retrieve the particle-level variable function.
                var_name = "reco_particle_" + var_name;
                auto factory = VarFactoryRegistry<RParticleType>::instance().get(var_name);
                auto var_fn = factory(varPars);
                
                VarFn<RType> var_fn_with_selector = [var_fn, selector](const RType & e) -> double
                {
                    // Apply the selector to the event.
                    size_t idx = selector(e);
                    if(idx == kNoMatch) return kNoMatchValue; // No match found.
                    // Apply the variable function to the selected particle.
                    return var_fn(e.particles[idx]);
                };
                return std::make_pair(full_name, spill_multivar_helper<TType, RType, TParticleType, RType>(
                    true_cut,
                    reco_cut_functions.empty() ? std::nullopt : std::optional<CutFn<RType>>(reco_cut),
                    true_particle_cut,
                    var_fn_with_selector,
                    event_cut,
                    ismc));
            }
            else
            {
                var_name = "reco_" + var_name;
                auto factory = VarFactoryRegistry<RType>::instance().get(var_name);
                auto var_fn = factory(varPars);
                return std::make_pair(var_name, spill_multivar_helper<TType, RType, TParticleType, RType>(
                    true_cut,
                    reco_cut_functions.empty() ? std::nullopt : std::optional<CutFn<RType>>(reco_cut),
                    true_particle_cut,
                    var_fn,
                    event_cut,
                    ismc));
            }
        }
        else if(var_type == "mctruth")
        {
            var_name = "true_" + var_name;
            auto factory = VarFactoryRegistry<MCTruth>::instance().get(var_name);
            auto var_fn = factory(varPars);
            return std::make_pair(var_name, spill_multivar_helper<TType, RType, TParticleType, MCTruth>(
                true_cut,
                reco_cut_functions.empty() ? std::nullopt : std::optional<CutFn<RType>>(reco_cut),
                true_particle_cut,
                var_fn,
                event_cut,
                ismc));
        }
        else if(var_type == "true_particle")
        {
            var_name = "true_particle_" + var_name;
            auto factory = VarFactoryRegistry<TParticleType>::instance().get(var_name);
            auto var_fn = factory(varPars);
            return std::make_pair(var_name, spill_multivar_helper<TType, RType, TParticleType, TParticleType>(
                true_cut,
                reco_cut_functions.empty() ? std::nullopt : std::optional<CutFn<RType>>(reco_cut),
                true_particle_cut,
                var_fn,
                event_cut,
                ismc));
        }
        else if(var_type == "reco_particle")
        {
            var_name = "reco_particle_" + var_name;
            auto factory = VarFactoryRegistry<RParticleType>::instance().get(var_name);
            auto var_fn = factory(varPars);
            return std::make_pair(var_name, spill_multivar_helper<TType, RType, TParticleType, RParticleType>(
                true_cut,
                reco_cut_functions.empty() ? std::nullopt : std::optional<CutFn<RType>>(reco_cut),
                true_particle_cut,
                var_fn,
                event_cut,
                ismc));
        }
        else
        {
            throw std::runtime_error("Illegal variable type '" + var_type + "' for variable " + var_name);
        }
    }
    else if(exec_mode == Mode::Reco)
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

        if(var_type == "true" || (var.has_field("selector") && var_type == "true_particle"))
        {
            if(var.has_field("selector"))
            {
                // Full name for the variable.
                std::string full_name = "true_" + var.get_string_field("selector") + "_" + var_name;

                // Retrieve the selector function.
                std::string selector_name = "true_" + var.get_string_field("selector");
                auto selector_factory = SelectorFactoryRegistry<TType>::instance().get(selector_name);
                auto selector = selector_factory(std::vector<double>{});
                                
                // Retrieve the particle-level variable function.
                var_name = "true_particle_" + var_name;
                auto factory = VarFactoryRegistry<TParticleType>::instance().get(var_name);
                auto var_fn = factory(varPars);
                
                VarFn<TType> var_fn_with_selector = [var_fn, selector](const TType & e) -> double
                {
                    // Apply the selector to the event.
                    size_t idx = selector(e);
                    if(idx == kNoMatch) return kNoMatchValue; // No match found.
                    // Apply the variable function to the selected particle.
                    return var_fn(e.particles[idx]);
                };
                return std::make_pair(full_name, spill_multivar_helper<RType, TType, TParticleType, TType>(
                    reco_cut,
                    true_cut_functions.empty() ? std::nullopt : std::optional<CutFn<TType>>(true_cut),
                    true_particle_cut,
                    var_fn_with_selector,
                    event_cut,
                    ismc));
            }
            else
            {
                var_name = "true_" + var_name;
                auto factory = VarFactoryRegistry<TType>::instance().get(var_name);
                auto var_fn = factory(varPars);
                return std::make_pair(var_name, spill_multivar_helper<RType, TType, TParticleType, TType>(
                    reco_cut,
                    true_cut_functions.empty() ? std::nullopt : std::optional<CutFn<TType>>(true_cut),
                    true_particle_cut,
                    var_fn,
                    event_cut,
                    ismc));
            }
        }
        else if(var_type == "reco" || (var.has_field("selector") && var_type == "reco_particle"))
        {
            if(var.has_field("selector"))
            {
                // Full name for the variable.
                std::string full_name = "reco_" + var.get_string_field("selector") + "_" + var_name;

                // Retrieve the selector function.
                std::string selector_name = "reco_" + var.get_string_field("selector");
                auto selector_factory = SelectorFactoryRegistry<RType>::instance().get(selector_name);
                auto selector = selector_factory(std::vector<double>{});

                // Retrieve the particle-level variable function.
                var_name = "reco_particle_" + var_name;
                auto factory = VarFactoryRegistry<RParticleType>::instance().get(var_name);
                auto var_fn = factory(varPars);

                VarFn<RType> var_fn_with_selector = [var_fn, selector](const RType & e) -> double
                {
                    // Apply the selector to the event.
                    size_t idx = selector(e);
                    if(idx == kNoMatch) return kNoMatchValue; // No match found.
                    // Apply the variable function to the selected particle.
                    return var_fn(e.particles[idx]);
                };
                return std::make_pair(full_name, spill_multivar_helper<RType, TType, RParticleType, RType>(
                    reco_cut,
                    true_cut_functions.empty() ? std::nullopt : std::optional<CutFn<TType>>(true_cut),
                    reco_particle_cut,
                    var_fn_with_selector,
                    event_cut,
                    ismc));
            }
            else
            {
                var_name = "reco_" + var_name;
                auto factory = VarFactoryRegistry<RType>::instance().get(var_name);
                auto var_fn = factory(varPars);
                return std::make_pair(var_name, spill_multivar_helper<RType, TType, TParticleType, RType>(
                    reco_cut,
                    true_cut_functions.empty() ? std::nullopt : std::optional<CutFn<TType>>(true_cut),
                    true_particle_cut,
                    var_fn,
                    event_cut,
                    ismc));
            }
        }
        else if(var_type == "mctruth")
        {
            var_name = "true_" + var_name;
            auto factory = VarFactoryRegistry<MCTruth>::instance().get(var_name);
            auto var_fn = factory(varPars);
            return std::make_pair(var_name, spill_multivar_helper<RType, TType, TParticleType, MCTruth>(
                reco_cut,
                true_cut_functions.empty() ? std::nullopt : std::optional<CutFn<TType>>(true_cut),
                true_particle_cut,
                var_fn,
                event_cut,
                ismc));
        }
        else if(var_type == "true_particle")
        {
            var_name = "true_particle_" + var_name;
            auto factory = VarFactoryRegistry<TParticleType>::instance().get(var_name);
            auto var_fn = factory(varPars);
            return std::make_pair(var_name, spill_multivar_helper<RType, TType, RParticleType, TParticleType>(
                reco_cut,
                true_cut_functions.empty() ? std::nullopt : std::optional<CutFn<TType>>(true_cut),
                reco_particle_cut,
                var_fn,
                event_cut,
                ismc));
        }
        else if(var_type == "reco_particle")
        {
            var_name = "reco_particle_" + var_name;
            auto factory = VarFactoryRegistry<RParticleType>::instance().get(var_name);
            auto var_fn = factory(varPars);
            return std::make_pair(var_name, spill_multivar_helper<RType, TType, RParticleType, RParticleType>(
                reco_cut,
                true_cut_functions.empty() ? std::nullopt : std::optional<CutFn<TType>>(true_cut),
                reco_particle_cut,
                var_fn,
                event_cut,
                ismc));
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

        if(var_type == "event")
        {
            var_name = "event_" + var_name;
            auto factory = VarFactoryRegistry<EventType>::instance().get(var_name);
            auto var_fn = factory(varPars);
            return std::make_pair(var_name, spill_multivar_helper(event_cut, var_fn));
        }
        else
        {
            throw std::runtime_error("Illegal variable type '" + var_type + "' for variable " + var_name);
        }
    }
}

// Helper method for constructing a SpillMultiVar object.
template<typename CutsOn, typename CompsOn, typename PCutsOn, typename VarOn>
ana::SpillMultiVar spill_multivar_helper(
    const CutFn<CutsOn> & cuts,
    const std::optional<CutFn<CompsOn>> & comps,
    const CutFn<PCutsOn> & pcuts,
    const VarFn<VarOn> & var,
    const CutFn<EventType> & event_cut,
    const bool ismc
)
{
    return ana::SpillMultiVar([comps, cuts, pcuts, var, ismc, event_cut](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>
    {
        std::vector<double> values;

        // Check if this event passes the event cut.
        if(!event_cut(*sr)) return values;

        // Case: configuration parameter "mode" is set to "true."
        if constexpr (std::is_same_v<CutsOn, TType>)
        {
            // Case: the variable type is a particle type.
            // If the variable is of the "particle" type, we need to build a
            // lookup table of matches.
            std::map<caf::Proxy<int64_t>, const caf::Proxy<caf::SRParticleDLP> *> particles;
            if constexpr (std::is_same_v<VarOn, TParticleType> || std::is_same_v<VarOn, RParticleType>)
            {
                for(auto const& i : sr->dlp)
                {
                    for(auto const& j : i.particles)
                    {
                        particles.insert(std::make_pair(j.id, &j));
                    }
                }
            }

            // Iterate over the true interactions.
            for(auto const& i : sr->dlp_true)
            {
                // Check for match
                size_t match_id = (i.match_ids.size() > 0) ? (size_t)i.match_ids[0] : kNoMatch;

                if constexpr(std::is_same_v<VarOn, RType>)
                {
                    if(cuts(i) && (!comps || (match_id != kNoMatch && (*comps)(sr->dlp[match_id]))))                    
                    {
                        values.push_back(match_id != kNoMatch ? var(sr->dlp[match_id]) : kNoMatchValue);
                    }
                }
                else if constexpr(std::is_same_v<VarOn, TType>)
                {
                    if(cuts(i) && (!comps || (match_id != kNoMatch && (*comps)(sr->dlp[match_id]))))
                    {
                        values.push_back(var(i));
                    }
                }
                else if constexpr(std::is_same_v<VarOn, MCTruth>)
                {
                    if(cuts(i) && (!comps || (match_id != kNoMatch && (*comps)(sr->dlp[match_id]))))
                    {
                        values.push_back(i.nu_id >= 0 ? var(sr->mc.nu[i.nu_id]) : kNoMatchValue);
                    }
                }
                else if constexpr(std::is_same_v<VarOn, TParticleType> || std::is_same_v<VarOn, RParticleType>)
                {
                    if(cuts(i) && (!comps || (match_id != kNoMatch && (*comps)(sr->dlp[match_id]))))
                    {
                        for(auto const & j : i.particles)
                        {
                            if(pcuts(j))
                            {
                                if constexpr(std::is_same_v<VarOn, TParticleType>)
                                    values.push_back(var(j));
                                else if constexpr(std::is_same_v<VarOn, RParticleType>)
                                {
                                    if(j.match_ids.size() > 0 && particles.find(j.match_ids[0]) != particles.end())
                                        values.push_back(var(*particles[j.match_ids[0]]));
                                    else
                                        values.push_back(kNoMatchValue); // No match found.
                                }
                            }
                        }
                    }
                }
            }
        }

        // Case: configuration parameter "mode" is set to "reco."
        else if constexpr(std::is_same_v<CutsOn, RType>)
        {
            // Case: the variable type is a particle type.
            // If the variable is of the "particle" type, we need to build a
            // lookup table of matches.
            std::map<caf::Proxy<int64_t>, const caf::Proxy<caf::SRParticleTruthDLP> *> particles;
            if constexpr (std::is_same_v<VarOn, TParticleType> || std::is_same_v<VarOn, RParticleType>)
            {
                for(auto const& i : sr->dlp_true)
                {
                    for(auto const& j : i.particles)
                    {
                        particles.insert(std::make_pair(j.id, &j));
                    }
                }
            }

            // Iterate over the reco interactions.
            for(auto const& i : sr->dlp)
            {
                // Check for match
                size_t match_id = (i.match_ids.size() > 0) ? (size_t)i.match_ids[0] : kNoMatch;

                if constexpr(std::is_same_v<VarOn, TType>)
                {
                    if(cuts(i) && (!comps || (match_id != kNoMatch && (*comps)(sr->dlp_true[match_id])) || !ismc))
                    {
                        values.push_back(ismc && match_id != kNoMatch ? var(sr->dlp_true[match_id]) : kNoMatchValue);
                    }
                }
                else if constexpr(std::is_same_v<VarOn, RType>)
                {
                    if(cuts(i) && (!comps || (match_id != kNoMatch && (*comps)(sr->dlp_true[match_id])) || !ismc))
                    {
                        values.push_back(var(i));
                    }
                }
                else if constexpr(std::is_same_v<VarOn, MCTruth>)
                {
                    if(cuts(i) && (!comps || (match_id != kNoMatch && (*comps)(sr->dlp_true[match_id])) || !ismc))
                    {
                        if(!ismc || match_id == kNoMatch)
                        {
                            values.push_back(kNoMatchValue);
                        }
                        else
                        {
                            int64_t nu_id = sr->dlp_true[match_id].nu_id;
                            values.push_back(nu_id >= 0 ? var(sr->mc.nu[nu_id]) : kNoMatchValue);
                        }
                    }
                }
                else if constexpr(std::is_same_v<VarOn, TParticleType> || std::is_same_v<VarOn, RParticleType>)
                {
                    if(cuts(i) && (!comps || (match_id != kNoMatch && (*comps)(sr->dlp_true[match_id]))))
                    {
                        for(auto const & j : i.particles)
                        {
                            if(pcuts(j))
                            {
                                if constexpr(std::is_same_v<VarOn, RParticleType>)
                                    values.push_back(var(j));
                                else if constexpr(std::is_same_v<VarOn, TParticleType>)
                                {
                                    if(j.match_ids.size() > 0 && particles.find(j.match_ids[0]) != particles.end())
                                        values.push_back(var(*particles[j.match_ids[0]]));
                                    else
                                        values.push_back(kNoMatchValue); // No match found.
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // Return the collected values.
        return values;
    });
}

// Helper method for constructing a SpillMultiVar object when run in the
// "event" mode.
ana::SpillMultiVar spill_multivar_helper(const CutFn<EventType> & cut, const VarFn<EventType> & var)
{
    return ana::SpillMultiVar([cut, var](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>
    {
        std::vector<double> values;
        if(cut(*sr))        
            values.push_back(var(*sr));
        return values;
    });
}

// Helper method for constructing a set of SpillMultiVar objects that track the
// exposure information for a given set of cuts.
std::vector<NamedSpillMultiVar> construct_exposure_vars(const std::vector<cfg::ConfigurationTable> & cuts)
{
    std::vector<NamedSpillMultiVar> exposure_vars;
    std::vector<CutFn<EventType>> cut_functions;
    std::vector<CutFn<SpillType>> spill_cut_functions;

    // Iterate over the cuts and construct the exposure variables.
    for(const auto & cut : cuts)
    {
        // Check if the cut decrements the exposure.
        if(cut.get_bool_field("decrements_exposure", false))
        {
            // Retrieve the cut parameters.
            std::vector<double> params;
            if(cut.has_field("parameters"))
                params = cut.get_double_vector("parameters");

            // Retrieve the cut function.
            std::string name = cut.get_string_field("name");

            if(cut.get_string_field("type") == "event")
            {
                name = "event_" + name;
                auto factory = CutFactoryRegistry<EventType>::instance().get(name);
                auto cut_fn = factory(params);
                cut_functions.push_back(cut_fn);
            }
            else if(cut.get_string_field("type") == "spill")
            {
                name = "spill_" + name;
                auto factory = CutFactoryRegistry<SpillType>::instance().get(name);
                auto cut_fn = factory(params);
                
                // We do not transform this to an event-level cut because we
                // need to apply it to each and every spill that contains
                // exposure that we want to track.
                spill_cut_functions.push_back(cut_fn);
            }
            else
            {
                throw std::runtime_error("Illegal cut type '" + cut.get_string_field("type") + "' for exposure cut " + name);
            }
        }
    }

    // Compose a common cut function.
    auto cut = [cut_functions](const EventType & e) -> bool {
        return std::all_of(cut_functions.begin(), cut_functions.end(), [&e](auto & f) { return f(e); });
    };
    
    // Compose a common spill cut function.
    auto spill_cut = [spill_cut_functions](const SpillType & s) -> bool {
        return std::all_of(spill_cut_functions.begin(), spill_cut_functions.end(), [&s](auto & f) { return f(s); });
    };

    // Compose the exposure variables
    auto livetime_var = [](const EventType & e) -> double {
        // Return the livetime for the event.
        if(e.hdr.ismc)
            return (e.hdr.first_in_subrun) ? (double)e.hdr.ngenevt : 0.0;
        else
            return e.hdr.bnbinfo.size() + e.hdr.numiinfo.size() + e.hdr.noffbeambnb + e.hdr.noffbeamnumi;
    };
    exposure_vars.push_back(std::make_pair("livetime", spill_multivar_helper(cut, livetime_var)));

    auto pot_var = [spill_cut](const EventType & e) -> double {
        
        if(e.hdr.ismc)
            return (e.hdr.first_in_subrun) ? (double)e.hdr.pot : 0.0;
        else
        {
            double tot(0);
            for(const auto & bnb : e.hdr.bnbinfo)
                tot += (spill_cut(bnb) ? (double)bnb.TOR875 : 0);
            return tot;
        }
    };
    exposure_vars.push_back(std::make_pair("pot", spill_multivar_helper(cut, pot_var)));

    return exposure_vars;
}

// Explicitly instantiate Registry for the factory types we use:
// Cut Registry
template class Registry<CutFactory<TType>>;
template class Registry<CutFactory<RType>>;
template class Registry<CutFactory<TParticleType>>;
template class Registry<CutFactory<RParticleType>>;
template class Registry<CutFactory<EventType>>;
template class Registry<CutFactory<SpillType>>;

// Var Registry
template class Registry<VarFactory<TType>>;
template class Registry<VarFactory<RType>>;
template class Registry<VarFactory<MCTruth>>;
template class Registry<VarFactory<TParticleType>>;
template class Registry<VarFactory<RParticleType>>;
template class Registry<VarFactory<EventType>>;

// Explicit instantiation for selector registries
template class Registry<SelectorFactory<TType>>;
template class Registry<SelectorFactory<RType>>;