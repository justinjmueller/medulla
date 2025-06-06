/**
 * @file utilities_electron2025.h
 * @brief Header file for definitions of utility functions for supporting
 * analysis variables and cuts.
 * @details This file contains definitions of utility functions which are used
 * to support the implementation of analysis variables and cuts. These functions
 * are intended to be used to simplify the implementation of variables and cuts
 * by providing common functionality which can be reused across multiple
 * variables and cuts.
 * @author jzettle@fnal.gov/mueller@fnal.gov
 */
#ifndef UTILITIES_ELECTRON2025_H
#define UTILITIES_ELECTRON2025_H

#include <vector>

#include "include/particle_variables.h"
#include "include/particle_cuts.h"

/**
 * @namespace utilities
 * @brief Namespace for organizing utility functions for supporting analysis
 * variables and cuts.
 * @details This namespace is intended to be used for organizing utility
 * functions which are used to support the implementation of analysis variables
 * and cuts. These functions are intended to be used to simplify the
 * implementation of variables and cuts by providing common functionality which
 * can be reused across multiple variables and cuts.
 * @note The namespace is intended to be used in conjunction with the
 * vars and cuts namespaces, which are used for organizing variables and cuts
 * which act on interactions.
 */
namespace utilities::electron2025
{
    /**
     * @brief Count the primaries of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to find the topology of.
     * @return the count of primaries of each particle type within the
     * interaction.
     */
    template<class T>
        std::vector<uint32_t> count_primaries(const T & obj)
        {
            std::vector<uint32_t> counts(5, 0);
            for(auto &p : obj.particles)
            {
                if(pcuts::final_state_signal(p))
                    ++counts[p.pid];
            }
            return counts;
        }

     template<class T>
        std::vector<uint32_t> count_primaries_ee(const T & obj)
        {
            std::vector<uint32_t> counts(5, 0);
            for(auto &p : obj.particles)
            {
                if(pcuts::final_state_signal_elec(p))
                    ++counts[p.pid];
            }
            return counts;
        }
    
    /**
     * @brief Finds the index corresponding to the leading particle of the specifed
     * particle type.
     * @details The leading particle is defined as the particle with the highest
     * kinetic energy. If the interaction is a true interaction, the initial kinetic
     * energy is used instead of the CSDA kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @param pid of the particle type.
     * @return the index of the leading particle (highest KE). 
     */
    template <class T>
        size_t leading_particle_index(const T & obj, uint16_t pid)
        {
            double leading_ke(0);
            size_t index(0);
            for(size_t i(0); i < obj.particles.size(); ++i)
            {
                const auto & p = obj.particles[i];
                double energy(p.csda_ke);
                if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                    energy = pvars::ke(p);
                if(p.pid == pid && energy > leading_ke)
                {
                    leading_ke = energy;
                    index = i;
                }
            }
            return index;
        }

    template <class T>
        size_t leading_particle_index_shower(const T & obj, uint16_t pid1, uint16_t pid2)
        {
            double leading_ke(0);
            size_t index(0);
            for(size_t i(0); i < obj.particles.size(); ++i)
            {
                const auto & p = obj.particles[i];
                double energy(p.csda_ke);
                if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                    energy = pvars::ke(p);
                if((p.pid == pid1 || p.pid == pid2) && energy > leading_ke && energy > 50)
                {
                    leading_ke = energy;
                    index = i;
                }
            }
            return index;
        }

    template <class T>
        std::vector<size_t> particle_indices(const T & obj, uint16_t pid1, uint16_t pid2)
        {
            //enable selection of any shower type (electron or photon) at the moment for benchmarking purposes
            double leading_ke(0);
            double subleading_ke(0);
            size_t index(0);
            size_t subindex(0);
            std::vector<size_t> indices;
            //if there is only one particle in the list (I expect this is the case for reconstructed only one shower), return 0 for both indices
            //figure out a smarter way, still has some that pass this and are added to the other plots, must not always be the case even in reco (are there tracks?)
            /*
            if(obj.particles.size() == 1)
            {
                indices.push_back(0);
                return indices;
            }
            */
            //if only one shower but possibly other particles, try this, otherwise it will return 0 for the second index even if only one shower
            /*
            std::vector<uint32_t> c(count_primaries_ee(obj));
            if((c[0] == 1 && c[1] == 0) || (c[0] == 0 && c[1] == 1))
            {
                indices.push_back(0);
                return indices;
            }
	         */
            //find leading particle of a given type
            for(size_t i(0); i < obj.particles.size(); ++i)
            {
                const auto & p = obj.particles[i];
                double energy(p.calo_ke);
                if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                    energy = pvars::ke(p);
                
                if((p.pid == pid1 || p.pid == pid2) && (energy > leading_ke && p.is_primary && energy > 50))
                {
                    leading_ke = energy;
                    index = i;
                }
            }
            
            indices.push_back(index);
            //find subleading particle of any shower type (electron, photon)
            for(size_t i(0); i < obj.particles.size(); ++i)
            {
                const auto & p = obj.particles[i];
                double energy(p.calo_ke);
                if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                    energy = pvars::ke(p);
                if((p.pid == pid1 || p.pid == pid2) && (energy > subleading_ke && p.is_primary && energy > 50) && i != index)
                {
                    subleading_ke = energy;
                    subindex = i;
                }
            }
           
            indices.push_back(subindex);

            return indices;
        }

    /**
     * @brief Finds the index corresponding to the leading muon.
     * @details The leading muon is defined as the muon with the highest
     * kinetic energy. If the interaction is a true interaction, the initial
     * kinetic energy is used instead of the CSDA kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading muon (highest KE).
     */
    template<class T>
        size_t leading_muon_index(const T & obj)
        {
            return leading_particle_index(obj, 2);
        }
    
    /**
     * @brief Finds the index corresponding to the leading proton.
     * @details The leading proton is defined as the proton with the highest
     * kinetic energy. If the interaction is a true interaction, the initial
     * kinetic energy is used instead of the CSDA kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading proton (highest KE).
     */
    template<class T>
        size_t leading_proton_index(const T & obj)
        {
            return leading_particle_index(obj, 4);
        }

    template<class T>
        size_t leading_shower_index(const T & obj)
        {
            return leading_particle_index_shower(obj, 1, 0);
        }
    template<class T>
        size_t subleading_shower_index(const T & obj)
        {
            std::vector<size_t> indices = utilities::electron2025::particle_indices(obj, 1, 0);
            return indices[1];
        }
}
#endif // UTILITIES_ELECTRON2025_H
