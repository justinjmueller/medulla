/**
 * @file selectors.h
 * @brief Header file for the selectors used in the SPINE analysis framework.
 * @details This file contains the definitions of selectors which can be used
 * to select a single particle (by index) within an interaction. This is a
 * useful feature for reducing down a collection of particles in a final state
 * to just a single one, which allows a user to broadcast a particle-level
 * variable "upwards" to the interaction level (i.e., it can be placed in a
 * branch of a tree that is otherwise filled with interaction-level variables).
 * @author mueller@fnal.gov
 */
#ifndef SELECTORS_H
#define SELECTORS_H
#include <vector>

#include "framework.h"
#include "include/particle_cuts.h"
#include "include/particle_variables.h"

/**
 * @namespace selectors
 * @brief Namespace for organizing selectors which act on interactions.
 * @details This namespace is intended to be used for organizing selectors
 * which act on interactions. Each selector is implemented as a function
 * which takes an interaction object as an argument and returns the index
 * of the selected particle. The function should be templated on the type
 * of interaction object if the selector is intended to be used on both
 * true and reconstructed interactions.
 */
namespace selectors
{
    /**
     * @brief Finds the index corresponding to the leading particle of the
     * specified particle type.
     * @details The leading particle is defined as the particle with the
     * highest kinetic energy. The method of calculating kinetic energy is
     * inherited by the @ref pvars::ke function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @param pid of the particle type.
     * @return the index of the leading particle (highest KE). 
     */
    template <class T>
    size_t leading_particle_index(const T & obj, uint16_t pid)
    {
        double leading_ke(0);
        size_t index(kNoMatch);
        for(size_t i(0); i < obj.particles.size(); ++i)
        {
            const auto & p = obj.particles[i];
            double energy(pvars::ke(p));
            if(pvars::pid(p) == pid && energy > leading_ke)
            {
                leading_ke = energy;
                index = i;
            }
        }
        return index;
    }

    /**
     * @brief Finds the index corresponding to the longest track.
     * @details The longest track is defined as the track with the longest
     * length, which is calculated upstream in SPINE. The particle instance is
     * required to have a semantic type of 1 (track) and have a start point
     * within 6 cm of the interaction vertex.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the longest track.
     */
    template<class T>
    size_t longest_track(const T & obj)
    {
        double longest_length(0);
        size_t index(kNoMatch);
        for(size_t i(0); i < obj.particles.size(); ++i)
        {
            const auto & p = obj.particles[i];

            // Distance between interaction vertex and particle start.
            double vertex_distance = std::sqrt(
                std::pow(pvars::start_x(p) - obj.vertex[0], 2) +
                std::pow(pvars::start_y(p) - obj.vertex[1], 2) +
                std::pow(pvars::start_z(p) - obj.vertex[2], 2)
            );

            // Skip particles that are not tracks or are too far from the
            // interaction vertex.
            if(pvars::semantic_type(p) != 1 || vertex_distance >= 6)
                continue;
            
            // Update the longest length and index if the current particle
            // is longer than the longest found so far.
            if(pvars::length(p) > longest_length)
            {
                longest_length = pvars::length(p);
                index = i;
            }
        }
        return index;
    }
    REGISTER_SELECTOR(longest_track, longest_track);

    /**
     * @brief Finds the index corresponding to the second longest track.
     * @details The second longest track is defined as the track with the
     * second longest length, which is calculated upstream in SPINE. The
     * particle instance is required to have a semantic type of 1 (track) and
     * have a start point within 6 cm of the interaction vertex.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the second longest track.
    */
    template<class T>
    size_t second_longest_track(const T & obj)
    {
        double longest_length(0);
        double second_longest_length(0);
        size_t index(kNoMatch), second_index(kNoMatch);
        for(size_t i(0); i < obj.particles.size(); ++i)
        {
            const auto & p = obj.particles[i];

            // Distance between interaction vertex and particle start.
            double vertex_distance = std::sqrt(
                std::pow(pvars::start_x(p) - obj.vertex[0], 2) +
                std::pow(pvars::start_y(p) - obj.vertex[1], 2) +
                std::pow(pvars::start_z(p) - obj.vertex[2], 2)
            );

            // Skip particles that are not tracks or are too far from the
            // interaction vertex.
            if(pvars::semantic_type(p) != 1 || vertex_distance >= 6)
                continue;

            // Check if the current particle is longer than the longest found
            // so far. If so, update the longest and second longest lengths.
            if(pvars::length(p) > longest_length)
            {
                second_longest_length = longest_length;
                longest_length = pvars::length(p);
                second_index = index;
                index = i;
            }

            // If the current particle is not longer than the longest but
            // is longer than the second longest, update the second longest.
            else if(pvars::length(p) > second_longest_length)
            {
                second_longest_length = pvars::length(p);
                second_index = i;
            }
        }
        return index;
    }
    REGISTER_SELECTOR(second_longest_track, second_longest_track);

    /**
     * @brief Finds the index corresponding to the leading photon.
     * @details The leading photon is defined as the photon with the highest
     * kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading photon (highest KE).
     */
    template<class T>
    size_t leading_photon(const T & obj)
    {
        return leading_particle_index(obj, 0);
    }
    REGISTER_SELECTOR(leading_photon, leading_photon);

    /**
     * @brief Finds the index corresponding to the leading electron.
     * @details The leading electron is defined as the electron with the highest
     * kinetic energy. If the interaction is a true interaction, the initial
     * kinetic energy is used instead of the CSDA kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading electron (highest KE).
     */
    template<class T>
    size_t leading_electron(const T & obj)
    {
        return leading_particle_index(obj, 1);
    }
    REGISTER_SELECTOR(leading_electron, leading_electron);

    /**
     * @brief Finds the index corresponding to the leading muon.
     * @details The leading muon is defined as the muon with the highest
     * kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading muon (highest KE).
     */
    template<class T>
    size_t leading_muon(const T & obj)
    {
        return leading_particle_index(obj, 2);
    }
    REGISTER_SELECTOR(leading_muon, leading_muon);

    /**
     * @brief Finds the index corresponding to the leading pion.
     * @details The leading pion is defined as the pion with the highest
     * kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading pion (highest KE).
     */
    template<class T>
    size_t leading_pion(const T & obj)
    {
        return leading_particle_index(obj, 3);
    }
    REGISTER_SELECTOR(leading_pion, leading_pion);
    
    /**
     * @brief Finds the index corresponding to the leading proton.
     * @details The leading proton is defined as the proton with the highest
     * kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading proton (highest KE).
     */
    template<class T>
    size_t leading_proton(const T & obj)
    {
        return leading_particle_index(obj, 4);
    }
    REGISTER_SELECTOR(leading_proton, leading_proton);
}
#endif // SELECTORS_H