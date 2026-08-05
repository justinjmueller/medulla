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
     * @brief Finds the index corresponding to the leading primary particle of
     * the specified particle type.
     * @details Like leading_particle_index, but restricted to particles
     * satisfying pvars::primary_classification.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @param pid the particle type.
     * @return the index of the leading primary particle (highest KE).
     */
    template <class T>
    size_t leading_primary_particle_index(const T & obj, uint16_t pid)
    {
        double leading_ke(0);
        size_t index(kNoMatch);
        for(size_t i(0); i < obj.particles.size(); ++i)
        {
            const auto & p = obj.particles[i];
            double energy(pvars::ke(p));
            if(pvars::pid(p) == pid && pvars::primary_classification(p) && energy > leading_ke)
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
     * required to be a primary particle with a semantic type of 1 (track).
     * No requirement is made on the particle's proximity to the interaction
     * vertex.
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

            // Skip particles that are not primary tracks.
            if(pvars::semantic_type(p) != 1 || !pvars::primary_classification(p))
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
     * particle instance is required to be a primary particle with a semantic
     * type of 1 (track). No requirement is made on the particle's proximity
     * to the interaction vertex.
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

            // Skip particles that are not primary tracks.
            if(pvars::semantic_type(p) != 1 || !pvars::primary_classification(p))
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
        return second_index;
    }
    REGISTER_SELECTOR(second_longest_track, second_longest_track);

    /**
     * @brief Finds the index corresponding to the third longest track.
     * @details The third longest track is defined as the track with the
     * third longest length, which is calculated upstream in SPINE. The
     * particle instance is required to be a primary particle with a semantic
     * type of 1 (track). No requirement is made on the particle's proximity
     * to the interaction vertex.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the third longest track.
    */
    template<class T>
    size_t third_longest_track(const T & obj)
    {
        double longest_length(0);
        double second_longest_length(0);
        double third_longest_length(0);
        size_t index(kNoMatch), second_index(kNoMatch), third_index(kNoMatch);
        for(size_t i(0); i < obj.particles.size(); ++i)
        {
            const auto & p = obj.particles[i];

            // Skip particles that are not primary tracks.
            if(pvars::semantic_type(p) != 1 || !pvars::primary_classification(p))
                continue;

            // Check if the current particle is longer than the longest found
            // so far. If so, shift the longest and second longest down into
            // the second and third longest slots.
            if(pvars::length(p) > longest_length)
            {
                third_longest_length = second_longest_length;
                third_index = second_index;
                second_longest_length = longest_length;
                second_index = index;
                longest_length = pvars::length(p);
                index = i;
            }

            // If the current particle is not longer than the longest but is
            // longer than the second longest, shift the second longest down
            // into the third longest slot.
            else if(pvars::length(p) > second_longest_length)
            {
                third_longest_length = second_longest_length;
                third_index = second_index;
                second_longest_length = pvars::length(p);
                second_index = i;
            }

            // If the current particle is not longer than the second longest
            // but is longer than the third longest, update the third longest.
            else if(pvars::length(p) > third_longest_length)
            {
                third_longest_length = pvars::length(p);
                third_index = i;
            }
        }
        return third_index;
    }
    REGISTER_SELECTOR(third_longest_track, third_longest_track);

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
        return leading_particle_index(obj, pvars::kPhoton);
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
        return leading_particle_index(obj, pvars::kElectron);
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
        return leading_particle_index(obj, pvars::kMuon);
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
        return leading_particle_index(obj, pvars::kPion);
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
        return leading_particle_index(obj, pvars::kProton);
    }
    REGISTER_SELECTOR(leading_proton, leading_proton);

    /**
     * @brief Finds the index corresponding to the leading primary photon.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading primary photon (highest KE).
     */
    template<class T>
    size_t leading_primary_photon(const T & obj)
    {
        return leading_primary_particle_index(obj, pvars::kPhoton);
    }
    REGISTER_SELECTOR(leading_primary_photon, leading_primary_photon);

    /**
     * @brief Finds the index corresponding to the leading primary electron.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading primary electron (highest KE).
     */
    template<class T>
    size_t leading_primary_electron(const T & obj)
    {
        return leading_primary_particle_index(obj, pvars::kElectron);
    }
    REGISTER_SELECTOR(leading_primary_electron, leading_primary_electron);

    /**
     * @brief Finds the index corresponding to the leading primary muon.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading primary muon (highest KE).
     */
    template<class T>
    size_t leading_primary_muon(const T & obj)
    {
        return leading_primary_particle_index(obj, pvars::kMuon);
    }
    REGISTER_SELECTOR(leading_primary_muon, leading_primary_muon);

    /**
     * @brief Finds the index corresponding to the leading primary pion.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading primary pion (highest KE).
     */
    template<class T>
    size_t leading_primary_pion(const T & obj)
    {
        return leading_primary_particle_index(obj, pvars::kPion);
    }
    REGISTER_SELECTOR(leading_primary_pion, leading_primary_pion);

    /**
     * @brief Finds the index corresponding to the leading primary proton.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the leading primary proton (highest KE).
     */
    template<class T>
    size_t leading_primary_proton(const T & obj)
    {
        return leading_primary_particle_index(obj, pvars::kProton);
    }
    REGISTER_SELECTOR(leading_primary_proton, leading_primary_proton);

    /**
     * @brief Finds the index corresponding to the target Michel.
     * @details The target Michel is defined as the Michel with the most
     * depositions in the interaction. 
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return the index of the target Michel (largest).
     */
    template<class T>
    size_t target_michel(const T & obj)
    {
        double largest_size(0);
        size_t index(kNoMatch);
        for(size_t i(0); i < obj.particles.size(); ++i)
        {
            const auto & p = obj.particles[i];
            double size(p.size);
            if(pvars::semantic_type(p) == 2 && size > largest_size)
            {
                largest_size = size;
                index = i;
            }
        }
        return index;
    }
    REGISTER_SELECTOR(target_michel, target_michel);

    /**
     * @brief Selects the leading and subleading photon forming the best pi0 candidate.
     * Helper function. Not registered.
     * @details
     * Reco branch: iterates all ordered primary-photon pairs above a 25 MeV
     * per-shower threshold, computes the diphoton invariant mass using the
     * vertex-to-shower-start opening angle, and selects the pair whose mass
     * is closest to PI0_MASS (135 MeV). Leading photon has higher calo KE.
     *
     * True branch: groups photon daughters by parent pi0 track ID via
     * utilities::get_true_pi0s, requires exactly two photon daughters, and
     * orders them by true KE.
     *
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return pair of indices {leading_photon, subleading_photon}, or
     *         {kNoMatch, kNoMatch} if no valid pair is found.
     */
    template<class T>
    std::pair<size_t, size_t> pi0_photon_pair(const T & obj)
    {
        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
        {
            auto true_primary_pi0s = utilities::get_true_pi0s(obj, true);

            int num_photon_daughters = 0;
            std::vector<size_t> daughter_indices;
            for(const auto & entry : true_primary_pi0s)
                for(size_t idx : entry.second)
                {
                    daughter_indices.push_back(idx);
                    if(obj.particles[idx].pid == 0) ++num_photon_daughters;
                }

            if(num_photon_daughters != 2)
                return {kNoMatch, kNoMatch};

            const auto & d0 = obj.particles[daughter_indices[0]];
            const auto & d1 = obj.particles[daughter_indices[1]];
            if(d0.ke > d1.ke)
                return {daughter_indices[0], daughter_indices[1]};
            else
                return {daughter_indices[1], daughter_indices[0]};
        }
        else
        {
            constexpr double threshold = 25.0;
            double vx = obj.vertex[0], vy = obj.vertex[1], vz = obj.vertex[2];

            std::vector<std::pair<std::pair<size_t,size_t>, double>> candidates;
            for(size_t i = 0; i < obj.particles.size(); ++i)
            {
                const auto & p = obj.particles[i];
                if(!(p.is_primary && p.pid == 0)) continue;

                double dx0 = p.start_point[0] - vx;
                double dy0 = p.start_point[1] - vy;
                double dz0 = p.start_point[2] - vz;
                double r0  = std::sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);

                for(size_t j = 0; j < obj.particles.size(); ++j)
                {
                    if(j == i) continue;
                    const auto & q = obj.particles[j];
                    if(!(q.is_primary && q.pid == 0)) continue;

                    double ke_p = pvars::calo_ke(p), ke_q = pvars::calo_ke(q);
                    double leading_ke    = (ke_p > ke_q) ? ke_p : ke_q;
                    double subleading_ke = (ke_p > ke_q) ? ke_q : ke_p;
                    if(leading_ke < threshold || subleading_ke < threshold) continue;

                    double dx1 = q.start_point[0] - vx;
                    double dy1 = q.start_point[1] - vy;
                    double dz1 = q.start_point[2] - vz;
                    double r1  = std::sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);

                    double costheta = (dx0*dx1 + dy0*dy1 + dz0*dz1) / (r0 * r1);
                    double mass = std::sqrt(2.0 * leading_ke * subleading_ke * (1.0 - costheta));
                    
                    // TODO: this currently results in a "no-op" when sorting
                    // candidate pairs. This is needed to match Lane's thesis
                    // analysis, but should be revisited in the future to see
                    // if it can be improved by using a more sophisticated
                    // metric for selecting the best pi0 candidate.
                    candidates.push_back({{i, j}, PLACEHOLDERVALUE});
                }
            }

            if(candidates.empty()) return {kNoMatch, kNoMatch};

            std::sort(candidates.begin(), candidates.end(),
                [](const auto & a, const auto & b) {
                    return std::abs(a.second - PI0_MASS) < std::abs(b.second - PI0_MASS);
                });

            auto [idx0, idx1] = candidates[0].first;
            double ke0 = pvars::calo_ke(obj.particles[idx0]);
            double ke1 = pvars::calo_ke(obj.particles[idx1]);
            return (ke0 > ke1) ? std::make_pair(idx0, idx1)
                               : std::make_pair(idx1, idx0);
        }
    }

    /**
     * @brief Find the index corresponding to the pi0 leading shower.
     * @details The leading shower is described as the pi0 daughter
     * shower with the highest kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     */
    template<class T>
    size_t pi0_leading_shower(const T & obj)
    {
        return pi0_photon_pair(obj).first;
    }
    REGISTER_SELECTOR(pi0_leading_shower, pi0_leading_shower);

    /**
     * @brief Find the index corresponding to the pi0 subleading shower.
     * @details The leading shower is described as the pi0 daughter
     * shower with the lowest kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     */
    template<class T>
    size_t pi0_subleading_shower(const T & obj)
    {
        return pi0_photon_pair(obj).second;
    }
    REGISTER_SELECTOR(pi0_subleading_shower, pi0_subleading_shower);
}
#endif // SELECTORS_H
