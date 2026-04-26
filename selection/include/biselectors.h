/**
 * @file biselectors.h
 * @brief Header file for biselectors used in the SPINE analysis framework.
 * @details This file contains definitions of biselectors which can be used
 * to select a pair of particles (by index) within an interaction. This is a
 * useful feature for reducing down a collection of particles in a final state
 * to just a pair, which allows a user to broadcast a two-particle variable
 * "upwards" to the interaction level.
 * @author mueller@fnal.gov
 */
#ifndef BISELECTORS_H
#define BISELECTORS_H
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>

#include "framework.h"
#include "include/selectors.h"
#include "include/utilities.h"

/**
 * @namespace biselectors
 * @brief Namespace for organizing biselectors which act on interactions.
 * @details This namespace is intended to be used for organizing biselectors
 * which act on interactions. Each biselector is implemented as a function
 * which takes an interaction object as an argument and returns a pair of
 * indices corresponding to the two selected particles. The function should
 * be templated on the type of interaction object if the biselector is
 * intended to be used on both true and reconstructed interactions.
 */
namespace biselectors
{
    /**
     * @brief Selects the leading muon and leading proton.
     * @details The leading muon and proton are defined as the particles with
     * the highest kinetic energy of their respective types.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return pair of indices: {leading_muon, leading_proton}.
     */
    template<class T>
    std::pair<size_t, size_t> muon_proton(const T & obj)
    {
        return { selectors::leading_muon(obj), selectors::leading_proton(obj) };
    }
    REGISTER_BISELECTOR(muon_proton, muon_proton);

    /**
     * @brief Selects the two longest tracks in the interaction.
     * @details The longest and second longest tracks are defined by their
     * track length as calculated upstream in SPINE.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return pair of indices: {longest_track, second_longest_track}.
     */
    template<class T>
    std::pair<size_t, size_t> two_longest_tracks(const T & obj)
    {
        return { selectors::longest_track(obj), selectors::second_longest_track(obj) };
    }
    REGISTER_BISELECTOR(two_longest_tracks, two_longest_tracks);

    /**
     * @brief Selects the leading muon and leading pion.
     * @details The leading muon and pion are defined as the particles with
     * the highest kinetic energy of their respective types.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to operate on.
     * @return pair of indices: {leading_muon, leading_pion}.
     */
    template<class T>
    std::pair<size_t, size_t> muon_pion(const T & obj)
    {
        return { selectors::leading_muon(obj), selectors::leading_pion(obj) };
    }
    REGISTER_BISELECTOR(muon_pion, muon_pion);

    /**
     * @brief Selects the leading and subleading photon forming the best pi0 candidate.
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
    REGISTER_BISELECTOR(pi0_photon_pair, pi0_photon_pair);
}
#endif // BISELECTORS_H
