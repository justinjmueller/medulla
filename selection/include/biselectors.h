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
        return selectors::pi0_photon_pair(obj);
    }
    REGISTER_BISELECTOR(pi0_photon_pair, pi0_photon_pair);
}
#endif // BISELECTORS_H
