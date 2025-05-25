/**
 * @file cuts.h
 * @brief Header file for definitions of analysis cuts.
 * @details This file contains definitions of analysis cuts which can be used
 * to select interactions. Each cut is implemented as a function which takes an
 * interaction object as an argument and returns a boolean. These are the
 * building blocks for defining more complex selections.
 * @author mueller@fnal.gov
*/
#ifndef CUTS_H
#define CUTS_H
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "utilities.h"
#include "framework.h"

/**
 * @namespace cuts
 * @brief Namespace for organizing generic cuts which act on interactions.
 * @details This namespace is intended to be used for organizing cuts which act
 * on interactions. Each cut is implemented as a function which takes an
 * interaction object as an argument and returns a boolean. The function should
 * be templated on the type of interaction object if the cut is intended to be
 * used on both true and reconstructed interactions.
 */
namespace cuts
{   
    /**
     * @brief Apply a cut on the validity of the flash match.
     * @details A "valid" flash match is defined as a flash-interaction
     * association with a flash time that is not NaN and a flash match
     * status of 1. The upstream flash matching algorithm (OpT0Finder) has a
     * flash filter that restricts candidate flashes to near the beam window,
     * which means that the majority of cosmogenic interactions are not
     * flash matched. If no flash match is found, the flash time is NaN. This
     * cut is intended to be applied as a preselection cut to reduce comparisons
     * to NaN values, which tend to be noisy on stderr.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction on which to place the flash validity cut.
     * @return true if the interaction is flash matched and the time is valid.
     */
    template<class T>
    bool valid_flashmatch(const T & obj)
    {
        return obj.flash_times.size() > 0 && obj.is_flash_matched == 1 && !std::isnan(obj.flash_times[0]);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, valid_flashmatch, valid_flashmatch);

    /**
     * @brief Apply no cut; all interactions passed.
     * @details This is a placeholder function for a cut which does not apply
     * any selection criteria. It is intended to be used in cases where a cut
     * function is required, but no selection is desired.
     * @tparam T the type of object (true or reco).
     * @param obj the interaction to select on.
     * @return true (always).
     */
    template<class T>
    bool no_cut(const T & obj) { return true; }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, no_cut, no_cut);

    /**
     * @brief Apply a cut to select neutrinos.
     * @details This function applies a cut to select neutrinos. This cut
     * makes use of the is_neutrino flag in the true interaction object and is
     * intended to be used to identify signal neutrinos.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction is a neutrino.
     * @note This cut is intended to be used for identifying neutrinos in
     * truth, which is useful for making signal definitions.
     */
    template<class T>
    bool neutrino(const T & obj) { return obj.nu_id >= 0; }
    REGISTER_CUT_SCOPE(RegistrationScope::True, neutrino, neutrino);

    /**
     * @brief Apply a cut to select cosmogenic interactions.
     * @details This function applies a cut to select cosmogenic interactions.
     * This cut makes use of the is_neutrino flag in the true interaction
     * object and is intended to be used to identify cosmogenic interactions.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction is a cosmogenic interaction.
     * @note This cut is intended to be used for identifying cosmogenic
     * interactions in truth, which is useful for making background definitions.
     */
    template<class T>
    bool cosmic(const T & obj) { return !neutrino(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, cosmic, cosmic);

    /**
     * @brief Apply a cut to select charged current interactions.
     * @details This function applies a cut to select charged current
     * interactions. This cut makes use of the `current_type` attribute in the
     * true interaction object and is intended to be used to identify charged
     * current interactions.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction is a charged current interaction.
     */
    template<class T>
    bool iscc(const T & obj) { return obj.current_type == 0; }
    REGISTER_CUT_SCOPE(RegistrationScope::True, iscc, iscc);

    /**
     * @brief Apply a fiducial volume cut; the interaction vertex must be
     * reconstructed within the fiducial volume.
     * @details The fiducial volume cut is applied on the reconstructed
     * interaction vertex upstream in SPINE. The fiducial volume is defined
     * (in a SPINE post-processor) as a 25 cm border around the x and y
     * detector faces, a 50 cm border around the downstream (+) z face, and a
     * 30 cm border around the upstream (-) z face. The fiducial volume is
     * intended to reduce the impact of detector edge effects on the analysis.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the vertex is in the fiducial volume.
     */
    template<class T>
    bool fiducial_cut(const T & obj)
    {
        return obj.is_fiducial && !(obj.vertex[0] > 210.215 && obj.vertex[1] > 60 && (obj.vertex[2] > 290 && obj.vertex[2] < 390));
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, fiducial_cut, fiducial_cut);
    
    /**
     * @brief Apply a containment cut on the entire interaction.
     * @details The containment cut is applied on the entire interaction. The
     * interaction is considered contained if all particles and all spacepoints
     * are contained within 5cm of the detector edges (configured in a SPINE 
     * post-processor). Additionally, no spacepoints are allowed to be
     * reconstructed in a TPC that did not create it. This is an unphysical
     * condition that can occur when a cosmic muon is moved according to an
     * assumed t0 that is very out-of-time.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the vertex is contained.
     */
    template<class T>
    bool containment_cut(const T & obj) { return obj.is_contained; }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, containment_cut, containment_cut);

    /**
     * @brief Apply a flash time cut on the interaction.
     * @details The flash time cut is applied on the interaction. The flash time
     * is required to be within the beam window, which is expected to be
     * [0 us, 1.6 us] for BNB and [0 us, 9.6 us] for NuMI. This cut is intended
     * to reduce the impact of cosmogenic interactions on analyses.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params the parameters for the cut.
     * @return true if the interaction has been matched to an in-time flash.
     * @note The switch to the NuMI beam window is applied by the definition of
     * a preprocessor macro (BEAM_IS_NUMI).
     * @note The cut window has been widened to reconcile the beam window as
     * observed in data and simulation.
     */
    template<class T>
    bool flash_cut(const T & obj, std::vector<double> params={})
    {
        if(!valid_flashmatch(obj))
            return false;
        else if(params.size() == 2 && obj.flash_times[0] >= params[0] && obj.flash_times[0] <= params[1])
            return true;
        else if(params.size() !=2)
            return true;
        else
            return false;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, flash_cut, flash_cut);

    /**
     * @brief Apply a cut to select interactions with no primary charged pions.
     * @details This function applies a cut to select interactions with no
     * primary charged pions in the final state as defined by the
     * @ref utilities::count_primaries function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has zero charged pions.
     */
    template<class T>
    bool no_charged_pions(const T & obj)
    {
        std::vector<uint32_t> c(utilities::count_primaries(obj));
        return c[3] == 0;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, no_charged_pions, no_charged_pions);

    /**
     * @brief Apply a cut to select interactions with no primary showers.
     * @details This function applies a cut to select interactions with no
     * primary showers in the final state as defined by the
     * @ref utilities::count_primaries function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has zero showers.
     */
    template<class T>
    bool no_showers(const T & obj)
    {
        std::vector<uint32_t> c(utilities::count_primaries(obj));
        return c[0] == 0 && c[1] == 0;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, no_showers, no_showers);

    /**
     * @brief Apply a cut to select interactions with a single primary muon.
     * @details This function applies a cut to select interactions with a
     * single primary muon in the final state as defined by the
     * @ref utilities::count_primaries function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a single primary muon.
     */
    template<class T>
    bool has_single_muon(const T & obj)
    {
        std::vector<uint32_t> c(utilities::count_primaries(obj));
        return c[2] == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, has_single_muon, has_single_muon);

    /**
     * @brief Apply a cut to select interactions with a single primary proton.
     * @details This function applies a cut to select interactions with a
     * single primary proton in the final state as defined by the
     * @ref utilities::count_primaries function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a single primary proton.
     */
    template<class T>
    bool has_single_proton(const T & obj)
    {
        std::vector<uint32_t> c(utilities::count_primaries(obj));
        return c[4] == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, has_single_proton, has_single_proton);

    /**
     * @brief Apply a cut to select interactions with nonzero primary protons.
     * @details This function applies a cut to select interactions with
     * nonzero primary protons in the final state as defined by the
     * @ref utilities::count_primaries function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has nonzero primary protons.
     */
    template<class T>
    bool has_nonzero_protons(const T & obj)
    {
        std::vector<uint32_t> c(utilities::count_primaries(obj));
        return c[4] > 0;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, has_nonzero_protons, has_nonzero_protons);

    /**
     * @brief Apply a cut to select interactions with at least one primary
     * photon.
     * @details This function applies a cut to select interactions with at
     * least one primary photon in the final state as defined by the
     * @ref utilities::count_primaries function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has at least one primary photon.
     */
    template<class T>
    bool has_photon(const T & obj)
    {
        std::vector<uint32_t> c(utilities::count_primaries(obj));
        return c[0] > 0;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, has_photon, has_photon);

    /**
     * @brief Apply a cut to select interactions with at least one primary
     * electron.
     * @details This function applies a cut to select interactions with at
     * least one primary electron in the final state as defined by the
     * @ref utilities::count_primaries function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has at least one primary electron.
     */
    template<class T>
    bool has_electron(const T & obj)
    {
        std::vector<uint32_t> c(utilities::count_primaries(obj));
        return c[1] > 0;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, has_electron, has_electron);
}
#endif
