/**
 * @file cuts_muon2024.h
 * @brief Header file for definitions of analysis cuts specific to the muon2024
 * analysis.
 * @details This file contains definitions of analysis cuts which can be used
 * to select interactions specific to the muon2024 analysis. The cuts are
 * intended to be used in conjunction with the generic cuts defined in cuts.h.
 * @author mueller@fnal.gov
*/
#ifndef CUTS_MUON2024_H
#define CUTS_MUON2024_H
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "include/utilities.h"

/**
 * @namespace cuts::muon2024
 * @brief Namespace for organizing cuts specific to the muon2024 analysis.
 * @details This namespace is intended to be used for organizing cuts which act
 * on interactions specific to the muon2024 analysis. Each cut is implemented as
 * a function which takes an interaction object as an argument and returns a
 * boolean. The function should be templated on the type of interaction object if
 * the cut is intended to be used on both true and reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the cuts
 * namespace, which is used for organizing generic cuts which act on interactions.
 */
namespace cuts::muon2024
{
    /**
     * @brief Apply a 1mu1p topological (final state) cut.
     * @details The interaction must have a topology matching 1mu1p as defined by
     * the conditions in the @ref utilities::count_primaries() function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a 1mu1p topology.
     * @note This cut is intended to be used for the muon2024 analysis.
     */
    template<class T>
    bool topological_1mu1p_cut(const T & obj)
    {
        std::vector<uint32_t> c(utilities::count_primaries(obj));
        return c[0] == 0 && c[1] == 0 && c[2] == 1 && c[3] == 0 && c[4] == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, topological_1mu1p_cut, topological_1mu1p_cut);

    /**
     * @brief Apply a 1muNp topological (final state) cut.
     * @details The interaction must have a topology matching 1muNp as defined by
     * the conditions in the @ref utilities::count_primaries() function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a 1muNp topology.
     * @note This cut is intended to be used for the muon2024 analysis.
     */
    template<class T>
    bool topological_1muNp_cut(const T & obj)
    {
        std::vector<uint32_t> c(utilities::count_primaries(obj));
        return c[0] == 0 && c[1] == 0 && c[2] == 1 && c[3] == 0 && c[4] >= 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, topological_1muNp_cut, topological_1muNp_cut);
    
    /**
     * @brief Apply a 1muX topological (final state) cut.
     * @details The interaction must have a topology matching 1muX as defined by
     * the conditions in the @ref utilities::count_primaries() function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a 1muX topology.
     * @note This cut is intended to be used for the muon2024 analysis.
     */
    template<class T>
    bool topological_1muX_cut(const T & obj)
    {
        std::vector<uint32_t> c(utilities::count_primaries(obj));
        return c[2] == 1 && (c[0] > 0 || c[1] > 0 || c[3] > 0 || c[4] > 0);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, topological_1muX_cut, topological_1muX_cut);

    /**
     * @brief Apply a fiducial volume, containment, flash time (BNB), and 1mu1p
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, containment, flash time
     * (BNB), and 1mu1p topological cut on the interaction using the logical "and"
     * of each previously defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * flash time, and 1mu1p topological cut.
     * @note This cut is intended to be used for the muon2024 analysis.
     */
    template<class T>
    bool all_1mu1p_cut(const T & obj) { return fiducial_cut<T>(obj) && containment_cut<T>(obj) && flash_cut<T>(obj) && topological_1mu1p_cut<T>(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, all_1mu1p_cut, all_1mu1p_cut);

    /**
     * @brief Apply a fiducial volume, flash time (BNB), and 1mu1p
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, flash time (BNB), and
     * 1mu1p topological cut on the interaction using the logical "and" of each
     * previously defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, flash time,
     * and 1mu1p topological cut.
     * @note This cut is intended to be used for the muon2024 analysis.
     */
    template<class T>
    bool all_1mu1p_no_containment_cut(const T & obj) { return fiducial_cut<T>(obj) && flash_cut<T>(obj) && topological_1mu1p_cut<T>(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, all_1mu1p_no_containment_cut, all_1mu1p_no_containment_cut);

    /**
     * @brief Apply a fiducial volume, containment, flash time (BNB), and 1muNp
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, containment, flash time
     * (BNB), and 1muNp topological cut on the interaction using the logical "and"
     * of each previously defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * flash time, and 1muNp topological cut.
     * @note This cut is intended to be used for the muon2024 analysis.
     */
    template<class T>
    bool all_1muNp_cut(const T & obj) { return fiducial_cut<T>(obj) && containment_cut<T>(obj) && flash_cut<T>(obj) && topological_1muNp_cut<T>(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, all_1muNp_cut, all_1muNp_cut);

    /**
     * @brief Apply a fiducial volume, flash time (BNB), and 1muNp
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, flash time (BNB), and
     * 1muNp topological cut on the interaction using the logical "and" of each
     * previously defined cut.
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, flash time,
     * and 1muNp topological cut.
     * @note This cut is intended to be used for the muon2024 analysis.
     */
    template<class T>
    bool all_1muNp_no_containment_cut(const T & obj) { return fiducial_cut<T>(obj) && flash_cut<T>(obj) && topological_1muNp_cut<T>(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, all_1muNp_no_containment_cut, all_1muNp_no_containment_cut);

    /**
     * @brief Apply a fiducial volume, containment, flash time (BNB), and 1muX
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, containment, flash time
     * (BNB), and 1muX topological cut on the interaction using the logical "and"
     * of each previously defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * flash time, and 1muX topological cut.
     * @note This cut is intended to be used for the muon2024 analysis.
     */
    template<class T>
    bool all_1muX_cut(const T & obj) { return fiducial_cut<T>(obj) && containment_cut<T>(obj) && flash_cut<T>(obj) && topological_1muX_cut<T>(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, all_1muX_cut, all_1muX_cut);

    /**
     * @brief Apply a fiducial volume, flash time (BNB), and 1muX
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, flash time (BNB), and
     * 1muX topological cut on the interaction using the logical "and" of each
     * previously defined cut.
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, flash time,
     * and 1muX topological cut.
     * @note This cut is intended to be used for the muon2024 analysis.
     */
    template<class T>
    bool all_1muX_no_containment_cut(const T & obj) { return fiducial_cut<T>(obj) && flash_cut<T>(obj) && topological_1muX_cut<T>(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, all_1muX_no_containment_cut, all_1muX_no_containment_cut);

    /**
     * @brief Apply a cut to select the 1mu1p signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1mu1p signal.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1mu1p topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining the signal.
     */
    template<class T>
    bool signal_1mu1p(const T & obj) { return neutrino(obj) && fiducial_cut(obj) && containment_cut(obj) && topological_1mu1p_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, signal_1mu1p, signal_1mu1p);

    /**
     * @brief Apply a cut to select the 1mu1p signal.
     * @details This function applies a cut on the final state, and fiducial
     * volume, of the interaction. This is the "true" 1mu1p signal.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, and 1mu1p
     * topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining the signal.
     */
    template<class T>
    bool signal_1mu1p_no_containment(const T & obj) { return neutrino(obj) && fiducial_cut(obj) && topological_1mu1p_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, signal_1mu1p_no_containment, signal_1mu1p_no_containment);

    /**
     * @brief Apply a cut to select the 1mu1p non-signal.
     * @details This function applies a cut on the final state, and fiducial volume,
     * and containment of the interaction. This is the "true" 1mu1p non-signal 
     * (1mu1p topology, but not signal).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1mu1p topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining a complement to the signal.
     */
    template<class T>
    bool nonsignal_1mu1p(const T & obj) { return neutrino(obj) && !(fiducial_cut(obj) && containment_cut(obj)) && topological_1mu1p_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, nonsignal_1mu1p, nonsignal_1mu1p);

    /**
     * @brief Apply a cut to select the 1mu1p non-signal.
     * @details This function applies a cut on the final state, and fiducial
     * volume, of the interaction. This is the "true" 1mu1p non-signal 
     * (1mu1p topology, but not signal).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, and 1mu1p
     * topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining a complement to the signal.
     */
    template<class T>
    bool nonsignal_1mu1p_no_containment(const T & obj) { return neutrino(obj) && !fiducial_cut(obj) && topological_1mu1p_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, nonsignal_1mu1p_no_containment, nonsignal_1mu1p_no_containment);

    /**
     * @brief Apply a cut to select the 1muNp signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1muNp signal.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1muNp topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining the signal.
     */
    template<class T>
    bool signal_1muNp(const T & obj) { return neutrino(obj) && fiducial_cut(obj) && containment_cut(obj) && topological_1muNp_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, signal_1muNp, signal_1muNp);

    /**
     * @brief Apply a cut to select the 1muNp signal.
     * @details This function applies a cut on the final state, and fiducial
     * volume, of the interaction. This is the "true" 1muNp signal.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, and 1muNp
     * topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining the signal.
     */
    template<class T>
    bool signal_1muNp_no_containment(const T & obj) { return neutrino(obj) && fiducial_cut(obj) && topological_1muNp_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, signal_1muNp_no_containment, signal_1muNp_no_containment);

    /**
     * @brief Apply a cut to select the 1muNp non-signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1muNp non-signal
     * (1muNp topology, but not signal).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1muNp topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining a complement to the signal.
     */
    template<class T>
    bool nonsignal_1muNp(const T & obj) { return neutrino(obj) && !(fiducial_cut(obj) && containment_cut(obj)) && topological_1muNp_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, nonsignal_1muNp, nonsignal_1muNp);

    /**
     * @brief Apply a cut to select the 1muNp non-signal.
     * @details This function applies a cut on the final state, and fiducial
     * volume, of the interaction. This is the "true" 1muNp non-signal 
     * (1muNp topology, but not signal).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, and 1muNp
     * topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining a complement to the signal.
     */
    template<class T>
    bool nonsignal_1muNp_no_containment(const T & obj) { return neutrino(obj) && !fiducial_cut(obj) && topological_1muNp_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, nonsignal_1muNp_no_containment, nonsignal_1muNp_no_containment);

    /**
     * @brief Apply a cut to select the 1muX signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1muX signal
     * definition.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1muX topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining the signal.
     */
    template<class T>
    bool signal_1muX(const T & obj) { return neutrino(obj) && fiducial_cut(obj) && containment_cut(obj) && topological_1muX_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, signal_1muX, signal_1muX);

    /**
     * @brief Apply a cut to select the 1muX signal.
     * @details This function applies a cut on the final state, and fiducial
     * volume, of the interaction. This is the "true" 1muX signal definition.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, and 1muX
     * topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining the signal.
     */
    template<class T>
    bool signal_1muX_no_containment(const T & obj) { return neutrino(obj) && fiducial_cut(obj) && topological_1muX_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, signal_1muX_no_containment, signal_1muX_no_containment);

    /**
     * @brief Apply a cut to select the 1muX non-signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1muX non-signal
     * (1muX topology, but not signal).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1muX topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining a complement to the signal.
     */
    template<class T>
    bool nonsignal_1muX(const T & obj) { return neutrino(obj) && !(fiducial_cut(obj) && containment_cut(obj)) && topological_1muX_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, nonsignal_1muX, nonsignal_1muX);

    /**
     * @brief Apply a cut to select the 1muX non-signal.
     * @details This function applies a cut on the final state, and fiducial
     * volume, of the interaction. This is the "true" 1muX non-signal 
     * (1muX topology, but not signal).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, and 1muX
     * topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining a complement to the signal.
     */
    template<class T>
    bool nonsignal_1muX_no_containment(const T & obj) { return neutrino(obj) && !fiducial_cut(obj) && topological_1muX_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, nonsignal_1muX_no_containment, nonsignal_1muX_no_containment);
}
#endif // CUTS_MUON2024_H
