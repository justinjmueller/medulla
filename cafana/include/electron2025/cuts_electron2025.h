/**
 * @file cuts_electron2025.h
 * @brief Header file for definitions of analysis cuts specific to the electron2025 benchmarking analysis.
 * @details This file  definitions of analysis cuts which can be used
 * to select interactions specific to the muon2024 analysis. The cuts are
 * intended to be used in conjunction with the generic cuts defined in cuts.h.
 * @author jzettle@fnal.gov/mueller@fnal.gov
*/
#ifndef CUTS_ELECTRON2025_H
#define CUTS_ELECTRON2025_H
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "include/utilities.h"
#include "include/electron2025/utilities_electron2025.h"

/**
 * @namespace cuts::electron2025
 * @brief Namespace for organizing cuts specific to the electron2025 analysis.
 * @details This namespace is intended to be used for organizing cuts which act
 * on interactions specific to the electron2025 analysis. Each cut is implemented as
 * a function which takes an interaction object as an argument and returns a
 * boolean. The function should be templated on the type of interaction object if
 * the cut is intended to be used on both true and reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the cuts
 * namespace, which is used for organizing generic cuts which act on interactions.
 */
namespace cuts::electron2025
{

  /**                                                                                                                                     
     * @brief Apply a 2e topological (final state) cut.
     * @details The interaction must have a topology matching 1mu1p as defined by
     * the conditions in the @ref utilities::count_primaries() function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a 2e topology.
     * @note This cut is intended to be used for the electron2025 benchmarking.     
     */
    template<class T>
        bool topological_2e_cut(const T & obj)
        {
            std::vector<uint32_t> c(utilities::electron2025::count_primaries_ee(obj));
            return c[0] == 0 && c[1] == 2 && c[2] == 0 && c[3] == 0 && c[4] == 0;
        }

    template<class T>
        bool topological_1e1gamma_cut(const T & obj)
        {
            std::vector<uint32_t> c(utilities::electron2025::count_primaries_ee(obj));
            return c[0] == 1 && c[1] == 1 && c[2] == 0 && c[3] == 0 && c[4] == 0;
        }

    template<class T>
        bool topological_1eNgamma_cut(const T & obj)
        {
            std::vector<uint32_t> c(utilities::electron2025::count_primaries_ee(obj));
            return c[0] > 1 && c[1] == 1 && c[2] == 0 && c[3] == 0 && c[4] == 0;
        }

    template<class T>
        bool topological_2gamma_cut(const T & obj)
        {
            std::vector<uint32_t> c(utilities::electron2025::count_primaries_ee(obj));
            return c[0] == 2 && c[1] == 0 && c[2] == 0 && c[3] == 0 && c[4] == 0;
        }
    
    template<class T>
        bool topological_gt2e_cut(const T & obj)
        {
            std::vector<uint32_t> c(utilities::electron2025::count_primaries_ee(obj));
            return c[0] == 0 && c[1] > 2 && c[2] == 0 && c[3] == 0 && c[4] == 0;
        }

    template<class T>
        bool topological_gt2gamma_cut(const T & obj)
        {
            std::vector<uint32_t> c(utilities::electron2025::count_primaries_ee(obj));
            return c[0] > 2 && c[1] == 0 && c[2] == 0 && c[3] == 0 && c[4] == 0;
        }

    template<class T>
        bool topological_1shower_cut(const T & obj)
        {
            std::vector<uint32_t> c(utilities::electron2025::count_primaries_ee(obj));
            return c[0] >= 1 || c[1] >= 1;
        }
    
    template<class T>
        bool topological_1showeronly_cut(const T & obj)
        {
            std::vector<uint32_t> c(utilities::electron2025::count_primaries_ee(obj));
            return ((c[0] == 1 && c[1] == 0) || (c[0] == 0 && c[1] == 1)) && (c[2] == 0 && c[3] == 0 && c[4] == 0);
        }
    template<class T>
        bool has_no_tracks(const T & obj)
        {
            std::vector<uint32_t> c(utilities::electron2025::count_primaries_ee(obj));
            return (c[2] == 0 && c[3] == 0 && c[4] == 0);
        }
    template<class T>
        bool has_two_showers(const T & obj)
        {
            std::vector<uint32_t> c(utilities::electron2025::count_primaries_ee(obj));
            return (c[0] + c[1] == 2);
        }
    template<class T>
        bool has_two_electrons(const T & obj)
        {
            std::vector<uint32_t> c(utilities::electron2025::count_primaries_ee(obj));
            return (c[1] == 2);
        }

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
            return c[2] == 1;
        }

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


    /**
     * @brief Apply a fiducial volume, containment, flash time (BNB), and 2e
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, containment, flash time
     * (BNB), and 2e topological cut on the interaction using the logical "and"
     * of each previously defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * flash time, and 2e topological cut.
     * @note This cut is intended to be used for the electron2025 benchmarking.
     */
    template<class T>
        bool all_2e_cut_bnb(const T & obj) { return fiducial_cut<T>(obj) && containment_cut<T>(obj) && flash_cut<T>(obj) && topological_2e_cut<T>(obj); }

    /**
     * @brief Apply a fiducial volume, containment, and 2e
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, containment, and 2e topological cut on the interaction using the logical "and"
     * of each previously defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * flash time, and 2e topological cut.
     * @note This cut is intended to be used for the electron2025 benchmarking.
     */

    template<class T>
	bool all_2e_cut(const T & obj) { return fiducial_cut<T>(obj) && topological_2e_cut<T>(obj); }

    template<class T>
	bool all_1shower_cut(const T & obj) { return fiducial_cut<T>(obj) && topological_1shower_cut<T>(obj); }

    template<class T>
	bool all_1showeronly_cut(const T & obj) { return fiducial_cut<T>(obj) && topological_1showeronly_cut<T>(obj); }

    template<class T>
	bool all_1e1gamma_cut(const T & obj) { return fiducial_cut<T>(obj) && topological_1e1gamma_cut<T>(obj); }

    template<class T>
	bool all_1eNgamma_cut(const T & obj) { return fiducial_cut<T>(obj) && topological_1eNgamma_cut<T>(obj); }
  
    template<class T>
	bool all_2gamma_cut(const T & obj) { return fiducial_cut<T>(obj) && topological_2gamma_cut<T>(obj); }

    template<class T>
	bool all_gt2e_cut(const T & obj) { return fiducial_cut<T>(obj) && topological_gt2e_cut<T>(obj); }

    template<class T>
	bool all_gt2gamma_cut(const T & obj) { return fiducial_cut<T>(obj) && topological_gt2gamma_cut<T>(obj); }

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

    /**
     * @brief Apply a cut to select the 1mu1p signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1mu1p signal.
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1mu1p topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining the signal.
     */
    bool signal_1mu1p(const caf::SRInteractionTruthDLPProxy & obj) { return neutrino(obj) && fiducial_cut(obj) && containment_cut(obj) && topological_1mu1p_cut(obj); }

    bool signal_1shower(const caf::SRInteractionTruthDLPProxy & obj) { return fiducial_cut(obj) && containment_cut(obj) && topological_1shower_cut(obj); }

    bool signal_2electrons(const caf::SRInteractionTruthDLPProxy & obj) { return fiducial_cut(obj) && containment_cut(obj) && topological_2e_cut(obj); }

    /**
     * @brief Apply a cut to select the 1mu1p non-signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1mu1p non-signal 
     * (1mu1p topology, but not signal).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1mu1p topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining a complement to the signal.
     */
    bool nonsignal_1mu1p(const caf::SRInteractionTruthDLPProxy & obj) { return neutrino(obj) && !(fiducial_cut(obj) && containment_cut(obj)) && topological_1mu1p_cut(obj); }

    /**
     * @brief Apply a cut to select the 1muNp signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1muNp signal.
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1muNp topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining the signal.
     */
    bool signal_1muNp(const caf::SRInteractionTruthDLPProxy & obj) { return neutrino(obj) && fiducial_cut(obj) && containment_cut(obj) && topological_1muNp_cut(obj); }

    /**
     * @brief Apply a cut to select the 1muNp non-signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1muNp non-signal
     * (1muNp topology, but not signal).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1muNp topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining a complement to the signal.
     */
    bool nonsignal_1muNp(const caf::SRInteractionTruthDLPProxy & obj) { return neutrino(obj) && !(fiducial_cut(obj) && containment_cut(obj)) && topological_1muNp_cut(obj); }

    /**
     * @brief Apply a cut to select the 1muX signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1muX signal
     * definition.
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1muX topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining the signal.
     */
    bool signal_1muX(const caf::SRInteractionTruthDLPProxy & obj) { return neutrino(obj) && fiducial_cut(obj) && containment_cut(obj) && topological_1muX_cut(obj); }

    /**
     * @brief Apply a cut to select the 1muX non-signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1muX non-signal
     * (1muX topology, but not signal).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1muX topological cut.
     * @note This cut is intended to be used for the muon2024 analysis for
     * defining a complement to the signal.
     */
    bool nonsignal_1muX(const caf::SRInteractionTruthDLPProxy & obj) { return neutrino(obj) && !(fiducial_cut(obj) && containment_cut(obj)) && topological_1muX_cut(obj); }
}
#endif // CUTS_ELECTRON2025_H
