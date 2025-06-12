/**
 * @file cuts_cosmics.h
 * @brief Definitions of cuts for cosmic muon studies.
 * @details This file contains definitions of cuts which can be used to select
 * cosmic muons in the ICARUS TPCs. The cuts are intended to be used in
 * conjunction with the generic cuts defined in cuts.h.
 * @author mueller@fnal.gov
 */
#ifndef CUTS_COSMICS_H
#define CUTS_COSMICS_H

#include "include/cuts.h"
#include "include/utilities.h"

/**
 * @namespace cuts::cosmics
 * @brief Namespace for organizing cuts specific to cosmic muon studies.
 * @details This namespace is intended to be used for organizing cuts which act
 * on cosmic muon interactions in the ICARUS TPCs. Each cut is implemented as a
 * function which takes an interaction object as an argument and returns a
 * boolean. 
 */
namespace cuts::cosmics
{
    /**
     * @brief Apply a cut to select cosmic muons.
     * @details This function applies a cut to select cosmic muons based on the
     * topology of the interaction. The interaction must have a topology matching
     * a cosmic muon as defined by the conditions in the @ref utilities::count_primaries()
     * function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a cosmic muon topology.
     */
    template<class T>
        bool single_cosmic_muon_cut(const T & obj)
        {
            std::vector<uint32_t> c(utilities::count_primaries(obj));
            return obj.nu_id < 0 && c[0] == 0 && c[1] == 0 && c[2] == 1 && c[3] == 0 && c[4] == 0;
        }
}
#endif // CUTS_COSMICS_H