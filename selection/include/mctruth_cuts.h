/**
 * @file mctruth_cuts.h
 * @brief Definitions of analysis cuts which can extract information from
 * the SRTrueInteraction object.
 * @details This file contains definitions of analysis cuts which can be
 * used to extract information from the SRTrueInteraction object. Each cut
 * is implemented as a function which takes an SRTrueInteraction object as an
 * argument and returns a bool. The association of an SRInteractionTruthDLP
 * object to an SRTrueInteraction object is handled upstream in the SpineVar
 * functions.
 * @author mueller@fnal.gov
 * @author rvizarr@fnal.gov
 */
#ifndef MCTRUTH_CUTS_H
#define MCTRUTH_CUTS_H
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"

#include "framework.h"

/**
 * @namespace mctruth
 * @brief Namespace for organizing cuts which act on true interactions.
 * @details This namespace is intended to be used for organizing cuts
 * which act on true interactions. Each cut is implemented as a function
 * which takes an SRTrueInteraction object as an argument and returns a bool.
 */
namespace mctruth
{
    /**
     * @brief Cut for charged current interactions at the generator (GENIE) level.
     * @details Distinct from the SPINE truth-level iscc cut. Uses obj.iscc
     * directly from the MCTruth object.
     * @tparam T the type of the object to apply the cut on.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @return true if the interaction is charged current.
     */
    template<typename T>
        bool iscc(const T & obj) { return obj.iscc; }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, iscc, iscc);

} // namespace mctruth
#endif