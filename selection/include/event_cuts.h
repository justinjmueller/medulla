/**
 * @file event_cuts.h
 * @brief Definitions of analysis cuts which can be applied to the
 * StandardRecord object.
 * @details This file contains definitions of analysis cuts which can be
 * applied to the StandardRecord object. Each cut is implemented as a
 * function which takes a StandardRecord object as an argument and returns
 * a boolean indicating whether the event passes the cut.
 * @author mueller@fnal.gov
 */
#ifndef EVENT_CUTS_H
#define EVENT_CUTS_H
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "framework.h"
#include "utilities.h"

/**
 * @namespace ecut
 * @brief Namespace for organizing cuts which act on events.
 * @details This namespace is intended to be used for organizing cuts
 * which act on events. Each cut is implemented as a function which takes
 * a StandardRecord object as an argument and returns a boolean.
 */
namespace ecut
{
    /**
     * @brief Apply no cut to the event.
     * @details This cut always returns true, indicating that the event passes
     * the cut. It is useful for testing purposes or when no cut is desired.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the cut on.
     * @return true, indicating that the event passes the cut.
     */
    template<typename T>
    bool no_cut(const T & sr) { return true; }
    REGISTER_CUT_SCOPE(RegistrationScope::Event, no_cut, no_cut);

    /**
     * @brief Apply a data quality cut on the event metadata.
     * @details This cut rejects events that belong to runs that are marked as
     * bad. This can happen due to various reasons, such as hardware issues,
     * non-physics configuration, or other data quality concerns.
     * @param sr the StandardRecord to apply the cut on.
     * @return true if the event is from a good run, false otherwise.
     */
    template<typename T>
    bool data_quality_cut(const T & sr)
    {
        if(sr.hdr.det == caf::kICARUS && !sr.hdr.ismc)
        {
            return utilities::is_icarus_good_run(sr.hdr.run);
        }
        else
        {
            // For other detectors, we assume all runs are good.
            // This can be extended to include other detectors' run quality checks.
            return true;
        }
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Event, data_quality_cut, data_quality_cut);
}

#endif // EVENT_CUTS_H