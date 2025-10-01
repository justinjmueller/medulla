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
#include "event_variables.h"

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
     * @brief Apply a cut enforcing that there are nonzero reco interactions.
     * @details This cut checks if the number of reco interactions in the event
     * is greater than zero. It is useful for ensuring that the file isn't
     * corrupted or empty.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the cut on.
     * @return true if there are nonzero reco interactions, false otherwise.
     */
    template<typename T>
    bool nonzero_reco_interactions(const T & sr)
    {
        // Check if the number of reco interactions is greater than zero.
        return (sr.ndlp > 0);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Event, nonzero_reco_interactions, nonzero_reco_interactions);

    /**
     * @brief Apply a cut enforcing a CRT-PMT veto.
     * @details This cut rejects events that do not have at least one flash in-
     * time with the beam that is not associated with a CRT hit.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the cut on.
     * @param params a vector of parameters for the cut, which can be used to
     * specify the time window for the flash association.
     * @return true if the event has at least one flash in-time with the beam
     * that is not associated with a CRT hit, false otherwise.
     */
    template<typename T>
    bool crtpmt_veto(const T & sr, std::vector<double> params={})
    {
        if(params.size() < 2)
        {
            throw std::invalid_argument("crtpmt_veto requires at least two parameters: time window start and end.");
        }
        bool crtpmt_matched(sr.ncrtpmt_matches == 0);
        for(auto const & c : sr.crtpmt_matches)
        {
            if(c.flashGateTime > params[0] && c.flashGateTime < params[1] && c.flashClassification == 0)
                crtpmt_matched = true;
        }
        return crtpmt_matched;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Event, crtpmt_veto, crtpmt_veto);

    /**
     * @brief Apply a cut on the global trigger time.
     * @details This cut checks if the global trigger time is within some
     * interval of time. This is useful for partitioning the dataset into
     * subsets based on the time of the global trigger. If no time interval is
     * specified, the cut will always return true, allowing all events to pass.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the cut on.
     * @return true if the global trigger time is within the specified interval,
     */
    template<typename T>
    bool global_trigger_time_cut(const T & sr, std::vector<double> params={})
    {
        if(params.empty())
        {
            // If no parameters are provided, we assume the cut passes for all events.
            return true;
        }
        else if(params.size() == 2)
        {
            // Check if the global trigger time is within the specified interval.
            double global_trigger_time = sr.hdr.triggerinfo.global_trigger_time;
            return (global_trigger_time >= params[0] && global_trigger_time <= params[1]);
        }
        else
        {
            throw std::invalid_argument("global_trigger_time_cut requires either no parameters or exactly two parameters.");
        }
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Event, global_trigger_time_cut, global_trigger_time_cut);

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

    /**
     * @brief Cut that implements trigger emulation on simulated data.
     * @details This cut emulates the trigger logic on simulated data by
     * checking if the event has a valid trigger time.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the cut on.
     * @return true if the event has a valid trigger time, false otherwise.
     */
    template<typename T>
    bool trigger_emulation_cut(const T & sr)
    {
        if(sr.hdr.ismc && std::isnan(sr.hdr.triggerinfo.global_trigger_det_time))
        {
            // If the event is simulated and the global trigger time is NaN,
            // we assume the trigger emulation cut fails.
            return false;
        }
        else
        {
            // For real data or if the trigger time is valid, we pass the cut.
            return true;
        }
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Event, trigger_emulation_cut, trigger_emulation_cut);

    /**
     * @brief A cut that places a threshold on the BNB Figure of Merit 2 (FoM2).
     * @details This cut checks if the FoM2 value is above a certain threshold,
     * which can be configured by the user. This is intended to be a handle on
     * the quality of the beam-target overlap of the spill that is associated
     * with the event
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the cut on.
     * @param params a vector of parameters for the cut, which can be used to
     * specify the threshold for the FoM2 value.
     * @return true if the FoM2 value is above the threshold, false otherwise.
     */
    /*template<typename T>
    bool bnb_fom2_cut(const T & sr, std::vector<double> params={})
    {
        if(params.empty())
        {
            throw std::invalid_argument("bnb_fom2_cut requires at least one parameter for the threshold (recommended 0.98).");
        }
        if(sr.hdr.ismc)
        {
            // If the event is simulated, do not apply the cut.
            return true;
        }
        else
        {
            double threshold = params[0];
            return (evar::bnb_fom2(sr) >= threshold);
        }
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Event, bnb_fom2_cut, bnb_fom2_cut);*/

    /**
     * @brief A cut that checks if the event is the first in the subrun.
     * @details This cut checks if the event is the first event in the subrun.
     * This can be relevant for certain attributes (e.g. spill info) that is
     * stored only once per subrun in the first event.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the cut on.
     * @return true if the event is the first in the subrun, false otherwise.
     */
    template<typename T>
    bool is_first_in_subrun_cut(const T & sr)
    {
        // This cut checks if the event is the first in the subrun.
        // It returns true if it is the first event in the subrun, false otherwise.
        return sr.hdr.first_in_subrun;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Event, is_first_in_subrun_cut, is_first_in_subrun_cut);
}

#endif // EVENT_CUTS_H