/**
 * @file event_variables.h
 * @brief Definitions of analysis variables which can extract information from
 * the StandardRecord object.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from the StandardRecord object. Each variable
 * is implemented as a function which takes a StandardRecord object as an
 * argument and returns a double.
 * @author mueller@fnal.gov
 */
#ifndef EVENT_VARIABLES_H
#define EVENT_VARIABLES_H
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "framework.h"

/**
 * @namespace evar
 * @brief Namespace for organizing variables which act on events.
 * @details This namespace is intended to be used for organizing variables
 * which act on events. Each variable is implemented as a function which takes
 * a StandardRecord object as an argument and returns a double.
 */
namespace evar
{
    /**
     * @brief Variable for the number of true SPINE interactions in the event.
     * @details This variable counts the number of true SPINE interactions in
     * the event.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @return the number of true SPINE interactions in the event.
     */
    template<typename T>
    double ntrue(const T & sr) { return sr.ndlp_true; }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, ntrue, ntrue);

    /**
     * @brief Variable for the number of reco SPINE interactions in the event.
     * @details This variable counts the number of reco SPINE interactions in
     * the event.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @return the number of reco SPINE interactions in the event.
     */
    template<typename T>
    double nreco(const T & sr) { return sr.ndlp; }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, nreco, nreco);

    /**
     * @brief Variable for the time of the global trigger.
     * @details This variable returns the time of the global trigger in Unix
     * epoch format (nanoseconds since 1970-01-01T00:00:00Z). This is useful
     * for dividing the dataset into different "epochs" based on the absolute
     * time of the global trigger.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @return the time of the global trigger in nanoseconds
     * since 1970-01-01T00:00:00Z.
     */
    template<typename T>
    double global_trigger_time(const T & sr) { return sr.hdr.triggerinfo.global_trigger_time; }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, global_trigger_time, global_trigger_time);

    /**
     * @brief Variable for the time of the beam gate in UTC
     * @details This variable returns the time of the beam gate in the absolute
     * time system in nanoseconds since 1970-01-01T00:00:00Z.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @return the time of the beam gate in nanoseconds since
     * 1970-01-01T00:00:00Z.
     */
    template<typename T>
    double beam_gate_time_abs(const T & sr) { return sr.hdr.triggerinfo.beam_gate_time_abs; }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, beam_gate_time_abs, beam_gate_time_abs);

    /**
     * @brief Variable for the time of the trigger within the beam gate.
     * @details This variable returns the time of the trigger within the beam
     * gate in microseconds.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @return the time of the trigger within the beam gate in microseconds.
     */
    template<typename T>
    double trigger_within_gate(const T & sr) { return sr.hdr.triggerinfo.trigger_within_gate; }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, trigger_within_gate, trigger_within_gate);

    /**
     * @brief Variable for the time of the beam gate in the detector time
     * system.
     * @details This variable returns the time of the beam gate in the
     * detector time system in microseconds.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @return the time of the beam gate in the detector time system in
     */
    template<typename T>
    double beam_gate_det_time(const T & sr) { return sr.hdr.triggerinfo.beam_gate_det_time; }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, beam_gate_det_time, beam_gate_det_time);

    /**
     * @brief Variable for the time of the global trigger in the detector time
     * system.
     * @details This variable returns the time of the global trigger in the
     * detector time system in microseconds.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @return the time of the global trigger in the detector time system in
     * microseconds.
     */
    template<typename T>
    double global_trigger_det_time(const T & sr) { return sr.hdr.triggerinfo.global_trigger_det_time; }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, global_trigger_det_time, global_trigger_det_time);

    /**
     * @brief Variable for time of the flash closest to the trigger time.
     * @details This variable is intended to provide the time of the flash
     * closest to the trigger time of the event. It is useful for producing a
     * "tophat"-style plot for locating the beam window and validating the
     * normalization.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @return the time of the flash closest to the trigger time.
     */
    template<typename T>
    double time_of_flash_closest_to_trigger(const T & sr)
    {
        double t0 = sr.hdr.triggerinfo.trigger_within_gate;
        double closest_flash_to_trigger = 10000;
        for(const auto & flash : sr.opflashes)
        {
            if(std::abs(flash.firsttime) < std::abs(closest_flash_to_trigger))
            {
                closest_flash_to_trigger = flash.firsttime;
            }
        }
        return closest_flash_to_trigger + t0;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, time_of_flash_closest_to_trigger, time_of_flash_closest_to_trigger);
}

#endif