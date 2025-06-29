/**
 * @file spill.h
 * @brief Definitions of analysis variables which can extract information from
 * the StandardRecord object.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from the StandardRecord object. Each variable
 * is implemented as a function which takes a StandardRecord object as an
 * argument and returns a double.
 * @author mueller@fnal.gov
 */
#ifndef SPILL_H
#define SPILL_H
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "framework.h"

/**
 * @namespace event
 * @brief Namespace for organizing variables which act on events.
 * @details This namespace is intended to be used for organizing variables
 * which act on events. Each variable is implemented as a function which takes
 * a StandardRecord object as an argument and returns a double.
 */
namespace event
{
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
            if(std::abs(flash.firsttime - t0) < std::abs(closest_flash_to_trigger - t0))
            {
                closest_flash_to_trigger = flash.firsttime;
            }
        }
        return closest_flash_to_trigger + t0;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, time_of_flash_closest_to_trigger, time_of_flash_closest_to_trigger);

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
    double global_trigger_time(const T & sr)
    {
        return sr.hdr.triggerinfo.global_trigger_time;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, global_trigger_time, global_trigger_time);
}

#endif