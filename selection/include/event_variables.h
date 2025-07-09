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

    template<typename T>
    double is_first_in_subrun(const T & sr)
    {
        // This variable returns 1 if the event is the first in the subrun,
        // otherwise it returns 0.
        return sr.hdr.first_in_subrun;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, is_first_in_subrun, is_first_in_subrun);

    /**
     * @brief Variable for the POT (Protons on Target) in the event.
     * @details This variable retrieves the POT (Protons on Target) in the
     * event by attaching to the pot variable in the header of the record.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @return the POT in the event.
     */
    template<typename T>
    double pot(const T & sr) { return sr.hdr.pot; }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, pot, pot);

    /**
     * @brief Variable for the POT (Protons on Target) from the spillinfo
     * vector in the header of the record.
     * @details This variable retrieves the POT (Protons on Target) from
     * the spillinfo vector in the header of the record. It sums up the
     * TOR875 values from all the spills in the BNBInfo vector.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @param params the parameters for the cut. This is used to apply a scale
     * factor to the POT if needed.
     * @return the total POT from the spillinfo vector in the header of the record.
     */
    template<typename T>
    double pot_from_spillinfo(const T & sr, std::vector<double> params={1.0})
    {
        double pot = 0;
        for(const auto & spill : sr.hdr.bnbinfo)
        {
            pot += params.at(0)*spill.TOR875;
        }
        return pot;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, pot_from_spillinfo, pot_from_spillinfo);

    /**
     * @brief Variable for the number of generated events (MC only) in the
     * event.
     * @details This variable retrieves the number of generated events in the
     * event by attaching to the ngenevt variable in the header of the record.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @return the number of generated events in the event.
     */
    template<typename T>
    double ngenevt(const T & sr) { return sr.hdr.ngenevt; }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, ngenevt, ngenevt);

    /**
     * @brief Variable for the number of BNB spills in the event.
     * @details This variable counts the number of BNB spills in the event by
     * checking the length of the BNBInfo vector in the header of the record.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @return the number of BNB spills in the event.
     */
    template<typename T>
    double nbnb(const T & sr) { return sr.hdr.bnbinfo.size(); }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, nbnb, nbnb);

    /**
     * @brief Variable for the number of NuMI spills in the event.
     * @details This variable counts the number of NuMI spills in the event by
     * checking the length of the NuMIInfo vector in the header of the record.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @return the number of NuMI spills in the event.
     */
    template<typename T>
    double nnumi(const T & sr) { return sr.hdr.numiinfo.size(); }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, nnumi, nnumi);

    /**
     * @brief Variable for the number of off-beam BNB gates in the event.
     * @details This variable retrieves the number of off-beam BNB gates in the
     * event by attaching to the noffbeambnb variable in the header of the
     * record.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @return the number of off-beam BNB gates in the event.
     */
    template<typename T>
    double noffbeambnb(const T & sr) { return sr.hdr.noffbeambnb; }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, noffbeambnb, noffbeambnb);

    /**
     * @brief Variable for the number of off-beam NuMI gates in the event.
     * @details This variable retrieves the number of off-beam NuMI gates in
     * the event by attaching to the noffbeamnumi variable in the header of the
     * record.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @return the number of off-beam NuMI gates in the event.
     */
    template<typename T>
    double noffbeamnumi(const T & sr) { return sr.hdr.noffbeamnumi; }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, noffbeamnumi, noffbeamnumi);

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
     * @param params The offset to subtract from the time of the flash in the
     * minimization process. The default value is 0.0, which means no offset
     * @return the time of the flash closest to the trigger time.
     */
    template<typename T>
    double time_of_flash_closest_to_trigger(const T & sr, std::vector<double> params={0.0})
    {
        if(params.size() < 1)
        {
            throw std::runtime_error("time_of_flash_closest_to_trigger requires at least one parameter for the offset.");
        }
        double t0 = sr.hdr.triggerinfo.trigger_within_gate;
        double closest_flash_to_trigger = 10000;
        for(const auto & flash : sr.opflashes)
        {
            if(std::abs(flash.firsttime - params[0]) < std::abs(closest_flash_to_trigger - params[0]))
            {
                closest_flash_to_trigger = flash.firsttime;
            }
        }
        return closest_flash_to_trigger + t0;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, time_of_flash_closest_to_trigger, time_of_flash_closest_to_trigger);
}

#endif