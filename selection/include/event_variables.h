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
#include "sbnanaobj/StandardRecord/SRBNBInfo.h"
#include "sbnana/SBNAna/Vars/BNBVars.h"

#include "framework.h"
#include "utilities.h"

/**
 * @brief Global vector to store BNB information across events.
 * @details This vector is used to store the BNB information across events,
 * which is necessary for certain calculations because the spill information is
 * only stored for the first event in the subrun. The vector is cleared when
 * the first event in the subrun is encountered, and it is filled with the
 * event number and the TOR875 value from the BNBInfo vector in the header of
 * the record.
 */
std::vector<std::tuple<uint32_t, double>> global_bnb_info;
size_t global_bnb_event_number = 0;

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
     * @brief Variable for the multiplicity of neutrino interactions in the
     * event.
     * @details This variable counts the number of neutrino interactions in the
     * event by checking how many interactions have a neutrino ID greater than
     * -1 (equivalent to the cuts::neutrino cut).
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @return double the multiplicity of neutrino interactions in the event.
     */
    template<typename T>
    double nnu(const T & sr)
    {
        size_t count = 0;
        for(const auto & interaction : sr.dlp_true)
        {
            if(interaction.nu_id > -1) ++count;
        }
        return count;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, nnu, nnu);

    /**
     * @brief Variable for the multiplicity of in-time interactions in the
     * event.
     * @details This variable counts the number of in-time interactions in the
     * event by checking how many interactions have a particle with a time
     * within the beam gate (i.e. the interaction creates activity in the beam
     * gate).
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @param params The beam gate window in microseconds. The default is
     * [0.0, 1.6].
     * @return double the multiplicity of in-time interactions in the event.
     */
    template<typename T>
    double nintime(const T & sr, std::vector<double> params={0.0, 1.6})
    {
        size_t count = 0;
        for(const auto & interaction : sr.dlp_true)
        {
            for(const auto & p : interaction.particles)
            {
                if(p.t >= params[0] && p.t <= params[1])
                {
                    ++count;
                    break; // Only count the interaction once
                }
            }
        }
        return count;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, nintime, nintime);

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
    double pot_from_spillinfo(const T & sr, std::vector<double> params={})
    {
        if(params.size() < 1)
            params.push_back(1.0); // Default scale factor if not provided
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
     * normalization. This variable uses the `firsttime` field of the optical
     * flash. The parameterized offset is used to account for the natural
     * offset of the reconstructed flash time from the trigger time, which is
     * not zero despite all systems being referenced to the trigger.
     * 
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
        size_t closest_flash_index = utilities::first_opflash_firsttime(sr, params[0]);
        if(closest_flash_index == kNoMatch)
            return kNoMatchValue;
        else
            return sr.opflashes[closest_flash_index].firsttime + t0;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, time_of_flash_closest_to_trigger, time_of_flash_closest_to_trigger);

    /**
     * @brief Variable for time of the flash closest to the trigger time.
     * @details This variable is intended to provide the time of the flash
     * closest to the trigger time of the event. It is useful for producing a
     * "tophat"-style plot for locating the beam window and validating the
     * normalization. This version uses the raw time of the flash instead of
     * the 'firsttime' field. The parameterized offset is used to account for
     * the natural offset of the reconstructed flash time from the trigger
     * time, which is not zero despite all systems being referenced to the
     * trigger.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @param params The offset to subtract from the time of the flash in the
     * minimization process. The default value is 0.0, which means no offset
     * @return the time of the flash closest to the trigger time.
     */
    template<typename T>
    double time_of_flash_closest_to_trigger_rawtime(const T & sr, std::vector<double> params={0.0})
    {
        if(params.size() < 1)
        {
            throw std::runtime_error("time_of_flash_closest_to_trigger_rawtime requires at least one parameter for the offset.");
        }
        double t0 = sr.hdr.triggerinfo.trigger_within_gate;
        size_t closest_flash_index = utilities::first_opflash_rawtime(sr, params[0]);
        if(closest_flash_index == kNoMatch)
            return kNoMatchValue;
        else
            return sr.opflashes[closest_flash_index].time + t0;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, time_of_flash_closest_to_trigger_rawtime, time_of_flash_closest_to_trigger_rawtime);

    /**
     * @brief Variable (wrapper) for the FoM2 (Figure of Merit 2) in the event.
     * @details This variable is a wrapper for the FoM2 variable, which is
     * defined as a SpillVar in the usual CAFAna parlance. It is used as a
     * metric that roughly characterizes the overlap of the beam with the
     * target and can be used as a cut to reject events that correspond to bad
     * beam conditions.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @return the FoM2 value for the event.
     */
    template<typename T>
    double bnb_fom2(const T & sr)
    {
        return ana::kSpillFoM2(&sr);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, bnb_fom2, bnb_fom2);

    /**
     * @brief Variable for the unfolded event POT (Protons on Target).
     * @details This variable retrieves the unfolded event POT by summing up
     * the TOR875 values from the BNBInfo vector in the header of the record.
     * This uses the stored global BNB info to ensure that the POT is bookkept
     * correctly for each event.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @return the unfolded event POT in the event.
     */
    template<typename T>
    double unfolded_event_pot(const T & sr)
    {
        // If this is the first event in the subrun, we need to reset the
        // global BNB info.
        if(sr.hdr.first_in_subrun && global_bnb_event_number != sr.hdr.evt)
        {
            global_bnb_info.clear();
            for(const auto & bnb_info : sr.hdr.bnbinfo)
            {
                // Store the event number and the TOR875 value.
                global_bnb_info.emplace_back((uint32_t)bnb_info.event, (double)bnb_info.TOR875);
            }
            global_bnb_event_number = sr.hdr.evt;
        }

        // Loop over the global BNB info and filter out the events that are not
        // this one.
        double pot = 0.0;
        for(const auto & bnb_info : global_bnb_info)
        {
            if(std::get<0>(bnb_info) == sr.hdr.evt)
            {
                // Add the TOR875 value to the pot.
                pot += std::get<1>(bnb_info);
            }
        }

        return pot;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, unfolded_event_pot, unfolded_event_pot);

    /**
     * @brief Variable for the number of unfolded BNB events in the event.
     * @details This variable counts the number of unfolded BNB events in the
     * event by checking the global BNB info vector. It is used to ensure that
     * the number of BNB events is bookkept correctly for each event.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @return the number of unfolded BNB events in the event.
     */
    template<typename T>
    double unfolded_event_nbnb(const T & sr)
    {
        // If this is the first event in the subrun, we need to reset the
        // global BNB info.
        if(sr.hdr.first_in_subrun && global_bnb_event_number != sr.hdr.evt)
        {
            global_bnb_info.clear();
            for(const auto & bnb_info : sr.hdr.bnbinfo)
            {
                // Store the event number and the TOR875 value.
                global_bnb_info.emplace_back((uint32_t)bnb_info.event, (double)bnb_info.TOR875);
            }
            global_bnb_event_number = sr.hdr.evt;
        }

        // Loop over the global BNB info and filter out the events that are not
        // this one.
        size_t nbnbs = 0;
        for(const auto & bnb_info : global_bnb_info)
        {
            if(std::get<0>(bnb_info) == sr.hdr.evt)
            {
                // Count the number of BNB events.
                ++nbnbs;
            }
        }

        return nbnbs;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Event, unfolded_event_nbnb, unfolded_event_nbnb);
}

#endif