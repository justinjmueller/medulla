/**
 * @file mctruth_cuts.h
 * @brief Definitions of analysis cuts applied at the GENIE generator level.
 * @details This file contains definitions of analysis cuts which operate
 * directly on the SRTrueInteraction object. Each cut is implemented as a
 * function which takes an SRTrueInteraction object as an argument and returns
 * a bool. The association of an SRInteractionTruthDLP object to an
 * SRTrueInteraction object is handled upstream in the framework.
 * @author mueller@fnal.gov
 * @author rvizarr@fnal.gov
 */
#ifndef MCTRUTH_CUTS_H
#define MCTRUTH_CUTS_H
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"

#include "framework.h"
#include "variables.h"

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

    /**
     * @brief Cut for exactly one true final state muon above threshold.
     * @details Applied at the GENIE generator level using obj.prim. The
     * kinetic energy is computed from the GENIE genE field.
     * @tparam T the type of the object to apply the cut on.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params KE threshold in MeV, defaults to 143.425 MeV.
     * @return true if exactly one muon above threshold.
     */
    template<typename T>
    bool single_muon(const T & obj, std::vector<double> params={143.425,})
    {
        int num_muons(0);
        for(const auto & p : obj.prim)
        {
            if(abs(p.pdg) == 13)
            {
                double ke = 1000. * (p.genE - (MUON_MASS/1000.));
                if(ke >= params[0])
                    num_muons++;
            }
        }
        return num_muons == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, single_muon, single_muon);

    /**
     * @brief Cut for zero true final state charged pions above threshold.
     * @details Applied at the GENIE generator level using obj.prim. The
     * kinetic energy is computed from the GENIE genE field.
     * @tparam T the type of the object to apply the cut on.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params KE threshold in MeV, defaults to 25 MeV.
     * @return true if no charged pions above threshold.
     */
    template<typename T>
    bool no_charged_pions(const T & obj, std::vector<double> params={25.0,})
    {
        for(const auto & p : obj.prim)
        {
            if(abs(p.pdg) == 211)
            {
                double ke = 1000. * (p.genE - (PION_MASS/1000.));
                if(ke >= params[0])
                    return false;
            }
        }
        return true;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_charged_pions, no_charged_pions);

    /**
     * @brief Cut for zero true final state neutral pions.
     * @details Applied at the GENIE generator level using obj.prim. Only
     * primary particles (start_process == 0) are considered.
     * @tparam T the type of the object to apply the cut on.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @return true if no neutral pions are present.
     */
    template<typename T>
    bool no_neutral_pions(const T & obj)
    {
        for(const auto & p : obj.prim)
        {
            if(p.start_process != 0) continue;
            if(p.pdg == 111)
                return false;
        }
        return true;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_neutral_pions, no_neutral_pions);

    /**
     * @brief Cut for zero true final state photons above threshold.
     * @details Applied at the GENIE generator level using obj.prim. The
     * energy is computed from the GENIE genE field.
     * @tparam T the type of the object to apply the cut on.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params energy threshold in MeV, defaults to 25 MeV.
     * @return true if no photons above threshold.
     */
    template<typename T>
    bool no_photons(const T & obj, std::vector<double> params={25.0,})
    {
        for(const auto & p : obj.prim)
        {
            if(p.pdg == 22)
            {
                double energy = 1000. * p.genE;
                if(energy >= params[0])
                    return false;
            }
        }
        return true;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_photons, no_photons);

    /**
     * @brief Cut for zero true final state electrons above threshold.
     * @details Applied at the GENIE generator level using obj.prim. The
     * kinetic energy is computed from the GENIE genE field.
     * @tparam T the type of the object to apply the cut on.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params KE threshold in MeV, defaults to 25 MeV.
     * @return true if no electrons above threshold.
     */
    template<typename T>
    bool no_electrons(const T & obj, std::vector<double> params={25.0,})
    {
        for(const auto & p : obj.prim)
        {
            if(abs(p.pdg) == 11)
            {
                double ke = 1000. * (p.genE - (ELECTRON_MASS/1000.));
                if(ke >= params[0])
                    return false;
            }
        }
        return true;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_electrons, no_electrons);

    /**
     * @brief Cut for exactly one true final state proton above threshold.
     * @details Applied at the GENIE generator level using obj.prim. The
     * kinetic energy is computed from the GENIE genE field.
     * @tparam T the type of the object to apply the cut on.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params KE threshold in MeV, defaults to 50 MeV.
     * @return true if exactly one proton above threshold.
     */
    template<typename T>
    bool single_proton(const T & obj, std::vector<double> params={50.0,})
    {
        int count(0);
        for(const auto & p : obj.prim)
        {
            if(p.pdg == 2212)
            {
                double ke = 1000. * (p.genE - (PROTON_MASS/1000.));
                if(ke >= params[0])
                    count++;
            }
        }
        return count == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, single_proton, single_proton);

    /**
     * @brief Apply a cut for a maximum energy transfer
     * @details This function applies a cut to select interactions below
     * a maximum energy transfer (aka, omega or q0).This cut is on the
     * true energy transfer from the neutrino to the hadronic system.
     * This is defined in the lab frame.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this is a
     * double representing an energy transfer in GeV
     * @return true if the energy transfer is below the given value
     */
    template<class T>
    bool is_below_max_e_transfer(const T & obj, std::vector<double> params={0.1,})
    {
      return obj.q0_lab<params[0];
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, is_below_max_e_transfer, is_below_max_e_transfer);

} // namespace mctruth
#endif
