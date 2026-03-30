/**
 * @file spill_vars.h
 * @brief Definitions of analysis variables which can be applied to individual
 * spills.
 * @details This file contains definitions of analysis variables which can be
 * applied to the individual spills stored in the StandardRecord header. Note
 * that these are NOT "SpillVars" in the parlance of CAFAna (which technically
 * act on entire events), but rather variables extracted from the *actual* spill
 * info objects (SRBNBInfo, SRNuMIInfo) stored in the header. Each variable is
 * implemented as a function which takes a spill info object as an argument and
 * returns a double.
 * @author mueller@fnal.gov
 */
#ifndef SPILL_VARS_H
#define SPILL_VARS_H
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRBNBInfo.h"
#include "sbnanaobj/StandardRecord/SRNuMIInfo.h"

#include "framework.h"

/**
 * @namespace svar
 * @brief Namespace for organizing variables which act on individual spills.
 * @details This namespace is intended to be used for organizing variables
 * which act on individual BNB or NuMI spill info objects. Each variable is
 * implemented as a function which takes a spill info object as an argument
 * and returns a double. BNB variables are registered under the
 * @ref RegistrationScope::BNBSpill scope; NuMI variables under
 * @ref RegistrationScope::NuMISpill.
 */
namespace svar
{
    //=========================================================================
    // BNB spill variables
    //=========================================================================

    /**
     * @brief Variable for the primary POT toroid (TOR875).
     * @details TOR875 is the primary toroid used for POT accounting in the BNB
     * beamline. It is located downstream of the focusing horn, close to the
     * target, and is the canonical quantity used when normalizing BNB data
     * samples to a POT exposure.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the TOR875 reading in protons.
     */
    template<typename T>
    double tor875(const T & spill) { return spill.TOR875; }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, tor875, tor875);

    /**
     * @brief Variable for the upstream POT toroid (TOR860).
     * @details TOR860 is located upstream of TOR875 and provides a secondary
     * measurement of the integrated proton beam intensity. Comparing TOR860
     * and TOR875 is useful for diagnosing beam losses between the two
     * monitors.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the TOR860 reading in protons.
     */
    template<typename T>
    double tor860(const T & spill) { return spill.TOR860; }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, tor860, tor860);

    /**
     * @brief Variable for the BNB Figure of Merit (FOM).
     * @details The FOM is a composite beam-quality metric for the BNB. See
     * [SBN DocDB 41901](https://sbn-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=41901)
     * for the definition. It is useful for quantifying the overall quality of
     * a spill beyond the individual monitor readings.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the FOM value for the spill.
     */
    template<typename T>
    double fom(const T & spill) { return spill.FOM; }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, fom, fom);

    /**
     * @brief Variable for the horn current (THCURR).
     * @details THCURR is the current applied to the focusing horn in
     * kiloAmperes. The nominal value for neutrino-mode running is +174 kA;
     * deviations from this value indicate abnormal horn conditions that can
     * alter the neutrino flux.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the horn current in kiloAmperes.
     */
    template<typename T>
    double thcurr(const T & spill) { return spill.THCURR; }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, thcurr, thcurr);

    /**
     * @brief Variable for the first loss monitor before the RWM (LM875A).
     * @details LM875A measures beam losses immediately before the resistive
     * wall monitor (RWM) in units of R/s. Elevated readings can indicate beam
     * halo or scraping and are used as part of the beam-quality selection.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the LM875A reading in R/s.
     */
    template<typename T>
    double lm875a(const T & spill) { return spill.LM875A; }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, lm875a, lm875a);

    /**
     * @brief Variable for the second loss monitor after the RWM (LM875B).
     * @details LM875B measures beam losses immediately after the resistive
     * wall monitor (RWM) in units of R/s. Together with LM875A and LM875C
     * it characterises the beam halo around the target.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the LM875B reading in R/s.
     */
    template<typename T>
    double lm875b(const T & spill) { return spill.LM875B; }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, lm875b, lm875b);

    /**
     * @brief Variable for the third loss monitor after the RWM (LM875C).
     * @details LM875C provides a third measurement of beam losses after the
     * RWM in units of R/s. It is used alongside LM875A and LM875B in the
     * beam-quality cut.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the LM875C reading in R/s.
     */
    template<typename T>
    double lm875c(const T & spill) { return spill.LM875C; }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, lm875c, lm875c);

    /**
     * @brief Variable for the horizontal beam position near Mag 875 (HP875).
     * @details HP875 is a horizontal position monitor located after Magnet
     * 875, in millimetres. It characterises the transverse centering of the
     * beam on the target and is useful for studies of beam-position-dependent
     * flux variations.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the HP875 reading in mm.
     */
    template<typename T>
    double hp875(const T & spill) { return spill.HP875; }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, hp875, hp875);

    /**
     * @brief Variable for the vertical beam position near Mag 875 (VP875).
     * @details VP875 is a vertical position monitor located after Magnet 875,
     * in millimetres. It is the vertical counterpart to @ref hp875 and is
     * used together with it to characterise the 2D beam spot position.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the VP875 reading in mm.
     */
    template<typename T>
    double vp875(const T & spill) { return spill.VP875; }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, vp875, vp875);

    /**
     * @brief Variable for the horizontal beam position at Target Station 2
     * (HPTG2).
     * @details HPTG2 is the horizontal position monitor at Target Station 2,
     * which is the station closest to the beryllium target, in millimetres.
     * It provides the most direct measurement of where the proton beam
     * strikes the target and is therefore the most useful position monitor
     * for flux systematics studies.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the HPTG2 reading in mm.
     */
    template<typename T>
    double hptg2(const T & spill) { return spill.HPTG2; }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, hptg2, hptg2);

    /**
     * @brief Variable for the vertical beam position at Target Station 2
     * (VPTG2).
     * @details VPTG2 is the vertical position monitor at Target Station 2,
     * the station closest to the beryllium target, in millimetres. It is the
     * vertical counterpart to @ref hptg2.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the VPTG2 reading in mm.
     */
    template<typename T>
    double vptg2(const T & spill) { return spill.VPTG2; }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, vptg2, vptg2);

    /**
     * @brief Variable for the target air temperature (BTJT2).
     * @details BTJT2 is the temperature of the air exiting the target
     * assembly in degrees Celsius. Anomalous temperature readings can
     * indicate target overheating or cooling-system faults and may be used
     * as a diagnostic for target integrity.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the BTJT2 reading in degrees Celsius.
     */
    template<typename T>
    double btjt2(const T & spill) { return spill.BTJT2; }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, btjt2, btjt2);

    //=========================================================================
    // NuMI spill variables
    //=========================================================================

    /**
     * @brief Variable for the primary POT toroid for NuMI (TRTGTD).
     * @details TRTGTD is the primary toroid used for POT accounting in the
     * NuMI beamline, located at the target. It is the canonical quantity used
     * when normalizing NuMI data samples to a POT exposure, analogous to
     * TOR875 for BNB.
     * @tparam T the NuMI spill container type.
     * @param spill the NuMI spill info object to apply the variable on.
     * @return the TRTGTD reading in protons.
     */
    template<typename T>
    double trtgtd(const T & spill) { return spill.TRTGTD; }
    REGISTER_VAR_SCOPE(RegistrationScope::NuMISpill, trtgtd, trtgtd);

    //=========================================================================
    // Shared spill variables
    //=========================================================================

    /**
     * @brief Variable for the event number matched to the spill.
     * @details Each spill is associated with a single event by grouping the
     * spills between the previous and current event boundaries. The spill
     * directly overlapping with the even is associated with that event.
     * @tparam T the spill container type (BNB or NuMI).
     * @param spill the spill info object to apply the variable on.
     * @return the event number associated with the spill.
     */
    template<typename T>
    double event(const T & spill) { return spill.event; }
    REGISTER_VAR_SCOPE(RegistrationScope::BothSpill, event, event);

} // namespace svar
#endif // SPILL_VARS_H
