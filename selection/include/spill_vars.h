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
#include <cmath>

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRBNBInfo.h"
#include "sbnanaobj/StandardRecord/SRNuMIInfo.h"
#include "sbnana/SBNAna/Vars/getBNBFoM.h"

#include "framework.h"

/**
 * @brief Guard a spill field against NaN/inf, substituting PLACEHOLDER.
 * @details Evaluates @p expr once as a double, then returns it unchanged if
 * finite or substitutes PLACEHOLDER if the value is NaN or infinite.
 * Intended for use as the sole statement in a spill variable function body.
 * @param expr the field expression to evaluate (e.g. @c spill.TOR875).
 */
#define SAFE_SPILL_VAR(expr)                                                \
    const double _v = static_cast<double>(expr);                            \
    return (std::isnan(_v) || std::isinf(_v)) ? PLACEHOLDERVALUE : _v

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
    double tor875(const T & spill) { SAFE_SPILL_VAR(spill.TOR875); }
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
    double tor860(const T & spill) { SAFE_SPILL_VAR(spill.TOR860); }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, tor860, tor860);

    /**
     * @brief Variable for the spill readout time (spill_time_sec).
     * @details The Unix timestamp of the spill readout in seconds. Used as
     * an input to the BNB Figure of Merit calculation to account for
     * time-dependent beam conditions.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the spill readout time in seconds.
     */
    template<typename T>
    double spill_time_sec(const T & spill) { SAFE_SPILL_VAR(spill.spill_time_sec); }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, spill_time_sec, spill_time_sec);

    /**
     * @brief Variable for the horizontal beam-width sigma at MWP 875 (M875HS).
     * @details M875HS is the horizontal beam-profile sigma measured by the
     * multi-wire profile monitor at station 875, in millimetres. It is
     * obtained by fitting a Gaussian to the wire-readout profile and
     * characterises the transverse beam size entering the target region.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the M875HS sigma in mm.
     */
    template<typename T>
    double m875hs(const T & spill) { SAFE_SPILL_VAR(spill.M875HS); }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, m875hs, m875hs);

    /**
     * @brief Variable for the vertical beam-width sigma at MWP 875 (M875VS).
     * @details M875VS is the vertical beam-profile sigma measured by the
     * multi-wire profile monitor at station 875, in millimetres. It is the
     * vertical counterpart to @ref m875hs.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the M875VS sigma in mm.
     */
    template<typename T>
    double m875vs(const T & spill) { SAFE_SPILL_VAR(spill.M875VS); }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, m875vs, m875vs);

    /**
     * @brief Variable for the horizontal beam-width sigma at MWP 876 (M876HS).
     * @details M876HS is the horizontal beam-profile sigma measured by the
     * multi-wire profile monitor at station 876, in millimetres. Station 876
     * is located immediately downstream of station 875 and provides a
     * complementary measurement of the beam width.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the M876HS sigma in mm.
     */
    template<typename T>
    double m876hs(const T & spill) { SAFE_SPILL_VAR(spill.M876HS); }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, m876hs, m876hs);

    /**
     * @brief Variable for the vertical beam-width sigma at MWP 876 (M876VS).
     * @details M876VS is the vertical beam-profile sigma measured by the
     * multi-wire profile monitor at station 876, in millimetres. It is the
     * vertical counterpart to @ref m876hs.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the M876VS sigma in mm.
     */
    template<typename T>
    double m876vs(const T & spill) { SAFE_SPILL_VAR(spill.M876VS); }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, m876vs, m876vs);

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
    double thcurr(const T & spill) { SAFE_SPILL_VAR(spill.THCURR); }
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
    double lm875a(const T & spill) { SAFE_SPILL_VAR(spill.LM875A); }
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
    double lm875b(const T & spill) { SAFE_SPILL_VAR(spill.LM875B); }
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
    double lm875c(const T & spill) { SAFE_SPILL_VAR(spill.LM875C); }
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
    double hp875(const T & spill) { SAFE_SPILL_VAR(spill.HP875); }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, hp875, hp875);

    /**
     * @brief Variable for the horizontal beam position at Target Station 1
     * (HPTG1).
     * @details HPTG1 is the horizontal position monitor at Target Station 1,
     * upstream of the beryllium target, in millimetres. It is the upstream
     * counterpart to @ref hptg2 and is used together with it in the BNB
     * Figure of Merit calculation.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the HPTG1 reading in mm.
     */
    template<typename T>
    double hptg1(const T & spill) { SAFE_SPILL_VAR(spill.HPTG1); }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, hptg1, hptg1);

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
    double vp875(const T & spill) { SAFE_SPILL_VAR(spill.VP875); }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, vp875, vp875);

    /**
     * @brief Variable for the vertical beam position near Mag 873 (VP873).
     * @details VP873 is a vertical position monitor located after Magnet 873,
     * in millimetres. It serves as a reliable replacement for VPTG1 and VPTG2
     * which are known to be unreliable for ICARUS and SBND runs, and is used
     * in the BNB Figure of Merit calculation.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the VP873 reading in mm.
     */
    template<typename T>
    double vp873(const T & spill) { SAFE_SPILL_VAR(spill.VP873); }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, vp873, vp873);

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
    double hptg2(const T & spill) { SAFE_SPILL_VAR(spill.HPTG2); }
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
    double vptg2(const T & spill) { SAFE_SPILL_VAR(spill.VPTG2); }
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
    double btjt2(const T & spill) { SAFE_SPILL_VAR(spill.BTJT2); }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, btjt2, btjt2);

    /**
     * @brief Computed BNB Figure of Merit (FOM) from individual beam monitors.
     * @details Computes the FOM from the individual beam monitoring device
     * readings using getBNBFoM(), taking into account both position/intensity
     * monitors and multi-wire beam-width measurements. This is the equivalent
     * of kSpillFoM in sbnana/SBNAna/Vars/BNBVars.cxx.
     * @tparam T the BNB spill container type.
     * @param spill the BNB spill info object to apply the variable on.
     * @return the computed FOM value for the spill.
     */
    template<typename T>
    double fom(const T & spill) {
        return getBNBFoM(
            spill_time_sec(spill),
            tor860(spill),
            tor875(spill),
            hp875(spill),
            hptg1(spill),
            hptg2(spill),
            vp873(spill),
            vp875(spill),
            m875hs(spill),
            m875vs(spill),
            m876hs(spill),
            m876vs(spill)
        );
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BNBSpill, fom, fom);

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
    double trtgtd(const T & spill) { SAFE_SPILL_VAR(spill.TRTGTD); }
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
    double event(const T & spill) { SAFE_SPILL_VAR(spill.event); }
    REGISTER_VAR_SCOPE(RegistrationScope::BothSpill, event, event);

} // namespace svar
#endif // SPILL_VARS_H
