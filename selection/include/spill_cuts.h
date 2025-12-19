/**
 * @file spill_cuts.h
 * @brief Definitions of analysis cuts which can be applied to individual
 * spills.
 * @details This file contains definitions of analysis cuts which can be
 * applied to the individual spills of the StandardRecord object. Note that
 * these are NOT "SpillCuts" in the parlance of CAFAna (which technically act
 * on entire events), but rather cuts which are applied to the *actual* spill
 * info in the StandardRecord object. Each cut is implemented as a function
 * which takes a StandardRecord object as an argument and returns a boolean
 * indicating whether the spill passes the cut.
 * @author mueller@fnal.gov
 */
#ifndef SPILL_CUTS_H
#define SPILL_CUTS_H
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/SBNAna/Vars/getBNBFoM.h"

#include "framework.h"

/**
 * @namespace scut
 * @brief Namespace for organizing cuts which act on spills.
 * @details This namespace is intended to be used for organizing cuts
 * which act on spills. Each cut is implemented as a function which takes
 * a StandardRecord object as an argument and returns a boolean.
 */
namespace scut
{
    /**
     * @brief Apply no cut to the spill.
     * @details This cut always returns true, indicating that the spill passes
     * the cut. It is useful for testing purposes or when no cut is desired.
     * @tparam T the spill container.
     * @param sr the StandardRecord to apply the cut on.
     * @return true, indicating that the spill passes the cut.
     */
    template<typename T>
    bool no_cut(const T & spill) { return true; }
    REGISTER_CUT_SCOPE(RegistrationScope::Spill, no_cut, no_cut);

    /**
     * @brief Apply a cut on various spill properties.
     * @details This cut checks if the spill properties are within some
     * interval of values centered at "nominal" values. This is intended to cut
     * out spills with low quality, and therefore with low relevance to
     * assumptions made in the analysis / central value simulation. This must
     * only be applied on data.
     * @tparam T the spill container type
     * @param sr the StandardRecord to apply the cut on.
     * @return true if the spill properties are within the specified intervals,
     * false otherwise.
     */
    template<typename T>
    bool beam_quality_cut(const T & spill)
    {
        // Sanity check to make sure no NaNs are present in the spill
        // properties.
        if(std::isnan(spill.TOR860) || std::isnan(spill.TOR875) ||
           std::isnan(spill.LM875A) || std::isnan(spill.LM875B) ||
           std::isnan(spill.LM875C) || std::isnan(spill.THCURR))
        {
            return false; // Spill properties must not be NaN
        }

        // Check if the spill properties are within the specified intervals.
        return (
            spill.TOR860 > +100e9 &&
            spill.TOR875 > +100e9 &&
            spill.LM875A > +1e-2 &&
            spill.LM875B > +1e-2 &&
            spill.LM875C > +1e-2 &&
            spill.THCURR > +173 && spill.THCURR < +175
        );
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Spill, beam_quality_cut, beam_quality_cut);

    /**
     * @brief Apply a cut on the beam FoM value of the spill.
     * @details This cut checks if the beam FoM value of the spill is above a
     * certain threshold, which can be configured by the user. This is a
     * wrapper around the tagged version implemented by the getBNBFoM 
     * function in `sbnana`. It is intended to be a handle on the overalap of
     * the beam with the target.
     * @tparam T the spill container type
     * @param sr the StandardRecord to apply the cut on.
     * @param params a vector of parameters for the cut, which can be used to
     * specify the threshold for the FoM value.
     * @return true if the FoM value is above the threshold, false otherwise.
     */
    template<typename T>
    bool bnb_fom_cut(const T & spill, std::vector<double> params={})
    {
        if(params.empty())
        {
            throw std::invalid_argument("bnb_fom_cut requires at least one parameter for the threshold (recommended 0.98).");
        }
        double threshold = params[0];
        return (getBNBFoM(
            spill.spill_time_sec,
            spill.TOR860,
            spill.TOR875,
            spill.HP875,
            spill.HPTG1,
            spill.HPTG2,
            spill.VP873,
            spill.VP875,
            spill.M875HS,
            spill.M875VS,
            spill.M876HS,
            spill.M876VS
        ) >= threshold);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Spill, bnb_fom_cut, bnb_fom_cut);

} // namespace scut
#endif // SPILL_CUTS_H