/**
 * @file cuts_pi0ana.h
 * @brief Header file for definitions of analysis cuts specific neutral pion-based analyses.
 * @details This file contains definitions of analysis cuts which can be used
 * to select interactions specific to the pi0ana analysis. The cuts are
 * intended to be used in conjunction with the generic cuts defined in cuts.h.
 * @author lkashur@colostate.edu
*/
#ifndef CUTS_PI0ANA_H
#define CUTS_PI0ANA_H
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "utilities_pi0ana.h"
#include "framework.h"
#include "selectors.h"

/**
 * @namespace cuts::pi0ana
 * @brief Namespace for organizing cuts specific to the ccpi0ana analysis.
 * @details This namespace is intended to be used for organizing cuts which act
 * on interactions specific to the ccpi0ana analysis. Each cut is implemented as
 * a function which takes an interaction object as an argument and returns a
 * boolean. The function should be templated on the type of interaction object if
 * the cut is intended to be used on both true and reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the cuts
 * namespace, which is used for organizing generic cuts which act on interactions.
 */
namespace cuts::pi0ana
{

  /**
   * @brief Apply cut on the number of photons in interaction (at least two).
   * @details This function applies a cut on the number of above-threshold photons
   * in the interaction, requiring at least two.
   * @tparam T the type of interaction (true or reco).
   * @param obj the interaction to select on.
   * @return true if the interaction contains at least two above-threshold photons.
   */
    template<class T>
    bool at_least_two_photons(const T & obj, std::vector<double> params={})
    {
        size_t count(0);
        for(const auto & p : obj.particles)
	      {
	          if(pvars::pid(p) == 0 && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
	          ++count;
	      }
        return count >= 2;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, at_least_two_photons, at_least_two_photons);

    /**
     * @brief Apply cut on the number of photons in interaction (less than four).
     * @details This function applies a cut on the number of above-threshold photons
     * in the interaction, requiring less than four.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction contains less than four above-threshold photons.
     */
    template<class T>
    bool less_than_four_photons(const T & obj, std::vector<double> params={})
    {
        size_t count(0);
        for(const auto & p : obj.particles)
	      {
	          if(pvars::pid(p) == 0 && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
	          ++count;
	      }
        return count < 4;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, less_than_four_photons, less_than_four_photons);

    /**
     * @brief Apply cut on the interaction's leading photon kinetic energy.
     * @details This function applies a cut on the interaction's leading photon,
     * requiring its kinetic energy to be above a specified threshold.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params the kinetic enregy threshold applied.
     * @return true if the interaction contains a photon with kinetic energy above specifed theshold.
     */
    template<class T>
    bool leading_photon_ke_cut(const T & obj, std::vector<double> params={})
    {
      
        bool passes(false);
        size_t phi = selectors::leading_photon(obj);
	      if(phi == kNoMatch)
	      {
	          return false;
	      }
	      else
	      {
	          auto & ph(obj.particles[phi]);
	          return pvars::calo_ke(ph) >= params[0];
	      }
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, leading_photon_ke_cut, leading_photon_ke_cut);
    
    /**
     * @brief Apply pi0 mass cut.
     * @details This function applies a cut on the invariant diphoton mass
     * of the interactions two most energetic photons.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the pi0 mass cut.
     */
    template<class T>
    bool valid_pi0_mass_cut(const T & obj, std::vector<double> params = {})
    {
        pi0 s = utilities_pi0ana::pi0_info(obj);
	      return s.mass < params[0];
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, valid_pi0_mass_cut, valid_pi0_mass_cut);

    /**
     * @brief Apply single pi0 cut in truth
     * @details This function applies a cut on the number of 
     * of primary neutral pions in the interaction at the truth level.
     * @param obj the kinetic energy threshold for the neutral pion.
     * @return true if number of true primary neutral pions above threshold is one.
     */
    template<class T>
    bool single_pi0(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params = {0.0,})
    {
        double num_primary_pi0s = utilities_pi0ana::true_primary_pi0_multiplicity(obj, params);
	      return num_primary_pi0s == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::True, single_pi0, single_pi0);
}
#endif // CUTS_PI0ANA_H
