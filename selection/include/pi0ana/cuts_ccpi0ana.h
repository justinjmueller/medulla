/**
 * @file cuts_ccpi0ana.h
 * @brief Header file for definitions of analysis cuts specific to the muonana
 * analysis.
 * @details This file contains definitions of analysis cuts which can be used
 * to select interactions specific to the pi0ana analysis. The cuts are
 * intended to be used in conjunction with the generic cuts defined in cuts.h.
 * @author lkashur@colostate.edu
*/
#ifndef CUTS_CCPI0ANA_H
#define CUTS_CCPI0ANA_H
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "utilities_ccpi0ana.h"
#include "framework.h"

/**
 * @namespace cuts::ccpi0ana
 * @brief Namespace for organizing cuts specific to the ccpi0ana analysis.
 * @details This namespace is intended to be used for organizing cuts which act
 * on interactions specific to the ccpi0ana analysis. Each cut is implemented as
 * a function which takes an interaction object as an argument and returns a
 * boolean. The function should be templated on the type of interaction object if
 * the cut is intended to be used on both true and reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the cuts
 * namespace, which is used for organizing generic cuts which act on interactions.
 */
namespace cuts::ccpi0ana
{
    /**
     * @brief Apply data cut.
     * @details This cut asserts the interaction is from data (not MC).
     * @param sr the Standard Record.
     * @return true if interaction is from data.
     * @note This cut is intended to be used for the ccpi0ana analysis.
     */
  //bool is_data(const caf::Proxy<caf::StandardRecord> sr) 
  //{
  //bool _is_data(sr.ndlp_true == 0);
  //return _is_data;
  //}
  
    /**
     * @brief Apply a 1mu0pi2gamma topological (final state) cut.
     * @details The interaction must have a topology matching 1mu0pi2gamma as defined by
     * the conditions in the @ref count_primaries() function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a 1mu0pi2gamma topology.
     * @note This cut is intended to be used for the ccpi0ana analysis.
     */
    template<class T>
        bool topological_1mu0pi2gamma_cut(const T & obj)
        {
	    std::vector<uint32_t> c(utilities_ccpi0ana::count_primaries(obj));
	    reco_inter s = utilities_ccpi0ana::reco_interaction_info(obj);
            return (utilities_ccpi0ana::reco_shower_criteria(obj) && s.pi0_momentum_mag >= MIN_PI0_MOMENTUM && c[2] == 1 && c[3] == 0);
        }

    template<class T>
        bool base_topology_cut(const T & obj)
        {
	    std::vector<uint32_t> c(utilities_ccpi0ana::count_primaries(obj));
	    return c[2] == 1 && c[3] == 0 && c[0] >= 2 && c[0] < 4;
	}
    REGISTER_CUT_SCOPE(RegistrationScope::Both, base_topology_cut, base_topology_cut);

    template<class T>
      bool leading_shower_energy_cut(const T & obj)
      {
	// default
	bool passes(false);

	size_t leading_shower_index(0);
	double max_shower_ke(-99999);

	// First loop to find leading shower
	for(size_t i(0); i < obj.particles.size(); ++i)
	  {
	    const auto & p = obj.particles[i];

	    // Primary particles                                                                                                                                                                       
	    if(!p.is_primary) continue;

	    // Leading shower
	    //if((PIDFUNC(p) == 0 || PIDFUNC(p) == 1) && p.ke > max_shower_ke0) // showers
	    if((PIDFUNC(p) == 0) && p.ke > max_shower_ke) // photons                                                                                                                                   
	      {
		max_shower_ke = p.ke;
		leading_shower_index = i;
	      }
	  }                                                                                                                                                                            

	if(max_shower_ke >= MIN_LEADING_SHOWER_ENERGY) passes = true;
	return passes;
      }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, leading_shower_energy_cut, leading_shower_energy_cut);

    template<class T>
        bool zero_charged_pions_cut(const T & obj)
        {
	    std::vector<uint32_t> c(utilities_ccpi0ana::count_primaries(obj));
	    return c[3] == 0;
	}

    template<class T>
        bool no_hadronic_activity_cut(const T & obj)
        {
	    bool passes(true);
	    TVector3 vertex(obj.vertex[0], obj.vertex[1], obj.vertex[2]);
	    // Loop over particle
	    for(size_t i(0); i < obj.particles.size(); ++i)
	    {
	        const auto & p = obj.particles[i];
		
		// Primaries
		if(!p.is_primary) continue;

		// Tracks (muons, pions, protons)
		if(p.pid < 2 || p.pid > 4) continue;
		
		// Get track end point (furthest from interaction vertex)
		TVector3 end_ptA(p.start_point[0], p.start_point[1], p.start_point[2]);
		TVector3 end_ptB(p.end_point[0], p.end_point[1], p.end_point[2]);
		TVector3 end_pt0;
		if((vertex - end_ptA).Mag() < (vertex - end_ptB).Mag())
		  end_pt0 = end_ptB;
		else
		  end_pt0 = end_ptA;

		// Second loop over nonprimaries
		for(size_t j(0); j < obj.particles.size(); ++j)
		{
		    const auto & q = obj.particles[j];
		    if(j == i) continue;

		    if(q.is_primary) continue;
		    if(q.pid < 2 || q.pid > 4) continue;

		    // Get track start point (closest to primary track end point)
		    TVector3 end_ptC(q.start_point[0], q.start_point[1], q.start_point[2]);
		    TVector3 end_ptD(q.end_point[0], q.end_point[1], q.end_point[2]);
		    TVector3 start_pt1;
		    if((end_pt0 - end_ptC).Mag() < (end_pt0 - end_ptD).Mag())
		      start_pt1 = end_ptC;
		    else
		      start_pt1 = end_ptD;
		    
		    // Check if first track end point is near second track start point
		    if((end_pt0 - start_pt1).Mag() < 1)
		      passes = false;

		} // end non primary track loop

	    } // end primary track loop

	    return passes;
        }

    template<class T>
        bool no_hadronic_michel_cut(const T & obj)
        {
	    bool passes(true);
	    for(size_t i(0); i < obj.particles.size(); ++i)
	    {
	        TVector3 vertex(obj.vertex[0], obj.vertex[1], obj.vertex[2]);
	      
	        const auto & p = obj.particles[i];
	        // primariy tracks
		if(!p.is_primary) continue;
		if(p.pid < 3 || p.pid > 4) continue;

		// Get hadron end point (furthest from interaction vertex)
                TVector3 end_ptA(p.start_point[0], p.start_point[1], p.start_point[2]);
                TVector3 end_ptB(p.end_point[0], p.end_point[1], p.end_point[2]);
                TVector3 end_pt0;
                if((vertex - end_ptA).Mag() < (vertex - end_ptB).Mag())
                  end_pt0 = end_ptB;
                else
                  end_pt0 = end_ptA;

		for(size_t j(0); j < obj.particles.size(); ++j)
	        {
		    if(j == 1) continue;
		    const auto & q = obj.particles[j];
		    if(q.is_primary) continue;
		    if(q.shape == 2)
		    {
		        // Is michel start point near hadronic endpoint
		        TVector3 michel_start(q.start_point[0], q.start_point[1], q.start_point[2]);
			if((michel_start - end_pt0).Mag() < 5)
			  passes = false;
		    }
		} // end nonprimary loop
	      
	    } // end primary loop

	    return passes;
        }


    template<class T>
        bool one_muon_cut(const T & obj)
        {
	    std::vector<uint32_t> c(utilities_ccpi0ana::count_primaries(obj));
	    return c[2] == 1;
        }

    template<class T>
      bool two_photons_cut(const T & obj)
      {
	std::vector<uint32_t> c(utilities_ccpi0ana::count_primaries(obj));
        return c[0] == 2;
      }

    template<class T>
        bool two_or_three_photons_cut(const T & obj)
        {
	    std::vector<uint32_t> c(utilities_ccpi0ana::count_primaries(obj));
	    return c[0] > 1 & c[0] < 4;
        }

      
    /**
     * @brief Apply pi0 mass cut.
     * @details This function applies a cut on the invariant diphoton mass
     * of the interactions two most energetic photons.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the pi0 mass cut.
     * @note This cut is intended to be used for the ccpi0ana analysis.
     */
    template<class T>
        bool valid_pi0_mass_cut(const T & obj)
        {
	  reco_inter s = utilities_ccpi0ana::reco_interaction_info(obj);
	  return s.pi0_mass < 400;
	}
    REGISTER_CUT_SCOPE(RegistrationScope::Both, valid_pi0_mass_cut, valid_pi0_mass_cut);

    /**
     * @brief Apply a fiducial volume, containment, flash time, 1mu0pi2gamma
     * topological, and pi0 mass cut (logical "and" of each).
     * @details This function applies a fiducial volume, containment, flash time,
     * 1mu0pi2gamma topological, and pi0 mass cut on the interaction using the logical "and"
     * of each previously defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * flash time, 1mu0pi2gamma topological, and pi0 mass cut.
     * @note This cut is intended to be used for the ccpi0ana analysis. 
     */
    template<class T>
    bool all_cut_icarus(const T & obj) {return fiducial_cut<T>(obj) && flash_cut<T>(obj, {-0.5, 1.6}) && base_topology_cut<T>(obj) && leading_shower_energy_cut<T>(obj) && valid_pi0_mass_cut<T>(obj);}
      //bool all_cut(const T & obj) {return fiducial_cut<T>(obj) && flash_cut<T>(obj) && base_topology_cut<T>(obj) && leading_shower_energy_cut<T>(obj);}

    /**
     * @brief Apply a fiducial volume, containment, flash time, 1mu0pi2gamma
     * topological, and pi0 mass cut (logical "and" of each).
     * @details This function applies a fiducial volume, containment, flash time,
     * 1mu0pi2gamma topological, and pi0 mass cut on the interaction using the logical "and"
     * of each previously defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * flash time, 1mu0pi2gamma topological, and pi0 mass cut.
     * @note This cut is intended to be used for the ccpi0ana analysis.
     */
    template<class T>
    bool all_cut_sbnd(const T & obj) {return fiducial_cut<T>(obj) && valid_flashmatch<T>(obj) && base_topology_cut<T>(obj) && leading_shower_energy_cut<T>(obj) && valid_pi0_mass_cut<T>(obj);}
    
    /**
     * @brief Apply a cut to select true CC 1mu interactions.
     * @details This function applies a cut on the final state
     * of the interaction.
     * @param obj the interaction to select on.
     * @return true if the interaction meets the CC 1mu criteria.
     */
    template<class T>
    bool cc1mu(const caf::SRInteractionTruthDLPProxy & obj)
    {
        truth_inter s = utilities_ccpi0ana::truth_interaction_info(obj);
	return s.is_cc && s.num_primary_muons_thresh == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::True, cc1mu, cc1mu);

    /**
     * @brief Apply a cut to select true 0pi interactions.
     * @details This function applies a cut on the final state
     * of the interaction.
     * @param obj the interaction to select on.
     * @return true if the interaction meets the 0pi criteria. 
     */
    template<class T>
    bool zero_pi(const caf::SRInteractionTruthDLPProxy & obj)
    {
        truth_inter s = utilities_ccpi0ana::truth_interaction_info(obj);
	return s.num_primary_pions_thresh == 0;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::True, zero_pi, zero_pi);

    /**
     * @brief Apply a cut to select true 1pi0 interactions.
     * @details This function applies a cut on the final state
     * of the interaction.
     * @param obj the interaction to select on.
     * @return true if the interaction meets the 1pi0 criteria.
     */
    template<class T>
    bool one_pi0(const caf::SRInteractionTruthDLPProxy & obj)
    {
        truth_inter s = utilities_ccpi0ana::truth_interaction_info(obj);
        return s.num_primary_pi0s_thresh == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::True, one_pi0, one_pi0);
    
    /**
     * @brief Apply a cut to select the 1mu0pi1pi0 signal.
     * @details This function applies a cut on the final state
     * of the interaction.
     * @param obj the interaction to select on.
     * @return true if the interaction meets the signal definition criteria.
     * @note This cut is intended to be used for the ccpi0ana analysis.
     */
    template<class T>
    bool signal_1mu0pi1pi0(const caf::SRInteractionTruthDLPProxy & obj)
    {
        truth_inter s = utilities_ccpi0ana::truth_interaction_info(obj);
        return s.num_primary_muons_thresh == 1 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.is_cc && s.is_neutrino;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::True, signal_1mu0pi1pi0, signal_1mu0pi1pi0);

    /**
     * @brief Apply a cut to select other (non 1mu0pi1pi0) neutrino interactions.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction.  This is the "true" non-signal cut.
     * (neutrino inteaction, but not 1mu0pi1pi0)
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment, 
     * and non-1mu0pi1pi0 topological cut.
     * @note This cut is intended to be used for the ccpi0ana analysis.
     */
    bool other_nu_1mu0pi1pi0(const caf::SRInteractionTruthDLPProxy & obj)
        {
	  truth_inter s = utilities_ccpi0ana::truth_interaction_info(obj);
	  return !(s.num_primary_muons_thresh == 1 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.is_cc) && s.is_neutrino;
        }

}
#endif // CUTS_CCPI0ANA_H
