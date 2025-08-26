/**
 * @file utilities_pi0ana.h
 * @brief Header file for definitions of utility functions for supporting
 * analysis variables and cuts.
 * @details This file contains definitions of utility functions which are used
 * to support the implementation of analysis variables and cuts. These functions
 * are intended to be used to simplify the implementation of variables and cuts
 * by providing common functionality which can be reused across multiple
 * variables and cuts.
 * @author lkashur@colostate.edu
 */
#ifndef UTILITIES_PI0ANA_H
#define UTILITIES_PI0ANA_H

#define PI0_MASS 134.9768
 
#include <iostream>
#include <vector>
#include <TVector3.h>
#include <unordered_map>
#include "include/framework.h"

struct pi0 {
  size_t leading_shower_index;
  double leading_shower_ke;
  double leading_shower_conv_dist;
  size_t subleading_shower_index;
  double subleading_shower_ke;
  double subleading_shower_conv_dist;
  double showers_costheta;
  double mass;
  double momentum;
  double beam_costheta;
};

/**
 * @namespace utilities_pi0ana
 * @brief Namespace for organizing utility functions for supporting analysis
 * variables and cuts.
 * @details This namespace is intended to be used for organizing utility
 * functions which are used to support the implementation of analysis variables
 * and cuts. These functions are intended to be used to simplify the
 * implementation of variables and cuts by providing common functionality which
 * can be reused across multiple variables and cuts.
 * @note The namespace is intended to be used in conjunction with the
 * vars and cuts namespaces, which are used for organizing variables and cuts
 * which act on interactions.
 */
namespace utilities_pi0ana
{
    static std::vector<double> reco_pi0_shower_thresholds = {20.0, 40.0};

    /**
     * @brief Group true pi0 daughters by parent ID
     * @param obj the interaction from which to count
     * @return an unordered map with parent ID as key and daughter particle indices as value
     */
    std::unordered_map< int, std::vector<size_t> > get_true_pi0s(const caf::SRInteractionTruthDLPProxy & obj, bool primaries = true, std::vector<double> params = {0.0,})
    {
        std::unordered_map< int, std::vector<size_t> > true_pi0s;
	
	/**
	 * Loop over particles, grouping photons
	 * and electrons that share a parent pi0.
	 */
	for(size_t i(0); i < obj.particles.size(); ++i)
	{
	    const auto & p = obj.particles[i];
	    

	    if(primaries == true and !p.is_primary) continue;
	    if(primaries == false and p.is_primary) continue;

	    if(p.parent_pdg_code == 111 && (p.pdg_code == 22 || p.pdg_code == 11 || p.pdg_code == -11))
	    {
	        true_pi0s[p.parent_track_id].push_back(i);
	    }
	}

	/**
	 * Remove pi0 candidates that are either
	 * subthreshold or have fewer than two daughters.
	 */
	std::vector<size_t> bad_pi0_ids;
        for(auto const & _pi0 : true_pi0s)
        {

	  int num_pi0_daughters(0);
	  TVector3 pi0_momentum(0,0,0);
	  double pi0_ke(0);
	  for(auto & _pidx : _pi0.second)
	    {
	      const auto & _p = obj.particles[_pidx];
	      TVector3 _p_momentum(_p.momentum[0], _p.momentum[1], _p.momentum[2]);
	      pi0_momentum += _p_momentum;
	      num_pi0_daughters++;
	    }
	    pi0_ke = std::sqrt(std::pow(PI0_MASS, 2) + std::pow(pi0_momentum.Mag(), 2)) - PI0_MASS;
            if(num_pi0_daughters < 2 || pi0_ke < params[0]) bad_pi0_ids.push_back(_pi0.first);
	}
	for(size_t i=0; i<bad_pi0_ids.size(); i++)
	{
            true_pi0s.erase(bad_pi0_ids[i]);
	}
	return true_pi0s;
    }

    /**
     * @brief Count the number of true primary neutral pions in interaction.
     * @param obj the interaction from which to count.
     * @param params the neutral pion kinetic energy threshold to consider.
     * @return the number of true primary neutral pions in interaction.
     */
    double true_primary_pi0_multiplicity(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params={0.0,})
    {
        std::unordered_map< int, std::vector<size_t> > true_primary_pi0s = get_true_pi0s(obj, true, params);
	return true_primary_pi0s.size();
    }

    /**
     * @brief Count the number of true nonprimary neutral pions in interaction.
     * @param obj the interaction from which to count.
     * @param params the neutral pion kinetic energy threshold to consider.
     * @return the number of true nonprimary neutral pions in interaction.
     */
    double true_nonprimary_pi0_multiplicity(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params={0.0,})
    {
        std::unordered_map< int, std::vector<size_t> > true_nonprimary_pi0s = get_true_pi0s(obj, false, params);
	return true_nonprimary_pi0s.size();
    }

    /**
     * @brief Extract information about neutral pion.
     * @details Extract information about neutral pion from its daughter particles, 
     * assuming only one neutral pion is present in interaction.
     * @param obj the interaction from which to count.
     * @param params the neutral pion kinetic energy threshold to consider.
     * @return the number of true nonprimary neutral pions in interaction.
     */
    template<class T>
    pi0 pi0_info(const T & obj)
    {
        /**
	 * Initialize output structure.
	 */
        pi0 s;
	TVector3 momentum(0,0,0);
	double beam_costheta(-5);
	size_t leading_shower_index(kNoMatch);
	double leading_shower_ke(-5);
	double leading_shower_conv_dist(-5);
        size_t subleading_shower_index(kNoMatch);
        double subleading_shower_ke(-5);
        double subleading_shower_conv_dist(-5);
	double showers_costheta(-5);
	double mass(-5);

	/**
	 * Neutrino beam direction and neutrino interaction vertex,
	 * to be used in later calculations.
	 */
	TVector3 beamdir(0,0,1);
	TVector3 vertex(obj.vertex[0], obj.vertex[1], obj.vertex[2]);

	/**
	 * True neutral pions
	 */
	if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
		       {
			   /**
			    * Retrieve map with true neutral pion daughter indices.
			    */
			   std::unordered_map< int, std::vector<size_t> > true_primary_pi0s = get_true_pi0s(obj, true);
			   
			   /**
			    * Calculate neutral pion momentum from daugher particles.
			    * Ideally the true neutral pion momentum would be stored
			    * in SPINE, but this is not currently the case.
			    */
			   int num_pi0_daughters(0);
			   int num_pi0_photon_daughters(0);
			   std::vector<size_t> pi0_daughter_indices;
			   for(auto const & _pi0 : true_primary_pi0s)
			   {
			       for(auto & _pidx : _pi0.second)
			       {
				   const auto & _p = obj.particles[_pidx];
				   TVector3 _p_momentum(_p.momentum[0], _p.momentum[1], _p.momentum[2]);
				   momentum += _p_momentum;
				   pi0_daughter_indices.push_back(_pidx);
				   num_pi0_daughters++;
				   if(_p.pid == 0) num_pi0_photon_daughters++;
			       }
			   }
			   
			   /**
			    * Calculate neutral pion angle with respect to neutrino beam.
			    */
			   beam_costheta = momentum.Unit().Dot(beamdir);
			   
			   /**
			    * If neutral pion decays to two photons (the usual case),
			    * retrieve information about leading and subleading photons.
			    */
			   if(num_pi0_photon_daughters == 2)
			   {
			       const auto & pi0_daughter0 = obj.particles[pi0_daughter_indices[0]];
			       const auto & pi0_daughter1 = obj.particles[pi0_daughter_indices[1]];

			       size_t _leading_shower_index;
			       size_t _subleading_shower_index;
			       if(pi0_daughter0.ke > pi0_daughter1.ke)
			       {
				   _leading_shower_index = pi0_daughter_indices[0];
				   _subleading_shower_index = pi0_daughter_indices[1];
			       }
			       else
			       {
				   _leading_shower_index = pi0_daughter_indices[1];
				   _subleading_shower_index = pi0_daughter_indices[0];
			       }
			       
			       /**
				* Leading photon information
				*/
			       const auto & leading_shower = obj.particles[_leading_shower_index];
			       leading_shower_index = _leading_shower_index;
			       leading_shower_ke = leading_shower.ke;
			       TVector3 leading_shower_start_point(leading_shower.start_point[0], leading_shower.start_point[1], leading_shower.start_point[2]);
			       TVector3 leading_shower_dir(leading_shower.momentum[0], leading_shower.momentum[1], leading_shower.momentum[2]);
			       leading_shower_dir = leading_shower_dir.Unit();
			       leading_shower_conv_dist = (vertex - leading_shower_start_point).Mag();
			       TVector3 leading_shower_momentum(leading_shower.momentum[0], leading_shower.momentum[1], leading_shower.momentum[2]);
			       
			       /**
				* Subleading photon information
				*/
			       const auto & subleading_shower = obj.particles[_subleading_shower_index];
                               subleading_shower_index = _subleading_shower_index;
                               subleading_shower_ke = subleading_shower.ke;
                               TVector3 subleading_shower_start_point(subleading_shower.start_point[0], subleading_shower.start_point[1], subleading_shower.start_point[2]);
                               TVector3 subleading_shower_dir(subleading_shower.momentum[0], subleading_shower.momentum[1], subleading_shower.momentum[2]);
                               subleading_shower_dir = subleading_shower_dir.Unit();
                               subleading_shower_conv_dist = (vertex - subleading_shower_start_point).Mag();
                               TVector3 subleading_shower_momentum(subleading_shower.momentum[0], subleading_shower.momentum[1], subleading_shower.momentum[2]);

			       /**
				* Information derived from leading and subleading photons.
				*/
			       showers_costheta = leading_shower_dir.Dot(subleading_shower_dir);
			       mass = sqrt(2*leading_shower_ke*subleading_shower_ke*(1-showers_costheta));
			   }
			   
			   /**
			    * Store output.
			    */
			   s.momentum = momentum.Mag()/1000.0;
			   s.beam_costheta = beam_costheta;
			   s.leading_shower_index = leading_shower_index;
			   s.leading_shower_ke = leading_shower_ke/1000.0;
			   s.leading_shower_conv_dist = leading_shower_conv_dist;
			   s.subleading_shower_index = subleading_shower_index;
                           s.subleading_shower_ke = subleading_shower_ke/1000.0;
                           s.subleading_shower_conv_dist = subleading_shower_conv_dist;
			   s.showers_costheta = showers_costheta;
			   s.mass = mass;
		       }
	
	/**
	 * Reco neutral pions
	 */
	else
	{
	    /**
	     * Loop over every possible pair of photons.
	     * Pair that has best agreement with true neutral pion mass
	     * is chosen as belonging to candidate neutral pion.
	     */
	    std::vector<std::pair< std::pair<size_t, size_t>, double> > showers_metric;
	    for(size_t i(0); i < obj.particles.size(); ++i)
	    {
	        /**
		 * Photon 0
		 */
	        const auto & p = obj.particles[i];
		if(!(p.is_primary && p.pid == 0)) continue;
		
		TVector3 sh0_start(p.start_point[0], p.start_point[1], p.start_point[2]);
		TVector3 sh0_start_dir = (sh0_start - vertex).Unit();
		TVector3 sh0_dir(p.start_dir[0], p.start_dir[1], p.start_dir[2]);

		for(size_t j(0); j < obj.particles.size(); ++j)
		{
		  if(j == i) continue;
		  
		  /**
		   * Photon 1
		   */
		  const auto & q = obj.particles[j];
		  if(!(q.is_primary && q.pid == 0)) continue;

		  TVector3 sh1_start(q.start_point[0], q.start_point[1], q.start_point[2]);
		  TVector3 sh1_start_dir = (sh1_start - vertex).Unit();
		  TVector3 sh1_dir(q.start_dir[0], q.start_dir[1], q.start_dir[2]);

		  double _leading_shower_ke(0);
		  double _subleading_shower_ke(0);

		  if(pvars::calo_ke(p) > pvars::calo_ke(q))
		    {
		      _leading_shower_ke = pvars::calo_ke(p);
		      _subleading_shower_ke = pvars::calo_ke(q);
		    }
		  else
		    {
		      _leading_shower_ke = pvars::calo_ke(q);
		      _subleading_shower_ke = pvars::calo_ke(p);
		    }
		  if(_leading_shower_ke < reco_pi0_shower_thresholds[1] || _subleading_shower_ke < reco_pi0_shower_thresholds[0]) continue;

		  /**
		   * Calculate invariant mass of photon pair.
		   */
		  double _showers_costheta = sh0_start_dir.Dot(sh1_start_dir);
		  double _mass = sqrt(2*_leading_shower_ke*_subleading_shower_ke*(1-_showers_costheta));
		  showers_metric.push_back(std::make_pair(std::make_pair(i, j), mass));
		  //std::sort(showers_metric.begin(), showers_metric.end(), [](const std::pair<std::pair<size_t, size_t>, double> &a, const std::pair<std::pair<size_t, size_t>, double> &b)
		  //{ return abs(a.second - 134.9768) < abs(b.second - 134.9768);});
		}
	    }
	    
	    /**
	     * Sort combinations of photon pairs by agreement with true neutral pion mass
	     * and identify leading and subleading photons.
	     */
	    if(!showers_metric.empty())
	    {
	        std::sort(showers_metric.begin(), showers_metric.end(), [](const std::pair<std::pair<size_t, size_t>, double> &a, const std::pair<std::pair<size_t, size_t>, double> &b)
			  { return abs(a.second - PI0_MASS) < abs(b.second - PI0_MASS);});
	      
	        std::pair<size_t, size_t> sh_pair_ids = showers_metric[0].first; // best metric agreement
		if(pvars::calo_ke(obj.particles[sh_pair_ids.first]) > pvars::calo_ke(obj.particles[sh_pair_ids.second]))
		{
		    leading_shower_index = sh_pair_ids.first;
		    subleading_shower_index = sh_pair_ids.second;
		}
		else
		{
		    leading_shower_index = sh_pair_ids.second;
		    subleading_shower_index = sh_pair_ids.first;
		}

		/**
		 * Leading photon information
		 */
		const auto & leading_shower = obj.particles[leading_shower_index];
		leading_shower_ke = pvars::calo_ke(leading_shower);
		TVector3 leading_shower_start_point(leading_shower.start_point[0], leading_shower.start_point[1], leading_shower.start_point[2]);
		TVector3 leading_shower_dir = leading_shower_start_point - vertex;
		leading_shower_dir = leading_shower_dir.Unit();
		leading_shower_conv_dist = (vertex - leading_shower_start_point).Mag();
		TVector3 leading_shower_momentum = leading_shower_ke * leading_shower_dir;

		/**
		 * Subleading photon information
		 */
		const auto & subleading_shower = obj.particles[subleading_shower_index];
                subleading_shower_ke = pvars::calo_ke(subleading_shower);
                TVector3 subleading_shower_start_point(subleading_shower.start_point[0], subleading_shower.start_point[1], subleading_shower.start_point[2]);
                TVector3 subleading_shower_dir = subleading_shower_start_point - vertex;
                subleading_shower_dir = subleading_shower_dir.Unit();
                subleading_shower_conv_dist = (vertex - subleading_shower_start_point).Mag();
		TVector3 subleading_shower_momentum = subleading_shower_ke * subleading_shower_dir;

		/**
		 * Information derived from leading and subleading photons
		 */
		showers_costheta = leading_shower_dir.Dot(subleading_shower_dir);
		mass = sqrt(2*leading_shower_ke*subleading_shower_ke*(1-showers_costheta));
		momentum = leading_shower_momentum + subleading_shower_momentum;
		beam_costheta = momentum.Unit().Dot(beamdir);

		/**
		 * Store output.
		 */
		s.momentum = momentum.Mag()/1000.0;
		s.beam_costheta = beam_costheta;
		s.leading_shower_index = leading_shower_index;
		s.leading_shower_ke = leading_shower_ke/1000.0;
		s.leading_shower_conv_dist = leading_shower_conv_dist;
		s.subleading_shower_index = subleading_shower_index;
		s.subleading_shower_ke = subleading_shower_ke/1000.0;
		s.subleading_shower_conv_dist = subleading_shower_conv_dist;
		s.showers_costheta = showers_costheta;
		s.mass = mass;
	    }
	    else
	    {
	        s.momentum = -5;
		s.beam_costheta = -5;
		s.leading_shower_index = kNoMatch;
		s.leading_shower_ke = -5;
		s.leading_shower_conv_dist = -5;
		s.subleading_shower_index = kNoMatch;
		s.subleading_shower_ke = -5;
		s.subleading_shower_conv_dist = -5;
		s.showers_costheta = -5;
		s.mass = -5;
	    }
	}

	return s;
    }
    
}
#endif // UTILITIES_PI0ANA_H
