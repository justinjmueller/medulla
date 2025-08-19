/**
 * @file vars_ncpi0ana.h
 * @brief Header file for definitions of analysis variables specific to the
 * ncpi0ana analysis.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from interactions specific to the ncpi0ana
 * analysis. Each variable is implemented as a function which takes an
 * interaction object as an argument and returns a double. These are the
 * building blocks for producing high-level plots of the selected interactions.
 * @author lkashur@colostate.edu
 */
#ifndef VARS_NCPI0ANA_H
#define VARS_NCPI0ANA_H

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"

#include <iostream>

#include "include/selectors.h"
#include "include/framework.h"

#include "include/cuts.h"
#include "include/pi0ana/cuts_ncpi0ana.h"

/**
 * @namespace vars::ncpi0ana
 * @brief Namespace for organizing variables specific to the ncpi0ana analysis.
 * @details This namespace is intended to be used for organizing variables which
 * act on interactions specific to the ncpi0ana analysis. Each variable is
 * implemented as a function which takes an interaction object as an argument
 * and returns a double. The function should be templated on the type of
 * interaction object if the variable is intended to be used on both true and
 * reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the vars
 * namespace, which is used for organizing generic variables which act on
 * interactions.
 */
namespace vars::ncpi0ana
{
    /**
     * @brief GUNDAM variable for enumerating interaction categories.        
     * @details This variable provides a basic categorization of interactions
     * using only signal, neutrino background, and cosmic background as the
     * three categories.
     * 1: Signal
     * 2: Signal (OOPS)
     * 3: Other nu
     * 4: Cosmic
     * @param obj the interaction to apply the variable on.
     * @return the enumerated category of the interaction. 
     */
    double is_signal_mc(const caf::SRInteractionTruthDLPProxy & obj)
    {
      truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);

      // Cosmic                                                                                                              
      uint16_t cat(4);

      // Nu                                                                                                                         
      if(s.is_neutrino)
        {
	  // Signal
	  if(s.num_primary_muons_thresh == 1 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.is_cc && s.is_fiducial) cat = 1;

	  // Signal (OOPS)
	  else if( (s.num_primary_muons == 1 && s.num_primary_pions == 0 && s.num_primary_pi0s == 1 && s.is_cc && s.is_fiducial) && (s.num_primary_muons_thresh != 1 || s.num_primary_pions_thresh != 0 || s.num_primary_pi0s_thresh != 1) ) cat = 2;

	  // Other nu
	  else cat = 3;
        }

      return cat;
    }


    /**
     * @brief GUNDAM variable for enumerating interaction categories.
     * @details This variable provides a basic categorization of interactions
     * using only signal, neutrino background, and cosmic background as the
     * three categories.
     * @param obj the interaction to apply the variable on.
     * @return the enumerated category of the interaction. 
     */
    double is_signal_data(const caf::SRInteractionTruthDLPProxy & obj)
    {
      int cat(-5);
      return cat;
    }


    /**
     * @brief Variable for enumerating interaction categories (simple).
     * @details This variable provides a basic categorization of interactions
     * using only signal pi0s, other pi0s, other nus, and cosmics.
     * @param obj the interaction to apply the variable on.
     * @return the enumerated category of the interaction.
     */
    /*
    template<class T>
    double category(const caf::SRInteractionTruthDLPProxy & obj)
    {
        truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);
	
	// Cosmic
	uint16_t cat(3);
	
	// Neutrino
	if(s.is_neutrino)
	{
	    // 1mu0pi1pi0 (in-, fiducial) 
	    if(s.num_primary_muons_thresh == 1 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.is_cc && s.is_fiducial) cat = 0;
	    // Other nu-induced pi0
	    else if(s.num_primary_pi0s >= 1) cat = 1;
	    // Other nu without pi0
	    else if(s.num_primary_pi0s == 0) cat = 2;
	}
	return cat;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, category, category);
    */
    
    /**
     * @brief Variable for enumerating interaction categories.
     * @details This variable provides a basic categorization of interactions
     * using the following categories:
     * 0: 1mu0pi1pi0 (in-, fiducial)
     * 1: 1mu0pi1pi0 (OOPS, fiducial)
     * 2: 1mu0pi1pi0 (OOFV)
     * 3: 1muNpi1pi0
     * 4: 1muNpi0pi0
     * 5: 1muNpi0
     * 6: NC 1pi0
     * 7: Other nu
     * 8: Cosmic
     * @param obj the interaction to apply the variable on.
     * @return the enumerated category of the interaction.
     */
    template<class T>
      double category_topology_v1(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params={})
      {
	truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);

	// Cosmic
	uint16_t cat(8);

	// Neutrino
	if(s.is_neutrino)
	  {

	    /////////////////////////////
	    /// Exactly params[1] protons
	    /////////////////////////////
	    if(params[0] == -2)
	      {
		// 0mu0pi1pi0 + specified proton count (in-phase, fiducial)
		if(s.num_primary_muons_thresh == 0 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.num_primary_protons_thresh == params[1] && !s.is_cc && s.is_fiducial) cat = 0;

		// 0mu0pi1pi0 + specified proton count (OOPS, fiducial)
		else if( (s.num_primary_muons == 0 && s.num_primary_pions == 0 && s.num_primary_pi0s == 1 && s.num_primary_protons == params[1] && !s.is_cc && s.is_fiducial) && (s.num_primary_muons_thresh != 0 || s.num_primary_pions_thresh != 0 || s.num_primary_protons_thresh != params[1] || s.num_primary_pi0s_thresh != 1) ) cat = 1;

		// 0mu0pi1pi0 + specified proton count (OOFV)
		else if( (s.num_primary_muons_thresh == 0 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.num_primary_protons_thresh == params[1] && !s.is_cc && !s.is_fiducial) ) cat = 2;
	 	
		// 0muNpi1pi0
		else if(s.num_primary_muons_thresh == 0 && s.num_primary_pions_thresh > 0 && s.num_primary_pi0s_thresh == 1 && !s.is_cc && s.is_fiducial) cat = 3;

		// 0muNpi0pi0
		else if(s.num_primary_muons_thresh == 0 && s.num_primary_pions_thresh > 0 && s.num_primary_pi0s_thresh == 0 && !s.is_cc && s.is_fiducial) cat = 4;

		// 0muNpi0
		else if(s.num_primary_muons_thresh == 0 && s.num_primary_pi0s_thresh > 1 && !s.is_cc && s.is_fiducial) cat = 5;
		
		// CC1pi0
		else if(s.num_primary_muons_thresh == 1 && s.num_primary_pi0s_thresh == 1 && s.is_cc && s.is_fiducial) cat = 6;
		
		// Other nu
		else cat = 7;	
	      }
	    //////////////////////////////
	    /// At least params[1] protons
	    //////////////////////////////
	    else if(params[0] == -1)
	      {
		// 0mu0pi1pi0 + specified proton count (in-phase, fiducial)
		if(s.num_primary_muons_thresh == 0 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.num_primary_protons_thresh >= params[1] && !s.is_cc && s.is_fiducial) cat = 0;
		
		// 0mu0pi1pi0 + specified proton count (OOPS, fiducial)
		else if( (s.num_primary_muons == 0 && s.num_primary_pions == 0 && s.num_primary_pi0s == 1 && s.num_primary_protons >= params[1] && !s.is_cc && s.is_fiducial) && (s.num_primary_muons_thresh != 0 || s.num_primary_pions_thresh != 0 || s.num_primary_protons_thresh < params[1] || s.num_primary_pi0s_thresh != 1) ) cat = 1;

		// 0mu0pi1pi0 + specified proton count (OOFV)
		else if( (s.num_primary_muons_thresh == 0 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.num_primary_protons_thresh >= params[1] && !s.is_cc && !s.is_fiducial) ) cat = 2;

		// 0muNpi1pi0
                else if(s.num_primary_muons_thresh == 0 && s.num_primary_pions_thresh > 0 && s.num_primary_pi0s_thresh == 1 && !s.is_cc && s.is_fiducial) cat = 3;

                // 0muNpi0pi0
                else if(s.num_primary_muons_thresh == 0 && s.num_primary_pions_thresh > 0 && s.num_primary_pi0s_thresh == 0 && !s.is_cc && s.is_fiducial) cat = 4;

                // 0muNpi0
                else if(s.num_primary_muons_thresh == 0 && s.num_primary_pi0s_thresh > 1 && !s.is_cc && s.is_fiducial) cat = 5;

                // CC1pi0
                else if(s.num_primary_muons_thresh == 1 && s.num_primary_pi0s_thresh == 1 && s.is_cc && s.is_fiducial) cat = 6;

                // Other nu
                else cat = 7;
	      }
	    else cat = -5;	    
	  }
	return cat;
      }
    REGISTER_VAR_SCOPE(RegistrationScope::True, category_topology_v1, category_topology_v1);

    
    template<class T>
    double category_topology_v2(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params={})
    {
        truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);
	
	// Cosmic
	uint16_t cat(10);
	
	// Neutrino
	if(s.is_neutrino)
	{
	  /////////////////////////////
	  /// Exactly params[1] protons
	  /////////////////////////////
	  if(params[0] == -2)
	  {
	    // 0mu0pi1pi0 (in-phase, fiducial)
	    if(s.num_primary_muons_thresh == 0 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.num_primary_protons_thresh == params[1] && !s.is_cc && s.is_fiducial) cat = 0;
	    // Other nu-induced pi0
	    else if(s.num_primary_pi0s >= 1) cat = 1;

	    // Other nu without pi0
	    else if(s.num_primary_pi0s == 0) cat = 2;
	  }
	  //////////////////////////////
	  /// At least params[1] protons
	  //////////////////////////////
	  else if(params[0] == -1)
	  {
	    // 0mu0pi1pi0 (in-phase, fiducial) 
            if(s.num_primary_muons_thresh == 0 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && s.num_primary_protons_thresh >= params[1] && !s.is_cc && s.is_fiducial) cat = 0;
            // Other nu-induced pi0
            else if(s.num_primary_pi0s >= 1) cat = 1;

            // Other nu without pi0
            else if(s.num_primary_pi0s == 0) cat = 2;
	  }

	}
	return cat;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, category_topology_v2, category_topology_v2);

    
    template<class T>
    double category_topology_v3(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params={})
    {
        truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);

        // Cosmic
	uint16_t cat(10);

	// Neutrino
	if(s.is_neutrino)
	{
	  
	  // 0mu 0pi 1pi0 (in-phase, fiducial)
	  if(s.num_primary_muons_thresh == 0 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && !s.is_cc && s.is_fiducial) cat = 0;

	  // 0mu 0pi 1pi0 (OOPS, fiducial)
	  else if( (s.num_primary_muons == 0 && s.num_primary_pions == 0 && s.num_primary_pi0s == 1 && !s.is_cc && s.is_fiducial) && (s.num_primary_muons_thresh != 0 || s.num_primary_pions_thresh != 0 || s.num_primary_pi0s_thresh != 1) ) cat = 1;

	  // 0mu 0pi 1pi0 (OOFV)
	  else if(s.num_primary_muons_thresh == 0 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 1 && !s.is_cc && !s.is_fiducial) cat = 2;

	  // 0mu 0pi (2+ pi0)
	  else if(s.num_primary_muons_thresh == 0 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh >= 2 && !s.is_cc && s.is_fiducial) cat = 3;

	  // 0mu 0pi 0pi0 Npi0_nonprim
	  else if(s.num_primary_muons_thresh == 0 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 0 && s.num_nonprimary_pi0s >= 1 && !s.is_cc && s.is_fiducial) cat = 4;

	  // 0mu 0pi 0pi0 0pi0_nonprim
	  else if(s.num_primary_muons_thresh == 0 && s.num_primary_pions_thresh == 0 && s.num_primary_pi0s_thresh == 0 && s.num_nonprimary_pi0s == 0 && !s.is_cc && s.is_fiducial) cat = 5;

	  // 0mu Npi Xpi0
	  else if(s.num_primary_muons_thresh == 0 && s.num_primary_pions_thresh >= 1 && !s.is_cc && s.is_fiducial) cat = 6;

	  // 1mu Npi0                                                                                                                                                                                    
	  else if(s.num_primary_muons_thresh == 1 && s.num_primary_pi0s_thresh >= 1 && s.is_cc && s.is_fiducial) cat = 7;

	  // Other nu
	  else cat = 8;	  
	  
	}
	
	return cat;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, category_topology_v3, category_topology_v3);

    /**
     * @brief Dummy GUNDAM variables.
     * @details "cut_type" specifies a signal or sideband cut.
     * "is_data" specifies data or MC.  "is_nu" specifies neutrino or cosmic.
     * @param obj the interaction to apply the variable on.
     * @return the dummy GUNDAM variable. 
     */
    template<class T>
        double cut_type(const T & obj)
        {
	    // Signal
	    double cat(1);

	    // Sideband
	    //cat = 2;
	  
	    return cat;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, cut_type, cut_type);

    template<class T>
        double is_data(const T & obj)
	{
            double cat(-5);
            return cat;
	}
    REGISTER_VAR_SCOPE(RegistrationScope::Both, is_data, is_data);

    template<class T>
        double is_nu(const T & obj)
        {
	    double cat(-5);
	    return cat;
	}
    REGISTER_VAR_SCOPE(RegistrationScope::Both, is_nu, is_nu);

    /**
     * @brief Variable for the base topology status of the interaction.
     * @details This variable holds the status of whether or not the
     * the interaction passes the base topology cut defined in cuts_ncpi0ana. 
     */
    template<class T>
    double base_topology_satisfied(const T & obj) {return cuts::ncpi0ana::base_topology_cut(obj);}
    REGISTER_VAR_SCOPE(RegistrationScope::Reco, base_topology_satisfied, base_topology_satisfied);

    template<class T>
    double num_protons_satisfied(const T & obj, std::vector<double> params={}) {return cuts::ncpi0ana::num_protons_cut(obj, params);}
    REGISTER_VAR_SCOPE(RegistrationScope::Reco, num_protons_satisfied, num_protons_satisfied);  

    /**
     * @brief Variable for the leading shower energy threshold status of the interaction.
     * @details This variable holds the status of whether or not the
     * the interaction passes the leading shower energy cut defined in cuts_ncpi0ana.
     */
    template<class T>
    double leading_shower_energy_satisfied(const T & obj) {return cuts::ncpi0ana::leading_shower_energy_cut(obj);}
    REGISTER_VAR_SCOPE(RegistrationScope::Reco, leading_shower_energy_satisfied, leading_shower_energy_satisfied);

    /**
     * @brief Variable for the pi0 mass cut status of the interaction.
     * @details This variable holds the status of whether or not the
     * the interaction passes the pi0 mass cut defined in cuts_ncpi0ana.
     */
    template<class T>
    double valid_pi0_mass_satisfied(const T & obj) {return cuts::ncpi0ana::valid_pi0_mass_cut(obj);}
    REGISTER_VAR_SCOPE(RegistrationScope::Reco, valid_pi0_mass_satisfied, valid_pi0_mass_satisfied);

    /**
     * @brief Variable for the status of all interaction cuts.
     * @details This variable holds the status of whether or not the
     * the interaction passes the "all" cut defined in cuts_ncpi0ana.
     */
    template<class T>
    double all_cut_icarus_satisfied(const T & obj) {return cuts::ncpi0ana::all_cut_icarus(obj);}
    REGISTER_VAR_SCOPE(RegistrationScope::Reco, all_cut_icarus_satisfied, all_cut_icarus_satisfied);

    template<class T>
    double all_cut_sbnd_satisfied(const T & obj, std::vector<double> params={}) {return cuts::ncpi0ana::all_cut_sbnd(obj, params);}
    REGISTER_VAR_SCOPE(RegistrationScope::Reco, all_cut_sbnd_satisfied, all_cut_sbnd_satisfied);

    /**
     * @brief Variable for pi0 leading photon energy.
     * @details Variable for pi0 leading photon energy
     * [MeV], as calculated with p.calo_ke attribute.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the pi0 leading photon energy.
     */
    template<class T>
        double pi0_leading_photon_energy(const T & obj)
        {
	    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
			   {
			       truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);
			       return s.pi0_leading_photon_energy;
			   }
            else
	    {
		reco_inter s = utilities_ncpi0ana::reco_interaction_info(obj);
		return s.pi0_leading_photon_energy;
	    }
	}
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pi0_leading_photon_energy, pi0_leading_photon_energy);

    template<class T>
        double pi0_leading_photon_start_dedx(const T & obj)
        {
	    reco_inter s = utilities_ncpi0ana::reco_interaction_info(obj);
	    return s.pi0_leading_photon_start_dedx;
	}

    /**
     * @brief Variable for pi0 leading photon conversion distance.
     * @details Variable for pi0 leading photon conversion distance
     * [cm], as calculated using interaction vertex and shower start point.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the pi0 leading photon conversion distance.
     */
    template<class T>
      double pi0_leading_photon_conv_dist(const T & obj)
      {
	  if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
			 {
			     truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);
			     return s.pi0_leading_photon_conv_dist;
			 }
	  else
	  {
	      reco_inter s = utilities_ncpi0ana::reco_interaction_info(obj);
	      return s.pi0_leading_photon_conv_dist;
	  }
      }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pi0_leading_photon_conv_dist, pi0_leading_photon_conv_dist);

    /**
     * @brief Variable for angle between pi0 leading photon cluster direction and vertex direction.
     * @details Variable for angle between pi0 leading photon cluster direction, calculated as the normalized
     * mean direction of the cluster, and the the pi0 leading photon vertex direction, calculated from the
     * vector between the interaction vertex and shower start point.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the angle between the pi0 leading photon cluster direction and vertex direction. 
     */
    template<class T>
      double pi0_leading_photon_cosphi(const T & obj)
      {
	if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
		       {
			 truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);
			 return -5;
			 //return s.pi0_leading_photon_cosphi;
		       }
	else
	  {
	    reco_inter s = utilities_ncpi0ana::reco_interaction_info(obj);
	    return s.pi0_leading_photon_cosphi;
	  }
      }

    template<class T>
        double pi0_leading_photon_ip(const T & obj)
        {
	    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
			   {
			       truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);
			       return -5;
			   }
	    else
	    {
	        reco_inter s = utilities_ncpi0ana::reco_interaction_info(obj);
		return s.pi0_leading_photon_ip;
	    }
	}

    /**
     * @brief Variable for pi0 subleading photon energy.
     * @details Variable for pi0 subleading photon energy
     * [MeV], as calculated with p.calo_ke attribute.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the pi0 subleading photon energy. 
     */
    template<class T>
        double pi0_subleading_photon_energy(const T & obj)
	{
	    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                             {
			       truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);
			       return s.pi0_subleading_photon_energy;
                             }
            else
	    {
	        reco_inter s = utilities_ncpi0ana::reco_interaction_info(obj);
                return s.pi0_subleading_photon_energy;
	    }
	}
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pi0_subleading_photon_energy, pi0_subleading_photon_energy);

    template<class T>
        double pi0_subleading_photon_start_dedx(const T & obj)
        {
	    reco_inter s = utilities_ncpi0ana::reco_interaction_info(obj);
	    return s.pi0_subleading_photon_start_dedx;
	}

    /**
     * @brief Variable for pi0 subleading photon conversion distance.
     * @details Variable for pi0 subleading photon conversion distance
     * [cm], as calculated using interaction vertex and shower start point.
     * @tparam T the type of interaction (true or reco). 
     * @param obj the interaction to apply the variable on.
     * @return the pi0 subleading photon conversion distance. 
     */
    template<class T>
      double pi0_subleading_photon_conv_dist(const T & obj)
      {
	  if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
		         {
			     truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);
			     return s.pi0_subleading_photon_conv_dist;
		         }
	  else
          {
	      reco_inter s = utilities_ncpi0ana::reco_interaction_info(obj);
	      return s.pi0_subleading_photon_conv_dist;
          }
      }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pi0_subleading_photon_conv_dist, pi0_subleading_photon_conv_dist);

    template<class T>
      double pi0_subleading_photon_ip(const T & obj)
      {
	if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
		       {
			 truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);
			 return -5;
		       }
	else
	  {
	    reco_inter s = utilities_ncpi0ana::reco_interaction_info(obj);
	    return s.pi0_subleading_photon_ip;
	  }
      }


    /**
     * @brief Variable for angle between pi0 subleading photon cluster direction and vertex direction.
     * @details Variable for angle between pi0 subleading photon cluster direction, calculated as the normalized
     * mean direction of the cluster, and the the pi0 subleading photon vertex direction, calculated from the 
     * vector between the interaction vertex and shower start point.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the angle between the pi0 subleading photon cluster direction and vertex direction.
     */
    template<class T>
        double pi0_subleading_photon_cosphi(const T & obj)
        {
	    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
			   {
			       truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);
			       return -5;
			       //return s.pi0_subleading_photon_cosphi;
			   }
	    else
	    {
	      reco_inter s = utilities_ncpi0ana::reco_interaction_info(obj);
	      return s.pi0_subleading_photon_cosphi;
	    }
        }

    /**
     * @brief Variable for neutral pion momentum magnitude.
     * @details Variable for momentum of the neutral pion
     * candidate [MeV/c].
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the neutral pion momentum.
     */
    template<class T>
        double pi0_momentum_mag(const T & obj)
        {
	    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
		       {
			   truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);
			   return s.pi0_momentum_mag;
		       }
	    else
	    {
	        reco_inter s = utilities_ncpi0ana::reco_interaction_info(obj);
		return s.pi0_momentum_mag;
	    } 
        }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pi0_momentum_mag, pi0_momentum_mag);

    template<class T>
        double pi0_photons_avg_ip(const T & obj)
        {
	    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
			   {
			       truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);
			       return -5;
			   }
            else
	    {
                reco_inter s = utilities_ncpi0ana::reco_interaction_info(obj);
                return s.pi0_photons_avg_ip;
	    }
	}

    /**
     * @brief Variable for neutral pion angle with beam.
     * @details Variable for the cosine of the angle between
     * the interaction's neutral pion and the beam.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the neutral pion angle w.r.t. beam.
     */
    template<class T>
        double pi0_beam_costheta(const T & obj)
        {
	    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
		       {
			   truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);
			   return s.pi0_beam_costheta;
		       }
	    else
	    {
	        reco_inter s = utilities_ncpi0ana::reco_interaction_info(obj);
		return s.pi0_beam_costheta;
	    }
        }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pi0_beam_costheta, pi0_beam_costheta);

    /**
     * @brief Variable for neutral pion (photons) opening angle.
     * @details Variable for the opening angle between neutral pion
     * photons, as calculated using the interaction vertex and shower
     * start points.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     */
    template<class T>
        double pi0_photons_costheta(const T & obj)
        {
	    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
			   {
			       truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);
			       return s.pi0_photons_costheta;
			   }
	    else
	    {
	        reco_inter s = utilities_ncpi0ana::reco_interaction_info(obj);
		return s.pi0_photons_costheta;
	    }

	}
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pi0_photons_costheta, pi0_photons_costheta);

    /**
     * @brief Variable for neutral pion mass.
     * @details Variable for neutral pion mass [MeV/c^2], as calculated
     * with photon energies and opening angle.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the neutral pion mass.
     */
    template<class T>
        double pi0_mass(const T & obj)
        {
	    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                           {
			       truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);
                               return s.pi0_mass;
                           }
	    else
	    {
	        reco_inter s = utilities_ncpi0ana::reco_interaction_info(obj);
	        return s.pi0_mass;
	    }
	}
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pi0_mass, pi0_mass);

    template<class T>
    double num_primary_protons_thresh(const T & obj)
    {
        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
		       {
			   truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);
			   return s.num_primary_protons_thresh;
		       }
	else
	{
	    reco_inter s = utilities_ncpi0ana::reco_interaction_info(obj);
	    return s.num_primary_protons_thresh;
	}
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, num_primary_protons_thresh, num_primary_protons_thresh);

    template<class T>
    double num_primary_protons(const T & obj)
    {
        if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
		       {
			   truth_inter s = utilities_ncpi0ana::truth_interaction_info(obj);
			   return s.num_primary_protons;
		       }
	else
	{
	    reco_inter s = utilities_ncpi0ana::reco_interaction_info(obj);
	    return s.num_primary_protons;
	}
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, num_primary_protons, num_primary_protons);


    /**
     * @brief Variable for total visible energy of interaction.        
     * @details This function calculates the total visible energy of the   
     * interaction by summing the energy of all particles that are identified                  
     * as counting towards the final state of the interaction.                        
     * @tparam T the type of interaction (true or reco).                         
     * @param obj interaction to apply the variable on.                     
     * @return the total visible energy of the interaction.                                                                         
     */
    template<class T>
        double visible_energy(const T & obj)
        {
	    double energy(0);
	    for(const auto & p : obj.particles)
	    {
	        if(utilities_ncpi0ana::final_state_signal(p))
		{
		    energy += pvars::energy(p);
		    if(PIDFUNC(p) == 4) energy -= pvars::mass(p) - PROTON_BINDING_ENERGY;
		}
	    }
	    return energy/1000.0;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, visible_energy, visible_energy);

    /**
     * @brief Variable for the transverse momentum of the interaction counting
     * only particles identified as contributing to the final state.
     * @details This function calculates the transverse momentum of the
     * interaction by summing the transverse momentum of all particles that are
     * identified as counting towards the final state of the interaction. The
     * neutrino direction is assumed to either be the BNB axis direction
     * (z-axis) or the unit vector pointing from the NuMI target to the
     * interaction vertex. See @ref utilities::transverse_momentum for details
     * on the extraction of the transverse momentum.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the transverse momentum of the primary particles.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
        double dpT(const T & obj)
        {
	    utilities::three_vector pt = {0, 0, 0};
	    for(const auto & p : obj.particles)
	    {
	        if(utilities_ncpi0ana::final_state_signal(p))
		{
		    // Sum up the transverse momentum of all final state particles                                         
		    utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
		    utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
		    utilities::three_vector this_pt = utilities::transverse_momentum(momentum, vtx);
		    pt = utilities::add(pt, this_pt);
		}
	    }
	    return utilities::magnitude(pt)/1000.0;
	}
    REGISTER_VAR_SCOPE(RegistrationScope::Both, dpT, dpT);

    /**
     * @brief Variable for phi_T of the interaction.                                     
     * @details phi_T is a transverse kinematic imbalance variable defined                                
     * using the transverse momentum of the leading muon and the total hadronic                                        
     * system. This variable is sensitive to the presence of F.S.I. The                                       
     * neutrino direction is assumed to either be the BNB axis direction                                       
     * (z-axis) or the unit vector pointing from the NuMI target to the                                 
     * interaction vertex. See @ref utilities::transverse_momentum for details                                                           
     * on the extraction of the transverse momentum.                                                       
     * @tparam T the type of interaction (true or reco).                                                    
     * @param obj the interaction to apply the variable on.                                                                   
     * @return the phi_T of the interaction.                                  
     * @note The switch to the NuMI beam direction instead of the BNB axis is                                       
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).                                      
     */
    template<class T>
        double dphiT(const T & obj)
        {
	    utilities::three_vector lepton_pt = {0, 0, 0};
	    utilities::three_vector hadronic_pt = {0, 0, 0};
	    for(const auto & p : obj.particles)
	    {
	        if(utilities_ncpi0ana::final_state_signal(p))
		{
		    // There should only be one lepton, so replace the lepton                                        
		    // transverse momentum if the particle is a lepton.                                
		    utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
		    utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
		    utilities::three_vector this_pt = utilities::transverse_momentum(momentum, vtx);
		    if(PIDFUNC(p) == 1 || PIDFUNC(p) == 2)
		      lepton_pt = this_pt;
		    // The total hadronic system is treated as a single object.
		    else if(PIDFUNC(p) > 2)
		      hadronic_pt = utilities::add(hadronic_pt, this_pt);
		}
	    }
	    return std::acos(-1 * utilities::dot_product(lepton_pt, hadronic_pt) / (utilities::magnitude(lepton_pt) * utilities::magnitude(hadronic_pt)));
	}
    REGISTER_VAR_SCOPE(RegistrationScope::Both, dphiT, dphiT);

    /**
     * @brief Variable for alpha_T of the interaction.
     * @details alpha_T is a transverse kinematic imbalance variable defined
     * using the transverse momentum of the total hadronic system and the
     * outgoing lepton. The neutrino direction is assumed to either be the BNB
     * axis direction (z-axis) or the unit vector pointing from the NuMI target
     * to the interaction vertex. See @ref utilities::transverse_momentum for
     * details on the extraction of the transverse momentum.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the alpha_T of the interaction.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
        double dalphaT(const T & obj)
        {
	    utilities::three_vector lepton_pt = {0, 0, 0};
	    utilities::three_vector total_pt = {0, 0, 0};
	    for(const auto & p : obj.particles)
	    {
	        if(utilities_ncpi0ana::final_state_signal(p))
		{
		    // There should only be one lepton, so replace the lepton
		    // transverse momentum if the particle is a lepton.
		    utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
		    utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
		    utilities::three_vector this_pt = utilities::transverse_momentum(momentum, vtx);
		    if(PIDFUNC(p) == 1 || PIDFUNC(p) == 2)
		      lepton_pt = this_pt;
		    total_pt = utilities::add(total_pt, this_pt);
		}
	    }
	    return std::acos(-1 * utilities::dot_product(total_pt, lepton_pt) / (utilities::magnitude(total_pt) * utilities::magnitude(lepton_pt)));
	}
    REGISTER_VAR_SCOPE(RegistrationScope::Both, dalphaT, dalphaT);

}
#endif // VARS_NCPI0ANA_H
