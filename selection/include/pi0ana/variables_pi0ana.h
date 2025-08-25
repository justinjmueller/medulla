/**
 * @file variables_pi0ana.h
 * @brief Header file for definitions of analysis variables specific to the
 * ccpi0 analysis.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from interactions specific to the ccpi0
 * analysis. Each variable is implemented as a function which takes an
 * interaction object as an argument and returns a double. These are the
 * building blocks for producing high-level plots of the selected interactions.
 * @author lkashur@colostate.edu
 */
#ifndef VARS_PI0ANA_H
#define VARS_PI0ANA_H

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"

#include <iostream>

#include "include/selectors.h"
#include "include/framework.h"

#include "include/cuts.h"
#include "include/pi0ana/cuts_pi0ana.h"

/**
 * @namespace vars::pi0ana
 * @brief Namespace for organizing variables specific to the ccpi0ana analysis.
 * @details This namespace is intended to be used for organizing variables which
 * act on interactions specific to the ccpi0ana analysis. Each variable is
 * implemented as a function which takes an interaction object as an argument
 * and returns a double. The function should be templated on the type of
 * interaction object if the variable is intended to be used on both true and
 * reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the vars
 * namespace, which is used for organizing generic variables which act on
 * interactions.
 */
namespace vars::pi0ana
{
   
    /**
     * @brief Variable for enumerating interaction topologies.
     * @details This variable provides a basic categorization of interactions
     * using the following categories:
     * 0: 1mu0pi1pi0 (in-phase, fiducial)
     * 1: Other nu-induced pi0
     * 2: Other nu without pi0
     * 4: Cosmic
     * @param obj the interaction to apply the variable on.
     * @return the enumerated topology of the interaction.
     */
    template<class T>
    double category_topology_ccpi0_simple1(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params={})
    {
	      double num_primary_photons_thresh = vars::photon_multiplicity(obj, {params[0]});
	      double num_primary_electrons_thresh = vars::electron_multiplicity(obj, {params[1]});
	      double num_primary_muons_thresh = vars::muon_multiplicity(obj, {params[2]});
	      double num_primary_pi0s_thresh = utilities_pi0ana::true_primary_pi0_multiplicity(obj, {params[3]});
	      double num_primary_pions_thresh = vars::pion_multiplicity(obj, {params[4]});
	      double num_primary_protons_thresh = vars::proton_multiplicity(obj, {params[5]});
	
	      double num_nonprimary_pi0s = utilities_pi0ana::true_nonprimary_pi0_multiplicity(obj, {0});

	      double num_primary_photons = vars::photon_multiplicity(obj, {0});
	      double num_primary_electrons = vars::electron_multiplicity(obj, {0});
        double num_primary_muons = vars::muon_multiplicity(obj, {0});
	      double num_primary_pi0s = utilities_pi0ana::true_primary_pi0_multiplicity(obj, {0.0});
	      double num_primary_pions = vars::pion_multiplicity(obj, {0});
	      double num_primary_protons = vars::proton_multiplicity(obj, {0});

        // Cosmic
	      uint16_t cat(10);

	      // Neutrino
	      if(cuts::neutrino(obj))
	      {
	          // 1mu0pi1pi0 (in-phase, fiducial)
	          if(num_primary_muons_thresh == 1 && num_primary_pions_thresh == 0 && num_primary_pi0s_thresh == 1 && cuts::iscc(obj) && cuts::fiducial_cut(obj)) cat = 0;
	          // Other nu-induced pi0
	          else if(num_primary_pi0s_thresh >= 1) cat = 1;
	          // Other nu without pi0
	          else if(num_primary_pi0s_thresh == 0) cat = 2;
	      }
        return cat;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, category_topology_ccpi0_simple1, category_topology_ccpi0_simple1);

    /**
     * @brief Variable for enumerating interaction topologies.
     * @details This variable provides a basic categorization of interactions
     * using the following categories:
     * 0: 1mu0pi1pi0 (in-phase, fiducial)
     * 1: 1mu 0pi (2+ pi0)
     * 2: 1mu Npi Xpi0
     * 3: 0mu Npi0
     * 4: Other nu
     * 10: Cosmic
     * @param obj the interaction to apply the variable on.
     * @return the enumerated topology of the interaction.
     */
    template<class T> 
    double category_topology_ccpi0_simple2(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params={})
    {
        double num_primary_photons_thresh = vars::photon_multiplicity(obj, {params[0]});
	      double num_primary_electrons_thresh = vars::electron_multiplicity(obj, {params[1]});
	      double num_primary_muons_thresh = vars::muon_multiplicity(obj, {params[2]});
	      double num_primary_pi0s_thresh = utilities_pi0ana::true_primary_pi0_multiplicity(obj, {params[3]});
	      double num_primary_pions_thresh = vars::pion_multiplicity(obj, {params[4]});
	      double num_primary_protons_thresh = vars::proton_multiplicity(obj, {params[5]});
	
	      double num_nonprimary_pi0s = utilities_pi0ana::true_nonprimary_pi0_multiplicity(obj, {0});
      
	      double num_primary_photons = vars::photon_multiplicity(obj, {0});
	      double num_primary_electrons = vars::electron_multiplicity(obj, {0});
	      double num_primary_muons = vars::muon_multiplicity(obj, {0});
	      double num_primary_pi0s = utilities_pi0ana::true_primary_pi0_multiplicity(obj, {0.0});
	      double num_primary_pions = vars::pion_multiplicity(obj, {0});
	      double num_primary_protons = vars::proton_multiplicity(obj, {0});

        // Cosmic
	      uint16_t cat(10);

	      // Neutrino
	      if(cuts::neutrino(obj))
	      {
	          // 1mu 0pi 1pi0 (in-phase, fiducial)
	          if(num_primary_muons_thresh == 1 && num_primary_pions_thresh == 0 && num_primary_pi0s_thresh == 1 && cuts::iscc(obj) && cuts::fiducial_cut(obj)) cat = 0;
	          // 1mu 0pi (2+ pi0)
	          else if(num_primary_muons_thresh == 1 && num_primary_pions_thresh == 0 && num_primary_pi0s_thresh >= 2 && cuts::iscc(obj) && cuts::fiducial_cut(obj)) cat = 1;
	          // 1mu Npi Xpi0
	          else if(num_primary_muons_thresh == 1 && num_primary_pions_thresh >= 1 && cuts::iscc(obj) && cuts::fiducial_cut(obj)) cat = 2;
	          // 0mu Npi0
	          else if(num_primary_muons_thresh == 0 && num_primary_pi0s_thresh >= 1 && !cuts::iscc(obj) && cuts::fiducial_cut(obj)) cat = 3;
	          // Other nu
	          else cat = 4;
	      }
	      return cat;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, category_topology_ccpi0_simple2, category_topology_ccpi0_simple2);
    
    /**
     * @brief Variable for enumerating interaction topologies.
     * @details This variable provides a basic categorization of interactions
     * using the following categories:
     * 0: 1mu0pi1pi0 (in-phase, fiducial)
     * 1: To-do...
     * 2: To-do...
     * 3: To-do...
     * @param obj the interaction to apply the variable on.
     * @return the enumerated topology of the interaction.
     */
    template<class T>
    double category_topology_ccpi0_complete(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params={})
    {
        double num_primary_photons_thresh = vars::photon_multiplicity(obj, {params[0]});
	      double num_primary_electrons_thresh = vars::electron_multiplicity(obj, {params[1]});
	      double num_primary_muons_thresh = vars::muon_multiplicity(obj, {params[2]});
	      double num_primary_pi0s_thresh = utilities_pi0ana::true_primary_pi0_multiplicity(obj, {params[3]});
	      double num_primary_pions_thresh = vars::pion_multiplicity(obj, {params[4]});
	      double num_primary_protons_thresh = vars::proton_multiplicity(obj, {params[5]});

	      double num_nonprimary_pi0s = utilities_pi0ana::true_nonprimary_pi0_multiplicity(obj, {0});
      
	      double num_primary_photons = vars::photon_multiplicity(obj, {0});
	      double num_primary_electrons = vars::electron_multiplicity(obj, {0});
	      double num_primary_muons = vars::muon_multiplicity(obj, {0});
	      double num_primary_pi0s = utilities_pi0ana::true_primary_pi0_multiplicity(obj, {0.0});
	      double num_primary_pions = vars::pion_multiplicity(obj, {0});
	      double num_primary_protons = vars::proton_multiplicity(obj, {0});

	      // Cosmic
	      uint16_t cat(10);

	      // Neutrino
	      if(cuts::neutrino(obj))
	      {
	          // 1mu 0pi 1pi0 (in-phase, fiducial)
	          if(num_primary_muons_thresh == 1 && num_primary_pions_thresh == 0 && num_primary_pi0s_thresh == 1 && cuts::iscc(obj) && cuts::fiducial_cut(obj)) cat = 0;
	          // 1mu 0pi 1pi0 (OOPS, fiducial)
	          else if( (num_primary_muons == 1 && num_primary_pions == 0 && num_primary_pi0s == 1 && cuts::iscc(obj) && cuts::fiducial_cut(obj)) && (num_primary_muons_thresh != 1 || num_primary_pions_thresh != 0 || num_primary_pi0s_thresh != 1) ) cat = 1;
	          // 1mu 0pi 1pi0 (OOFV)
	          else if(num_primary_muons_thresh == 1 && num_primary_pions_thresh == 0 && num_primary_pi0s_thresh == 1 && cuts::iscc(obj) && !cuts::fiducial_cut(obj)) cat = 2;
	          // 1mu 0pi (2+ pi0)
	          else if(num_primary_muons_thresh == 1 && num_primary_pions_thresh == 0 && num_primary_pi0s_thresh >= 2 && cuts::iscc(obj) && cuts::fiducial_cut(obj)) cat = 3;
	          // 1mu 0pi 0pi0 Npi0_nonprim
	          else if(num_primary_muons_thresh == 1 && num_primary_pions_thresh == 0 && num_primary_pi0s_thresh == 0 && num_nonprimary_pi0s >= 1 && cuts::iscc(obj) && cuts::fiducial_cut(obj)) cat = 4;
	          // 1mu 0pi 0pi0 0pi0_nonprim
	          else if(num_primary_muons_thresh == 1 && num_primary_pions_thresh == 0 && num_primary_pi0s_thresh == 0 && num_nonprimary_pi0s == 0 && cuts::iscc(obj) && cuts::fiducial_cut(obj)) cat = 5;
	          // 1mu Npi Xpi0
	          else if(num_primary_muons_thresh == 1 && num_primary_pions_thresh >= 1 && cuts::iscc(obj) && cuts::fiducial_cut(obj)) cat = 6;
	          // 0mu Npi0
	          else if(num_primary_muons_thresh == 0 && num_primary_pi0s_thresh >= 1 && !cuts::iscc(obj) && cuts::fiducial_cut(obj)) cat = 7;
	          // Other nu
	          else cat = 8;
	      }
	      return cat;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, category_topology_ccpi0_complete, category_topology_ccpi0_complete);

    /**
     * @brief Variable for enumerating interaction topologies.
     * @details This variable provides a basic categorization of interactions
     * using the following categories:
     * 0: 0mu0pi1pi0 (in-phase, fiducial)
     * 1: Other nu-induced pi0
     * 2: Other nu without pi0.
     * 10: Cosmic
     * @param obj the interaction to apply the variable on.
     * @return the enumerated topology of the interaction.
     */
    template<class T>
    double category_topology_ncpi0_simple1(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params={})
    {
        double num_primary_photons_thresh = vars::photon_multiplicity(obj, {params[0]});
	      double num_primary_electrons_thresh = vars::electron_multiplicity(obj, {params[1]});
	      double num_primary_muons_thresh = vars::muon_multiplicity(obj, {params[2]});
	      double num_primary_pi0s_thresh = utilities_pi0ana::true_primary_pi0_multiplicity(obj, {params[3]});
	      double num_primary_pions_thresh = vars::pion_multiplicity(obj, {params[4]});
	      double num_primary_protons_thresh = vars::proton_multiplicity(obj, {params[5]});
	
	      double num_nonprimary_pi0s = utilities_pi0ana::true_nonprimary_pi0_multiplicity(obj, {0.0});
	
	      double num_primary_photons = vars::photon_multiplicity(obj, {0.0});
	      double num_primary_electrons = vars::electron_multiplicity(obj, {0.0});
	      double num_primary_muons = vars::muon_multiplicity(obj, {0});
	      double num_primary_pi0s = utilities_pi0ana::true_primary_pi0_multiplicity(obj, {0.0});
	      double num_primary_pions = vars::pion_multiplicity(obj, {0});
	      double num_primary_protons = vars::proton_multiplicity(obj, {0});

	      // Cosmic
	      uint16_t cat(10);
	
	      // Neutrino
	      if(cuts::neutrino(obj))
	      {
            // 0mu0pi1pi0 (in-phase, fiducial)
	          if(num_primary_muons_thresh == 0 && num_primary_pions_thresh == 0 && num_primary_pi0s == 1 && !cuts::iscc(obj) && cuts::fiducial_cut(obj)) cat = 0;
	          // Other nu-induced pi0
	          else if(num_primary_pi0s >= 1) cat = 1;
	          // Other nu without pi0
	          else if(num_primary_pi0s == 0) cat = 2;
	      }
	      return cat;

    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, category_topology_ncpi0_simple1, category_topology_ncpi0_simple1);

    /**
     * @brief Dummy GUNDAM variable for cut type.
     * @details Cut type specifies a signal or sideband cut.
     * @param obj the interaction to apply the variable on.
     * @return the cut type for GUNDAM. 
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

    /**
     * @brief Dummy GUNDAM variable for sample classification (data or not).
     * @details "is_data" specifies whether or not sample comes from data.
     * @param obj the interaction to apply the variable on.
     * @return the sample classification (data or not.)
     */
    template<class T>
    double is_data(const T & obj)
    {
        double cat(-5);
	      return cat;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, is_data, is_data);

    /**
     * @brief Dummy GUNDAM variable for interaction classification (neutrino or not).
     * @details "is_nu" specifies wheter or not an interaction is from a neutrino.
     * @param obj the interactin to apply the variable on.
     * @return the interaction classification (nu or not.)
     */
    template<class T>
    double is_nu(const T & obj)
    {
        double cat(-5);
	      return cat;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, is_nu, is_nu);

    /**
     * @brief Variable for neutral pion momentum magnitude.
     * @details Momentum is calculated upstream in SPINE reconstruction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the momentum of the neutral pion in GeV.
     */
    template<class T>
    double pi0_momentum(const T & obj)
    {
        pi0 s = utilities_pi0ana::pi0_info(obj);
	      return s.momentum;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pi0_momentum, pi0_momentum);
    
    /**
     * @brief Variable for neutral pion angle to neutrino beam.
     * @details Angle is calculated using neutral pion momentum and beam direction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the neutral pion angle to the neutrino beam.
     */
    template<class T>
    double pi0_beam_costheta(const T & obj)
    {
        pi0 s = utilities_pi0ana::pi0_info(obj);
	      return s.beam_costheta;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pi0_beam_costheta, pi0_beam_costheta);

    /**
     * @brief Variable for pi0 leading photon conversion distance.
     * @details Variable for pi0 leading photon conversion distance
     * [cm], as calculated using interaction vertex and shower start point.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the pi0 leading photon conversion distance in cm.
     */
    template<class T>
    double pi0_leading_shower_conv_dist(const T & obj)
    {
        pi0 s = utilities_pi0ana::pi0_info(obj);
	      return s.leading_shower_conv_dist;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pi0_leading_shower_conv_dist, pi0_leading_shower_conv_dist);

    /**
     * @brief Variable for pi0 subleading photon conversion distance.
     * @details Variable for pi0 subleading photon conversion distance
     * [cm], as calculated using interaction vertex and shower start point.
     * @tparam T the type of interaction (true or reco). 
     * @param obj the interaction to apply the variable on.
     * @return the pi0 subleading photon conversion distance in cm. 
     */
    template<class T>
    double pi0_subleading_shower_conv_dist(const T & obj)
    {
        pi0 s = utilities_pi0ana::pi0_info(obj);
        return s.subleading_shower_conv_dist;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pi0_subleading_shower_conv_dist, pi0_subleading_shower_conv_dist);

    /**
    * @brief Variable for pi0 diphoton opening angle (cosine).
    * @Details Variable for pi0 diphoton opening angle, as calculated 
    * using the neutrino interaction vertex and shower start points in reco, 
    * and the true photon momenta in truth.
    * @tparam T the type of interaction (true or reco).
    * @param obj the interaction to apply the variable on.
    * @return the pi0 diphoton opening angle.
    */
    template<class T>
    double pi0_showers_costheta(const T & obj)
    {
	      pi0 s = utilities_pi0ana::pi0_info(obj);
	      return s.showers_costheta;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pi0_showers_costheta, pi0_showers_costheta);

    /**
     * @brief Variable for neutral pion mass.
     * @details Variable for neutral pion mass, as calculated
     * with photon energies and opening angle.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the neutral pion mass in MeV/c^2.
     */
    template<class T>
    double pi0_mass(const T & obj)
    {

        pi0 s = utilities_pi0ana::pi0_info(obj);
	      return s.mass;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pi0_mass, pi0_mass);
}
#endif // VARS_PI0ANA_H
