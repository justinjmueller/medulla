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
	      double num_primary_pi0s_thresh = utilities::true_primary_pi0_multiplicity(obj, {params[3]});
	      double num_primary_pions_thresh = vars::pion_multiplicity(obj, {params[4]});
	      double num_primary_protons_thresh = vars::proton_multiplicity(obj, {params[5]});
	
	      double num_nonprimary_pi0s = utilities::true_nonprimary_pi0_multiplicity(obj, {0});

	      double num_primary_photons = vars::photon_multiplicity(obj, {0});
	      double num_primary_electrons = vars::electron_multiplicity(obj, {0});
        double num_primary_muons = vars::muon_multiplicity(obj, {0});
	      double num_primary_pi0s = utilities::true_primary_pi0_multiplicity(obj, {0.0});
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
	      double num_primary_pi0s_thresh = utilities::true_primary_pi0_multiplicity(obj, {params[3]});
	      double num_primary_pions_thresh = vars::pion_multiplicity(obj, {params[4]});
	      double num_primary_protons_thresh = vars::proton_multiplicity(obj, {params[5]});
	
	      double num_nonprimary_pi0s = utilities::true_nonprimary_pi0_multiplicity(obj, {0});
      
	      double num_primary_photons = vars::photon_multiplicity(obj, {0});
	      double num_primary_electrons = vars::electron_multiplicity(obj, {0});
	      double num_primary_muons = vars::muon_multiplicity(obj, {0});
	      double num_primary_pi0s = utilities::true_primary_pi0_multiplicity(obj, {0.0});
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
	      double num_primary_pi0s_thresh = utilities::true_primary_pi0_multiplicity(obj, {params[3]});
	      double num_primary_pions_thresh = vars::pion_multiplicity(obj, {params[4]});
	      double num_primary_protons_thresh = vars::proton_multiplicity(obj, {params[5]});

	      double num_nonprimary_pi0s = utilities::true_nonprimary_pi0_multiplicity(obj, {0});
      
	      double num_primary_photons = vars::photon_multiplicity(obj, {0});
	      double num_primary_electrons = vars::electron_multiplicity(obj, {0});
	      double num_primary_muons = vars::muon_multiplicity(obj, {0});
	      double num_primary_pi0s = utilities::true_primary_pi0_multiplicity(obj, {0.0});
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
	      double num_primary_pi0s_thresh = utilities::true_primary_pi0_multiplicity(obj, {params[3]});
	      double num_primary_pions_thresh = vars::pion_multiplicity(obj, {params[4]});
	      double num_primary_protons_thresh = vars::proton_multiplicity(obj, {params[5]});
	
	      double num_nonprimary_pi0s = utilities::true_nonprimary_pi0_multiplicity(obj, {0.0});
	
	      double num_primary_photons = vars::photon_multiplicity(obj, {0.0});
	      double num_primary_electrons = vars::electron_multiplicity(obj, {0.0});
	      double num_primary_muons = vars::muon_multiplicity(obj, {0});
	      double num_primary_pi0s = utilities::true_primary_pi0_multiplicity(obj, {0.0});
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
     * @brief Variable for enumerating interaction topologies.
     * @details This variable provides a basic categorization of interactions
     * using the following categories:
     * 0: NCpi0 signal (in-phase, fiducial, 0mu0pi1pi0)
     * 1: NCpi0 non-signal background
     * 2: CCpi0
	 * 3: Other nu without pi0.
     * 10: Cosmic
     * @param obj the interaction to apply the variable on.
     * @return the enumerated topology of the interaction.
     */
    template<class T>
    double category_topology_ncpi0_simple2(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params={})
    {
		double num_primary_pi0s = utilities::true_primary_pi0_multiplicity(obj, {0.0});

	    // Cosmic
	    uint16_t cat(10);
	
	    // Neutrino
	    if(cuts::neutrino(obj))
	    {
            // 0mu0pi1pi0 (in-phase, fiducial)
	        if(cuts::pi0ana::single_pi0<caf::SRInteractionTruthDLPProxy>(obj, {params[4]}) && cuts::no_muons(obj, {params[2]}) && cuts::no_charged_pions(obj, {params[3]}) 
			&& cuts::two_photons(obj, {params[0]}) && cuts::pi0ana::leading_photon_ke_cut(obj, {params[6]})
			&& !cuts::iscc(obj) && cuts::fiducial_cut(obj)) cat = 0;
	        // NCpi0 non-signal background
	        else if(num_primary_pi0s > 0 && !cuts::iscc(obj)) cat = 1;
			// CCpi0
	        else if(num_primary_pi0s > 0 && cuts::iscc(obj)) cat = 2;
	        // Other nu without pi0
	        else if(num_primary_pi0s == 0) cat = 3;
	    }
	    return cat;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, category_topology_ncpi0_simple2, category_topology_ncpi0_simple2);

	/**
     * @brief Variable for enumerating interaction topologies.
     * @details This variable provides a basic categorization of interactions
     * using the following categories:
     * 0: NCpi0 signal (in-phase, fiducial, 0mu0pi1pi0)
     * 1: NCpi0 non-signal background
     * 2: CCpi0
	 * 3: Other nu without pi0.
     * 10: Cosmic
     * @param obj the interaction to apply the variable on.
     * @return the enumerated topology of the interaction.
     */
    template<class T>
    double category_topology_ncpi0_simple3(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params={})
    {
		double num_primary_pi0s = utilities::true_primary_pi0_multiplicity(obj, {0.0});

		bool showers_in_fid_vol(false);
        size_t phi = selectors::pi0_leading_shower(obj);
		size_t dex = selectors::pi0_subleading_shower(obj);

	    if(phi != kNoMatch && dex != kNoMatch) 
		{
			auto & lead(obj.particles[phi]);
        	auto & sublead(obj.particles[dex]);
			showers_in_fid_vol = (pcuts::containment_cut(lead) && pcuts::containment_cut(sublead));
		}

	    // Cosmic
	    uint16_t cat(10);
	
	    // Neutrino
	    if(cuts::neutrino(obj))
	    {
            // 0mu0pi1pi0 (in-phase, fiducial)
	        if(cuts::pi0ana::single_pi0<caf::SRInteractionTruthDLPProxy>(obj, {params[4]}) && cuts::no_muons(obj, {params[2]}) && cuts::no_charged_pions(obj, {params[3]}) 
			&& cuts::two_photons(obj, {params[0]}) && cuts::pi0ana::leading_photon_ke_cut(obj, {params[6]}) && showers_in_fid_vol
			&& !cuts::iscc(obj) && cuts::fiducial_cut(obj)) cat = 0;
	        // NCpi0 non-signal background
	        else if(num_primary_pi0s > 0 && !cuts::iscc(obj)) cat = 1;
			// CCpi0
	        else if(num_primary_pi0s > 0 && cuts::iscc(obj)) cat = 2;
	        // Other nu without pi0
	        else if(num_primary_pi0s == 0) cat = 3;
	    }
	    return cat;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, category_topology_ncpi0_simple3, category_topology_ncpi0_simple3);

	/**
     * @brief Variable for enumerating interaction topologies.
     * @details This variable provides a basic categorization of interactions
     * using the following categories:
     * 0: NCpi0 single shower signal
	 * 1: CCpi0 single shower signal
     * 2: NCpi0 non-signal background
     * 3: CCpi0 non-signal background
	 * 4: Other nu without pi0.
     * 10: Cosmic
     * @param obj the interaction to apply the variable on.
     * @return the enumerated topology of the interaction.
     */
    template<class T>
    double category_topology_single_shower(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params={})
    {
		double num_primary_pi0s = utilities::true_primary_pi0_multiplicity(obj, {0.0});

	    // Cosmic
	    uint16_t cat(10);
	
	    // Neutrino
	    if(cuts::neutrino(obj))
	    {
            // NCpi0 single shower signal
	        if(cuts::no_muons(obj, {params[1]}) && cuts::no_charged_pions(obj, {params[2]}) && cuts::pi0ana::single_shower(obj, {params[0]}) 
			&& cuts::fiducial_cut(obj) && !cuts::iscc(obj)) cat = 0;
			// CCpi0 single shower signal
			else if(cuts::no_muons(obj, {params[1]}) && cuts::no_charged_pions(obj, {params[2]}) && cuts::pi0ana::single_shower(obj, {params[0]}) 
			&& cuts::fiducial_cut(obj) && cuts::iscc(obj)) cat = 1;
	        // NCpi0 non-signal background
	        else if(num_primary_pi0s > 0 && !cuts::iscc(obj)) cat = 2;
			// CCpi0
	        else if(num_primary_pi0s > 0 && cuts::iscc(obj)) cat = 3;
	        // Other nu without pi0
	        else if(num_primary_pi0s == 0) cat = 4;
	    }
	    return cat;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, category_topology_single_shower, category_topology_single_shower);

	/**
     * @brief Variable for enumerating interaction topologies.
     * @details This variable provides a basic categorization of interactions
     * using the following categories:
     * 0: NCpi0 single shower signal
	 * 1: CCpi0 single shower signal
     * 2: NCpi0 non-signal background
     * 3: CCpi0 non-signal background
	 * 4: Other nu without pi0.
     * 10: Cosmic
     * @param obj the interaction to apply the variable on.
     * @return the enumerated topology of the interaction.
     */
    template<class T>
    double category_topology_single_shower_v2(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params={})
    {
		double num_primary_pi0s = utilities::true_primary_pi0_multiplicity(obj, {0.0});

		bool shower_in_fid_vol(false);
        size_t phi = selectors::leading_primary_shower(obj);

	    if(phi != kNoMatch)
		{
			auto & lead(obj.particles[phi]);
			shower_in_fid_vol = pcuts::containment_cut(lead);
		}

	    // Cosmic
	    uint16_t cat(10);
	
	    // Neutrino
	    if(cuts::neutrino(obj))
	    {
            // NCpi0 single shower signal
	        if(cuts::no_muons(obj, {params[1]}) && cuts::no_charged_pions(obj, {params[2]}) && cuts::pi0ana::single_shower(obj, {params[0]}) && shower_in_fid_vol
			&& cuts::fiducial_cut(obj) && !cuts::iscc(obj)) cat = 0;
			// CCpi0 single shower signal
			else if(cuts::no_muons(obj, {params[1]}) && cuts::no_charged_pions(obj, {params[2]}) && cuts::pi0ana::single_shower(obj, {params[0]}) && shower_in_fid_vol
			&& cuts::fiducial_cut(obj) && cuts::iscc(obj)) cat = 1;
	        // NCpi0 non-signal background
	        else if(num_primary_pi0s > 0 && !cuts::iscc(obj)) cat = 2;
			// CCpi0
	        else if(num_primary_pi0s > 0 && cuts::iscc(obj)) cat = 3;
	        // Other nu without pi0
	        else if(num_primary_pi0s == 0) cat = 4;
	    }
	    return cat;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, category_topology_single_shower_v2, category_topology_single_shower_v2);

	/**
     * @brief Variable for enumerating the effect of truth level cuts for NCpi0 interactions.
     * @details This variable describes the effect of truth level cuts
     * using the following format:
	 * 5xxxxxxxx
	 * Where the numbers x in [0,3] represent the effect of each cut. 0 is failed individually
	 * and collectively, 1 is failed individually only (impossible), 2 is failed collectively
	 * only, 3 is passed.
     * @param obj the interaction to apply the variable on.
     * @return The enumerated effect of truth cuts.
     */
	template<class T>
    double ncpi0_truth_cut_study(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params={})
    {
		bool leading_shower_in_fid_vol(false);
		bool subleading_shower_in_fid_vol(false);
        size_t phi = selectors::pi0_leading_shower(obj);
		size_t dex = selectors::pi0_subleading_shower(obj);

	    if(phi != kNoMatch) 
		{
			auto & lead(obj.particles[phi]);
			leading_shower_in_fid_vol = pcuts::containment_cut(lead);
		}
		if(dex != kNoMatch)
		{
        	auto & sublead(obj.particles[dex]);
			subleading_shower_in_fid_vol = pcuts::containment_cut(sublead);
		}

	    // No cut
	    double cut_results = 5.0;
		bool run = true;
		// NCpi0
		bool ncpi0_cut = (cuts::pi0ana::single_pi0<caf::SRInteractionTruthDLPProxy>(obj, {params[4]}) && !cuts::iscc(obj));
		run = run && ncpi0_cut;
		cut_results = cut_results * 10 + ncpi0_cut*2 + run;
		// Fiducial
		bool fiducial_cut = cuts::fiducial_cut(obj);
		run = run && fiducial_cut;
		cut_results = cut_results * 10 + fiducial_cut*2 + run;
		// Leading in fiducial volume
		run = run && leading_shower_in_fid_vol;
		cut_results = cut_results * 10 + leading_shower_in_fid_vol*2 + run;
		// Subleading in fiducial volume
		run = run && subleading_shower_in_fid_vol;
		cut_results = cut_results * 10 + subleading_shower_in_fid_vol*2 + run;
		// No muons
		bool no_muons_cut = cuts::no_muons(obj, {params[1]});
		run = run && no_muons_cut;
		cut_results = cut_results * 10 + no_muons_cut*2 + run;
		// No charged pions
		bool no_charged_pions_cut = cuts::no_charged_pions(obj, {params[2]});
		run = run && no_charged_pions_cut;
		cut_results = cut_results * 10 + no_charged_pions_cut*2 + run;
		// Two photons
		bool two_photons_cut = cuts::two_photons(obj, {params[0]});
		run = run && two_photons_cut;
		cut_results = cut_results * 10 + two_photons_cut*2 + run;
		// Leading photon KE cut
		bool leading_photon_ke_cut = cuts::pi0ana::leading_photon_ke_cut(obj, {params[5]});
		run = run && leading_photon_ke_cut;
		cut_results = cut_results * 10 + leading_photon_ke_cut*2 + run;

		return cut_results;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, ncpi0_truth_cut_study, ncpi0_truth_cut_study);

}
#endif // VARS_PI0ANA_H