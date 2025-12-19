/**
 * @file mctruth.h
 * @brief Definitions of analysis variables which can extract information from
 * the SRTrueInteraction object.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from the SRTrueInteraction object. Each variable
 * is implemented as a function which takes an SRTrueInteraction object as an
 * argument and returns a double. The association of an SRInteractionTruthDLP
 * object to an SRTrueInteraction object is handled upstream in the SpineVar
 * functions.
 * @author mueller@fnal.gov
 */
#ifndef MCTRUTH_H
#define MCTRUTH_H
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"

#include "framework.h"

/**
 * @namespace mctruth
 * @brief Namespace for organizing variables which act on true interactions.
 * @details This namespace is intended to be used for organizing variables
 * which act on true interactions. Each variable is implemented as a function
 * which takes an SRTrueInteraction object as an argument and returns a double.
 */
namespace mctruth
{
    /**
     * @brief Variable for the true neutrino energy.
     * @details This variable is intended to provide the true energy of the
     * parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino energy.
     */
    template<typename T>
        double neutrino_energy(const T & obj) { return obj.E; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, neutrino_energy, neutrino_energy);

    /**
     * @brief Variable for momentum transfer squared.
     * @details This variable is intended to provide the true
     * square of the momentum transfered from the neutrino to 
     * the nucleus.
     * @tparam T the type of object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true momentum transfer squared.
     */
    template<typename T>
        double neutrino_Q2(const T & obj) { return obj.Q2; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, neutrino_Q2, neutrino_Q2);

    /**
     * @brief Variable for hadronic invariant mass.
     * @details This variable is intended to provide the true
     * hadronic invariant mass of the neutrino interaction.
     * @tparam T the type of object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true hadronic invariant mass.
     */
    template<typename T>
        double neutrino_W(const T & obj) { return obj.w; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, neutrino_W, neutrino_W);

    /**
     * @brief Variable for the true neutrino baseline.
     * @details This variable is intended to provide the true baseline of the
     * parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino baseline.
     */
    template<typename T>
        double baseline(const T & obj) { return obj.baseline; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, baseline, baseline);

    /**
     * @brief Variable for the true neutrino PDG code.
     * @details This variable is intended to provide the true PDG code of the
     * parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino PDG code.
     */
    template<typename T>
        double pdg(const T & obj) { return obj.pdg; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pdg, pdg);

    /**
     * @brief Variable for the PDG code of the parent of the neutrino.
     * @details This variable is intended to provide the PDG code of the
     * parent of the neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the PDG code of the parent of the neutrino.
     */
    template<typename T>
        double parent_pdg(const T & obj) { return obj.parent_pdg; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, parent_pdg, parent_pdg);

    /**
     * @brief Variable for the true neutrino current value.
     * @details This variable is intended to provide the true current value of
     * the parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino current value.
     */
    template<typename T>
        double cc(const T & obj) { return obj.iscc; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, cc, cc);

    /**
     * @brief Variable for the interaction mode of the interaction.
     * @details This variable is intended to provide the interaction mode of the
     * interaction. This is based on the GENIE interaction mode enumeration 
     * defined in the LArSoft MCNeutrino class.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the interaction mode.
     */
    template<typename T>
        double interaction_mode(const T & obj) { return obj.genie_mode; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, interaction_mode, interaction_mode);

    /**
     * @brief Variable for the interaction type of the interaction.
     * @details This variable is intended to provide the interaction type of the
     * interaction. This is based on the GENIE interaction type enumeration 
     * defined in the LArSoft MCNeutrino class.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the interaction type.
     */
    template<typename T>
        double interaction_type(const T & obj) { return obj.genie_inttype; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, interaction_type, interaction_type);

    /**
     * @brief Variable for the number of true final state neutral pions in the interaction.
     * @details Variable for the number of true primary neutral pions in each SRTrueInteraction.
     * @tparam T the type of object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the number of true final state neutral pions in the interaction.
     */
    template<typename T>
      double npi0s_srtruth(const T & obj, std::vector<double> params={0.0,})
      {
	  int num_pi0s(0);

	  // Loop over primary particles
	  for(const auto & p : obj.prim)
	  {
	      // Check pi0 pdg_code
	      if(p.pdg == 111)
	      {
		// Check KE threshold
		double ke(-5);
                ke = 1000. * (p.genE - (PI0_MASS/1000.)); // MeV
		if(ke >= params[0])
		  num_pi0s++;
	      }
	  }
	  return num_pi0s;
      }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, npi0s_srtruth, npi0s_srtruth);

    /**
     * @brief Variable for the number of true final muons in the interaction.
     * @details Variable for the number of true primary muons in each SRTrueInteraction.
     * @tparam T the type of object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the number of true final state muons.
     */
    template<typename T>
      double nmuons_srtruth(const T & obj, std::vector<double> params={0.0,})
      {
	  int num_muons(0);

	  // Loop over primary particles
	  for(const auto & p : obj.prim)
	  {
	      // Check muon pdg_code 
	    if(p.pdg == 13)
	    {
		// Check KE threshold
		double ke(-5);
		ke = 1000. * (p.genE - (MUON_MASS/1000.)); // MeV
		if(ke >= params[0])
		  num_muons++;
	    }
          }
	  return num_muons;
      }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nmuons_srtruth, nmuons_srtruth);
    
    /**
     * @brief Variable for the number of true final pions in the interaction.
     * @details Variable for the number of true primary pions in each SRTrueInteraction.
     * @tparam T the type of object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the number of true final state pions.
     */
    template<typename T>
      double npions_srtruth(const T & obj, std::vector<double> params={0.0,})
      {
	  int num_pions(0);

	  // Loop over primary particles 
	  for(const auto & p : obj.prim)
          {
	      // Check pion pdg_code
	      if(abs(p.pdg) == 211)
              {
                  // Check KE threshold
                  double ke(-5);
		  ke = 1000. * (p.genE - (PION_MASS/1000.)); // MeV
		  if(ke >= params[0])
		    num_pions++;
              }
          }
	  return num_pions;
      }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, npions_srtruth, npions_srtruth);

    /**
     * @brief Variable for the momentum of the true muon.
     * @details Variable for the momentum magnitude of the muon in true 1mu interactions.
     * @tparam T the type of object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the momentum of the true muon.
     */
    template<typename T>
      double muon_p(const T & obj, std::vector<double> params={0.0,})
      {
	  int num_muons = nmuons_srtruth(obj, params);

	  TVector3 mom(0,0,0);
	  double mom_mag(-5);
	  if(num_muons == 1)
          {
	      for(const auto & p : obj.prim)
              {
		  double ke(-5);
		  ke = 1000. * (p.genE - (MUON_MASS/1000.)); // MeV
		  if(p.pdg == 13 && ke >= params[0])
                  {
		      mom.SetX(p.genp.x);
		      mom.SetY(p.genp.y);
		      mom.SetZ(p.genp.z);
		      mom_mag = mom.Mag();
                  }
              }
          }
	  return mom_mag;
      }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, muon_p, muon_p);

   /**
     * @brief Variable for the momentum of the true neutral pion.
     * @details Variable for the momentum magnitude of the neutral pion in true 1pi0 interactions.
     * @tparam T the type of object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the momentum of the true neutral pion.
     */
    template<typename T>
      double pi0_p(const T & obj, std::vector<double> params={0.0,})
      {
	  int num_pi0s = npi0s_srtruth(obj, params);
	  
	  TVector3 mom(0,0,0);
	  double mom_mag(-5);
	  if(num_pi0s == 1)
	  {
	      for(const auto & p : obj.prim)
	      {

		  double ke(-5);
		  ke = 1000. * (p.genE - (PI0_MASS/1000.)); // MeV
		  if(p.pdg == 111 && ke >= params[0])
		  {
		      mom.SetX(p.genp.x);
		      mom.SetY(p.genp.y);
		      mom.SetZ(p.genp.z);

		      mom_mag = mom.Mag();
		  }
	      }
	  }
	  return mom_mag;
      }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pi0_p, pi0_p);
    
    /**
     * @brief Variable for the angle the true neutral pion makes with the beam.
     * @details Variable for the cosine of the true neutral pion angle with
     * respect to the neutrino beam.
     * @tparam T the type of object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the cosine of the angle between the true neutral pion and the beam.
     */
    template<typename T>
      double pi0_nu_costheta(const T & obj, std::vector<double> params={0.0,})
      {
	  int num_pi0s = npi0s_srtruth(obj, params);
	
	  TVector3 mom(0,0,0);
	  if(num_pi0s == 1)
	  {
	      for(const auto & p : obj.prim)
              {
		  double ke(-5);
		  ke = 1000. * (p.genE - (PI0_MASS/1000.)); // MeV 
		  if(p.pdg == 111 && ke >= params[0])
                  {
		      mom.SetX(p.genp.x);
		      mom.SetY(p.genp.y);
		      mom.SetZ(p.genp.z);
                  }
              }
	  }
	  
	  TVector3 nu_dir(obj.momentum.x, obj.momentum.y, obj.momentum.z);
	  double costheta = mom.Unit().Dot(nu_dir.Unit());
	  return costheta;
	
      }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pi0_nu_costheta, pi0_nu_costheta);

    /**
     * @brief Variable for the angle the true muon makes with the beam.
     * @details Variable for the cosine of the true muon angle with
     * respect to the neutrino beam.
     * @tparam T the type of object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the cosine of the angle between the true muon and the beam.
     */
    template<typename T>
      double muon_nu_costheta(const T & obj, std::vector<double> params={0.0,})
      {
	int num_muons = nmuons_srtruth(obj, params);

	TVector3 mom(0,0,0);
	if(num_muons == 1)
        {
	    for(const auto & p : obj.prim)
	    {
	        double ke(-5);
		ke = 1000. * (p.genE - (MUON_MASS/1000.)); // MeV 
		if(p.pdg == 13 && ke >= params[0])
                {
		    mom.SetX(p.genp.x);
		    mom.SetY(p.genp.y);
		    mom.SetZ(p.genp.z);
		}
	    }
	}

	TVector3 nu_dir(obj.momentum.x, obj.momentum.y, obj.momentum.z);
	double costheta = mom.Unit().Dot(nu_dir.Unit());
	return costheta;

      }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, muon_nu_costheta, muon_nu_costheta);

    /**
     * @brief Variable for true fiducial status of neutrino.
     * @details This variable states whether the true neutrino
     * interaction vertex is within the ICARUS fiducial volume.
     * @tparam T the type of object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the fiducial status of the neutrino.
     */
    template<typename T>
      double fiducial(const T & obj)
      {
	  bool is_fiducial(true);
	  float vtx_x = obj.position.x;
	  float vtx_y = obj.position.y;
	  float vtx_z = obj.position.z;
	  

	  // drift requirement, east
	  if( vtx_x < 0 && (vtx_x < -358.49 + 25 || vtx_x > -61.94 - 25) )
	    is_fiducial = false;
	  
	  // drift requirement, west
	  if( vtx_x > 0 && (vtx_x < 61.94 + 25 || vtx_x > 358.49 - 25) )
	    is_fiducial = false;

	  // verteical requirement
	  if( vtx_y < -181.86 + 25 || vtx_y > 134.96 - 25 )
	    is_fiducial = false;

	  // beam requirement
	  if( vtx_z < -894.9505 + 30 || vtx_z > 894.9505 - 50 )
	    is_fiducial = false;

	  // mystery z-gap
	  if( vtx_z > -100 && vtx_z < 100 )
	    is_fiducial = false;

	  // hanging cable
	  if( vtx_x > 210.215 && vtx_y > 60 && (vtx_z > 290 && vtx_z < 390) )
	    is_fiducial = false;

	  return is_fiducial;
      }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, fiducial, fiducial);

    /**
     * @brief Variable for enumerating interaction topologies.
     * @details This variable provides a basic categorization of interactions
     * using the following categories:
     * 0: 1mu0pi1pi0 (in-phase, fiducial)
     * 1: To-do...
     * 2: To-do...
     * 3: To-do...
     * @tparam T the type of object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the enumerated topology of the interaction.
     */
    template<typename T>
      double category_topology_ccpi0_simple2(const T & obj, std::vector<double> params={0.0,})
      {
	
	int nmuons_thresh = nmuons_srtruth(obj, {params[2]});
	int npi0s_thresh = npi0s_srtruth(obj, {params[3]});
	int npions_thresh = npions_srtruth(obj, {params[4]});

	uint16_t cat(10);

	// 1mu 0pi 1pi0 (in-phase, fiducial)
	if(cc(obj) && pdg(obj) == 14 && nmuons_thresh == 1 && npions_thresh == 0 && npi0s_thresh == 1 && fiducial(obj)) cat = 0;
	
	// 1mu 0pi (2+ pi0)
	else if(cc(obj) && pdg(obj) == 14 && nmuons_thresh == 1 && npions_thresh == 0 && npi0s_thresh >= 2 && fiducial(obj)) cat = 1;

	// 1mu Npi Xpi0
	else if(cc(obj) && pdg(obj) == 14 && nmuons_thresh == 1 && npions_thresh >= 1 && fiducial(obj)) cat = 2;

	// 0mu Npi0
	else if(!cc(obj) && nmuons_thresh == 0 && npi0s_thresh >= 1 && fiducial(obj)) cat = 3;

	// Other nu
	else cat = 4;

	return cat;
      }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, category_topology_ccpi0_simple2, category_topology_ccpi0_simple2);

} // namespace mctruth
#endif
