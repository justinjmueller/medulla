/**
 * @file vars_electron2025.h
 * @brief Header file for definitions of analysis variables specific to the
 * muon2024 analysis.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from interactions specific to the electron2025 benchmarking. 
 * Each variable is implemented as a function which takes an
 * interaction object as an argument and returns a double. These are the
 * building blocks for producing high-level plots of the selected interactions.
 * @author jzettle@fnal.gov
 */
#ifndef VARS_ELECTRON2025_H
#define VARS_ELECTRON2025_H

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"

#include "include/utilities.h"
#include "include/cuts.h"
#include "include/electron2025/cuts_electron2025.h"
#include "include/electron2025/utilities_electron2025.h"

#include <iostream>

/**
 * @namespace vars::electron2025
 * @brief Namespace for organizing variables specific to the muon2024 analysis.
 * @details This namespace is intended to be used for organizing variables which
 * act on interactions specific to the muon2024 analysis. Each variable is
 * implemented as a function which takes an interaction object as an argument
 * and returns a double. The function should be templated on the type of
 * interaction object if the variable is intended to be used on both true and
 * reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the vars
 * namespace, which is used for organizing generic variables which act on
 * interactions.
 */
namespace vars::electron2025
{
    /**
     * @brief Variable for enumerating interaction categories.
     * @details This variable provides a basic categorization of interactions
     * using only signal, neutrino background, and cosmic background as the
     * three categories.
     * 0: 1mu1p (contained and fiducial)
     * 1: 1mu1p (not contained or not fiducial)
     * 2: 1muNp (N > 1, contained and fiducial)
     * 3: 1muNp (N > 1, not contained or fiducial)
     * 4: 1muX (not 1muNp, contained and fiducial)
     * 5: 1muX (not 1muNp, not contained or fiducial)
     * 6: Other nu
     * 7: cosmic
     * @param obj The interaction to apply the variable on.
     * @return the enumerated category of the interaction.
    */

    double category(const caf::SRInteractionTruthDLPProxy & obj)
    {
        double cat(6);
        if(cuts::electron2025::all_2e_cut(obj)) cat = 0;
        else if(cuts::electron2025::all_1e1gamma_cut(obj)) cat = 1;
        else if(cuts::electron2025::all_1eNgamma_cut(obj)) cat = 2;
        else if(cuts::electron2025::all_2gamma_cut(obj)) cat = 3;
        else if(cuts::electron2025::all_gt2e_cut(obj)) cat = 4;
        else if(cuts::electron2025::all_gt2gamma_cut(obj)) cat = 5;
        //else if(cuts::neutrino(obj)) cat = 1;
        return cat;
    }

    double category_muons(const caf::SRInteractionTruthDLPProxy & obj)
    {
        double cat(7);
        if(cuts::electron2025::signal_1mu1p(obj)) cat = 0;
        else if(cuts::electron2025::nonsignal_1mu1p(obj)) cat = 1;
        else if(cuts::electron2025::signal_1muNp(obj)) cat = 2;
        else if(cuts::electron2025::nonsignal_1muNp(obj)) cat = 3;
        else if(cuts::electron2025::signal_1muX(obj)) cat = 4;
        else if(cuts::electron2025::nonsignal_1muX(obj)) cat = 5;
        else if(cuts::neutrino(obj)) cat = 6;
        return cat;
    }
    template<class T>
        double category_templated(const T & obj)
        {
            double cat(6);
            if(cuts::electron2025::all_2e_cut(obj)) cat = 0;
            else if(cuts::electron2025::all_1e1gamma_cut(obj)) cat = 1;
            else if(cuts::electron2025::all_1eNgamma_cut(obj)) cat = 2;
            else if(cuts::electron2025::all_2gamma_cut(obj)) cat = 3;
            else if(cuts::electron2025::all_gt2e_cut(obj)) cat = 4;
            else if(cuts::electron2025::all_gt2gamma_cut(obj)) cat = 5;
            //else if(cuts::neutrino(obj)) cat = 1;
            return cat;
        }

    /**
     * @brief Variable for the opening angle between leading muon and proton.
     * @details The leading muon and proton are defined as the particles with the
     * highest kinetic energy. The opening angle is defined as the arccosine of
     * the dot product of the momentum vectors of the leading muon and proton.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the opening angle between the leading muon and
     * proton.
     */
    template<class T>
        double opening_angle(const T & obj)
        {
            auto & m(obj.particles[utilities::leading_particle_index(obj, 2)]);
            auto & p(obj.particles[utilities::leading_particle_index(obj, 4)]);
            return std::acos(m.start_dir[0] * p.start_dir[0] + m.start_dir[1] * p.start_dir[1] + m.start_dir[2] * p.start_dir[2]);
        }

    template<class T>
        double opening_angle_ee(const T & obj)
        {
            std::vector<size_t> indices = utilities::electron2025::particle_indices(obj, 1, 0);
            if(indices[0] == indices[1]) return PLACEHOLDERVALUE; //return NaN if there is only one shower
            auto & e1(obj.particles[indices[0]]); //extract the leading shower
            auto & e2(obj.particles[indices[1]]); //extract the subleading shower 
            return std::acos(e1.start_dir[0] * e2.start_dir[0] + e1.start_dir[1] * e2.start_dir[1] + e1.start_dir[2] * e2.start_dir[2]);
        }

    template<class T>
        double visible_energy_ee(const T & obj)
        {
            double energy(0);
            for(const auto & p : obj.particles)
            {
                if(p.is_primary)
                {
                    if(pcuts::final_state_signal_elec(p)) //at the moment, we are only interested in electrons with energy > 25 MeV
                    {
                        energy += pvars::energy(p);
                    }
                }
            }
            return energy/1000.0;
        }

    template<class T>
        double leading_shower_energy(const T & obj)
        {
            std::vector<size_t> indices = utilities::electron2025::particle_indices(obj, 1, 0);
            auto & e1(obj.particles[indices[0]]); //extract the leading shower
            return pvars::energy(e1);
        }
    
    template<class T>
        double subleading_shower_energy(const T & obj)
        {
            std::vector<size_t> indices = utilities::electron2025::particle_indices(obj, 1, 0);
            if(indices[0] == indices[1]) return PLACEHOLDERVALUE; //return NaN if there is only one shower
            auto & e2(obj.particles[indices[1]]); //extract the subleading shower
            return pvars::energy(e2);
        }

    template<class T>
        double invariant_mass(const T & obj)
        {
            double inv_mass(0);
            std::vector<size_t> indices = utilities::electron2025::particle_indices(obj, 1, 0);
            if(indices[0] == indices[1]) return PLACEHOLDERVALUE; //return NaN if there is only one shower
            auto & e1(obj.particles[indices[0]]); //extract the leading shower
            auto & e2(obj.particles[indices[1]]); //extract the subleading shower
            double e1_energy(pvars::energy(e1));
            double e2_energy(pvars::energy(e2));
            double mag_1 = std::sqrt(std::pow(pvars::px(e1), 2) + std::pow(pvars::py(e1), 2) + std::pow(pvars::pz(e1), 2));
            double mag_2 = std::sqrt(std::pow(pvars::px(e2), 2) + std::pow(pvars::py(e2), 2) + std::pow(pvars::pz(e2), 2));
            double costheta = pvars::px(e1)/mag_1 * pvars::px(e2)/mag_2 + pvars::py(e1)/mag_1 * pvars::py(e2)/mag_2 + pvars::pz(e1)/mag_1 * pvars::pz(e2)/mag_2;
            //inv_mass = std::sqrt(2 * e1_energy * e2_energy * (1 - (costheta)));
            inv_mass = std::sqrt(e1.mass*e1.mass + e2.mass*e2.mass + 2 * ((e1_energy * e2_energy) - (pvars::px(e1) * pvars::px(e2) + pvars::py(e1) * pvars::py(e2) + pvars::pz(e1) * pvars::pz(e2))));
            return inv_mass;
        }

    /**
     * @brief Variable for the best-match IoU of the particle.
     * @details The best-match IoU is the intersection over union of the
     * points belonging to a pair of reconstructed and true particles. The
     * best-match IoU is calculated upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the best-match IoU of the particle.
     */
    template<class T>
        double iou(const T & p)
        {
            if(p.match_ids.size() > 0)
                return p.match_overlaps[0];
            else 
                return PLACEHOLDERVALUE;
        }

    template<class T>
        double nshowers(const T & obj)
        {
            std::vector<uint32_t> c(utilities::electron2025::count_primaries_ee(obj));
            return c[0]+c[1];
        }
    
    template<class T>
        double nelectrons(const T & obj)
        {
            std::vector<uint32_t> c(utilities::electron2025::count_primaries_ee(obj));
            return c[1];
        }

    template<class T>   
        double nphotons(const T & obj)
        {
            std::vector<uint32_t> c(utilities::electron2025::count_primaries_ee(obj));
            return c[0];
        }

}
#endif // VARS_ELECTRON2025_H