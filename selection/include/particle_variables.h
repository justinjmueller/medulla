/**
 * @file particle_variables.h
 * @brief Header file for definitions of variables which act on single
 * particles.
 * @details This file contains definitions of variables which act on single
 * particles. Each variable is implemented as a function which takes a particle
 * object as an argument and returns a double. These variables are intended to
 * be used to define more complex variables which act on interactions.
 * @author mueller@fnal.gov
*/
#ifndef PARTICLE_VARIABLES_H
#define PARTICLE_VARIABLES_H
#define ELECTRON_MASS 0.5109989461
#define MUON_MASS 105.6583745
#define PION_MASS 139.57039
#define PROTON_MASS 938.2720813

#include "include/particle_utilities.h"
#include "scorers.h"

/**
 * @namespace pvars
 * @brief Namespace for organizing generic variables which act on single
 * particles.
 * @details This namespace is intended to be used for organizing variables which
 * act on single particles. Each variable is implemented as a function which
 * takes a particle object as an argument and returns a double. The function
 * should be templated on the type of particle object if the variable is
 * intended to be used on both true and reconstructed particles.
 * @note The namespace is intended to be used in conjunction with the
 * vars namespace, which is used for organizing variables which act on
 * interactions.
 */
namespace pvars
{
    /**
     * @brief Variable for the particle's primary classification.
     * @details This variable returns the primary classification of the particle.
     * The primary classification is determined upstream in the SPINE
     * reconstruction and is based on the softmax scores of the particle.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the primary classification of the particle.
     */
    template<class T>
    double primary_classification(const T & p)
    {
        if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
            return p.is_primary ? 1 : 0;
        else
            return (*primfn)(p);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, primary_classification, primary_classification);
    
    /**
     * @brief Variable for the particle's PID.
     * @details This variable returns the PID of the particle. The PID is
     * determined by the softmax scores of the particle. This function uses the
     * configured PID function, which can be set by the user in the
     * configuration file.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the PID of the particle.
     */
    template<class T>
    double pid(const T & p)
    {
        if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
            return p.pid;
        else
            return (*pidfn)(p);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, pid, pid);

    /**
     * @brief Variable for the semantic type of the particle.
     * @details This variable returns the semantic type of the particle. The
     * semantic type is determined by majority-vote of the pixel-level semantic
     * segmentation of the particle. The semantic types are defined as follows:
     * 0: shower, 1: track, 2: Michel electron, 3: delta electron,
     * 4: low-energy, 5: ghost, and -1: unknown.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the semantic type of the particle.
     */
    template<class T>
    double semantic_type(const T & p)
    {
        return p.shape;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, semantic_type, semantic_type);

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
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, iou, iou);

    /**
     * @brief Variable for the containment status of the particle.
     * @details The containment status is determined upstream in the SPINE
     * reconstruction and is based on the set of all points in the particle,
     * which must be contained within the volume of the TPC that created them.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the containment status of the particle.
     */
    template<class T>
    double containment(const T & p) { return p.is_contained; }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, containment, containment);

    /**
     * @brief Variable for the mass of the particle.
     * @details The mass of the particle is determined by the PID of the
     * particle. This couples the PID to the mass of the particle, so it is
     * necessary to use the appropriate PID function rather than the in-built
     * PID attribute.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the mass of the particle.
     */
    template<class T>
    double mass(const T & p)
    {
        double mass(0);
        if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
        {
            mass = p.mass;
        }
        else
        {
            switch(int(pvars::pid(p)))
            {
                case 0:
                    mass = 0;
                    break;
                case 1:
                    mass = ELECTRON_MASS;
                    break;
                case 2:
                    mass = MUON_MASS;
                    break;
                case 3:
                    mass = PION_MASS;
                    break;
                case 4:
                    mass = PROTON_MASS;
                    break;
                default:
                    mass = PLACEHOLDERVALUE;
                    break;
            }
        }
        return mass;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, mass, mass);

    /**
     * @brief Variable for the CSDA kinetic energy of the particle.
     * @details The CSDA kinetic energy is calculated upstream in the SPINE
     * reconstruction by relating the length of the track to the energy loss
     * of the particle.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the CSDA kinetic energy of the particle.
     * @note This function is only valid for particles which are contained in
     * the detector.
     */
    template<class T>
    double csda_ke(const T & p)
    {
        double pidx(pvars::pid(p));
        if(pidx < 0 || std::isinf(p.csda_ke_per_pid[pidx]))
            return PLACEHOLDERVALUE;
        return p.csda_ke_per_pid[pidx];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, csda_ke, csda_ke);

    /**
     * @brief Variable for the MCS kinetic energy of the particle.
     * @details The MCS kinetic energy is calculated upstream in the SPINE
     * reconstruction by comparing the distribution of successive angles of
     * track segments to the expected distribution for multiple scattering.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the MCS kinetic energy of the particle.
     */
    template<class T>
    double mcs_ke(const T & p)
    {
        size_t pidx(pvars::pid(p));
        return std::isinf(p.mcs_ke_per_pid[pidx]) ? PLACEHOLDERVALUE : (double)p.mcs_ke_per_pid[pidx];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, mcs_ke, mcs_ke);

    /**
     * @brief Variable for the calorimetric kinetic energy of the particle.
     * @details The calorimetric kinetic energy is calculated upstream in the
     * SPINE reconstruction as the sum of the energy of each spacepoint in the
     * particle.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the calorimetric kinetic energy of the particle.
     */
    template<class T>
    double calo_ke(const T & p)
    {
        return p.calo_ke;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, calo_ke, calo_ke);

    /**
     * @brief Variable for true particle starting kinetic energy.
     * @details The starting kinetic energy is defined as the total energy
     * minus the rest mass energy of the particle.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the starting kinetic energy of the particle.
     */
    template<class T>
    double ke(const T & p)
    {
        double energy(0);
        if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
        {
            energy = p.energy_init - mass(p);
        }
        else
        {
            if(pvars::pid(p) < 2) [[likely]]
                energy += calo_ke(p);
            else
            {
                if(p.is_contained) energy += csda_ke(p);
                else energy += mcs_ke(p);
            }
        }
        return energy;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, ke, ke);

    /**
     * @brief Variable for the best estimate of the particle energy.
     * @details At the most basic decision level, this is based on the
     * shower/track designation. Showers can only be reconstructed
     * calorimetrically, while tracks can be reconstructed calorimetrically,
     * by range (if contained), or by multiple scattering (if exiting).
     * @tparam T the type of particle.
     * @param p the particle to apply the variable on.
     * @return the best estimate of the particle energy.
     */
    template<class T>
    double energy(const T & p)
    {
        double energy = ke(p) + mass(p);
        return energy;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, energy, energy);

    /**
     * @brief Variable for the length of the particle track.
     * @details The length of the track is calculated upstream in the SPINE
     * reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the length of the particle track.
     */
    template<class T>
    double length(const T & p)
    {
        return p.length;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, length, length);

    /**
     * @brief Variable for the angle of the particle in the x-wire plane for
     * a horizontal anode plane.
     * @details The angle is calculated as the angle of the particle start
     * direction within the plane defined by the x-axis and the wire direction
     * of a horizontal anode plane (e.g. the front induction plane of ICARUS).
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the angle of the particle in the x-wire plane.
     */
    template<class T>
    double theta_xw_horizontal(const T & p)
    {
        // Wire orientation
        double wx(0), wy(1), wz(0);
        return std::acos(p.start_dir[0] * wx + p.start_dir[1] * wy + p.start_dir[2] * wz);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, theta_xw_horizontal, theta_xw_horizontal);

    /**
     * @brief Variable for the angle of the particle in the x-wire plane for
     * an anode plane with wires oriented along an angle +60 degrees from the
     * horizontal.
     * @details The angle is calculated as the angle of the particle start
     * direction within the plane defined by the x-axis and the wire direction
     * of an anode plane with wires oriented along an angle of +60 degrees from
     * the horizontal (e.g. the second induction / collection plane of ICARUS).
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the angle of the particle in the x-wire plane.
     */
    template<class T>
    double theta_xw_p60(const T & p)
    {
        // Wire orientation
        double wx(0), wy(0.5 * std::sqrt(3)), wz(0.5);
        return std::acos(p.start_dir[0] * wx + p.start_dir[1] * wy + p.start_dir[2] * wz);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, theta_xw_p60, theta_xw_p60);
    
    /**
     * @brief Variable for the angle of the particle in the x-wire plane for
     * an anode plane with wires oriented along an angle -60 degrees from the
     * horizontal.
     * @details The angle is calculated as the angle of the particle start
     * direction within the plane defined by the x-axis and the wire direction
     * of an anode plane with wires oriented along an angle of -60 degrees from
     * the horizontal (e.g. the third induction / collection plane of ICARUS).
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the angle of the particle in the x-wire plane.
     */
    template<class T>
    double theta_xw_m60(const T & p)
    {
        // Wire orientation
        double wx(0), wy(-0.5 * std::sqrt(3)), wz(0.5);
        return std::acos(p.start_dir[0] * wx + p.start_dir[1] * wy + p.start_dir[2] * wz);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, theta_xw_m60, theta_xw_m60);
    
    /**
     * @brief Variable for the angle of the particle in the x-wire plane for
     * a vertical anode plane.
     * @details The angle is calculated as the angle of the particle start
     * direction within the plane defined by the x-axis and the wire direction
     * of a vertical anode plane (e.g. the collection plane of SBND).
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the angle of the particle in the x-wire plane.
     */
    template<class T>
    double theta_xw_vertical(const T & p)
    {
        // Wire orientation
        double wx(0), wy(1), wz(0);
        return std::acos(p.start_dir[0] * wx + p.start_dir[1] * wy + p.start_dir[2] * wz);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, theta_xw_vertical, theta_xw_vertical);

    /**
     * @brief Variable for the angle of the particle in the x-wire plane for
     * an anode plane with wires oriented along an angle +30 degrees from the
     * horizontal.
     * @details The angle is calculated as the angle of the particle start
     * direction within the plane defined by the x-axis and the wire direction
     * of an anode plane with wires oriented along an angle of +30 degrees from
     * the horizontal (e.g. the first induction plane of SBND).
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the angle of the particle in the x-wire plane.
     */
    template<class T>
    double theta_xw_p30(const T & p)
    {
        // Wire orientation
        double wx(0), wy(0.5), wz(0.5 * std::sqrt(3));
        return std::acos(p.start_dir[0] * wx + p.start_dir[1] * wy + p.start_dir[2] * wz);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, theta_xw_p30, theta_xw_p30);

    /**
     * @brief Variable for the angle of the particle in the x-wire plane for
     * an anode plane with wires oriented along an angle -30 degrees from the
     * horizontal.
     * @details The angle is calculated as the angle of the particle start
     * direction within the plane defined by the x-axis and the wire direction
     * of an anode plane with wires oriented along an angle of -30 degrees from
     * the horizontal (e.g. the second induction plane of SBND).
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the angle of the particle in the x-wire plane.
     */
    template<class T>
    double theta_xw_m30(const T & p)
    {
        // Wire orientation
        double wx(0), wy(-0.5), wz(0.5 * std::sqrt(3));
        return std::acos(p.start_dir[0] * wx + p.start_dir[1] * wy + p.start_dir[2] * wz);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, theta_xw_m30, theta_xw_m30);

    /**
     * @brief Variable for the x-coordinate of the particle starting point.
     * @details The starting point is the point at which the particle is created
     * and is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the x-coordinate of the particle starting point.
     */
    template<class T>
    double start_x(const T & p)
    {
        return p.start_point[0];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, start_x, start_x);

    /**
     * @brief Variable for the y-coordinate of the particle starting point.
     * @details The starting point is the point at which the particle is created
     * and is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the y-coordinate of the particle starting point.
     */
    template<class T>
    double start_y(const T & p)
    {
        return p.start_point[1];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, start_y, start_y);

    /**
     * @brief Variable for the z-coordinate of the particle starting point.
     * @details The starting point is the point at which the particle is created
     * and is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the z-coordinate of the particle starting point.
     */
    template<class T>
    double start_z(const T & p)
    {
        return p.start_point[2];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, start_z, start_z);

    /**
     * @brief Variable for the x-coordinate of the particle end point.
     * @details The end point is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the x-coordinate of the particle end point.
     */
    template<class T>
    double end_x(const T & p)
    {
        return std::isinf(p.end_point[0]) ? PLACEHOLDERVALUE : (double)p.end_point[0];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, end_x, end_x);

    /**
     * @brief Variable for the y-coordinate of the particle end point.
     * @details The end point is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the y-coordinate of the particle end point.
     */
    template<class T>
    double end_y(const T & p)
    {
        return std::isinf(p.end_point[1]) ? PLACEHOLDERVALUE : (double)p.end_point[1];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, end_y, end_y);
    
    /**
     * @brief Variable for the z-coordinate of the particle end point.
     * @details The end point is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the z-coordinate of the particle end point.
     */
    template<class T>
    double end_z(const T & p)
    {
        return std::isinf(p.end_point[2]) ? PLACEHOLDERVALUE : (double)p.end_point[2];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, end_z, end_z);

    /**
     * @brief Variable for the x-component of the particle start direction.
     * @details The start direction is predicted upstream in the SPINE
     * reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the x-component of the particle start direction.
     */
    template<class T>
    double start_dir_x(const T & p)
    {
        return p.start_dir[0];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, start_dir_x, start_dir_x);
    
    /**
     * @brief Variable for the y-component of the particle start direction.
     * @details The start direction is predicted upstream in the SPINE
     * reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the y-component of the particle start direction.
     */
    template<class T>
    double start_dir_y(const T & p)
    {
        return p.start_dir[1];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, start_dir_y, start_dir_y);

    /**
     * @brief Variable for the z-component of the particle start direction.
     * @details The start direction is predicted upstream in the SPINE
     * reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the z-component of the particle start direction.
     */
    template<class T>
    double start_dir_z(const T & p)
    {
        return p.start_dir[2];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, start_dir_z, start_dir_z);

    /**
     * @brief Variable for the x-component of the particle end direction.
     * @details The end direction is predicted upstream in the SPINE
     * reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the x-component of the particle end direction.
     */
    template<class T>
    double end_dir_x(const T & p)
    {
        return p.end_dir[0];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, end_dir_x, end_dir_x);

    /**
     * @brief Variable for the y-component of the particle end direction.
     * @details The end direction is predicted upstream in the SPINE
     * reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the y-component of the particle end direction.
     */
    template<class T>
    double end_dir_y(const T & p)
    {
        return p.end_dir[1];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, end_dir_y, end_dir_y);

    /**
     * @brief Variable for the z-component of the particle end direction.
     * @details The end direction is predicted upstream in the SPINE
     * reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the z-component of the particle end direction.
     */
    template<class T>
    double end_dir_z(const T & p)
    {
        return p.end_dir[2];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, end_dir_z, end_dir_z);

    /**
     * @brief Variable for the x-component of the particle momentum.
     * @details The momentum is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the x-component of the particle momentum.
     */
    template<class T>
    double px(const T & p)
    {
        return p.momentum[0];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, px, px);
    
    /**
     * @brief Variable for the y-component of the particle momentum.
     * @details The momentum is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the x-component of the particle momentum.
     */
    template<class T>
    double py(const T & p)
    {
        return p.momentum[1];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, py, py);

    /**
     * @brief Variable for the z-component of the particle momentum.
     * @details The momentum is predicted upstream in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the z-component of the particle momentum.
     */
    template<class T>
    double pz(const T & p)
    {
        return p.momentum[2];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, pz, pz);
    
    /**
     * @brief Variable for the transverse momentum of a particle.
     * @details This function calculates the transverse momentum of the
     * particle with respect to the assumed neutrino direction. The neutrino
     * direction is assumed to either be the BNB axis direction (z-axis) or the
     * unit vector pointing from the NuMI target to the interaction vertex.
     * See @ref utilities::transverse_momentum for details on the extraction of
     * the transverse momentum.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the transverse momentum of the particle.
     */
    template<class T>
    double dpT(const T & p)
    {
        utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
        utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
        utilities::three_vector pt = utilities::transverse_momentum(momentum, vtx);
        return utilities::magnitude(pt);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, dpT, dpT);

    /**
     * @brief Variable for the polar angle (w.r.t the z-axis) of the particle.
     * @details The polar angle is defined as the arccosine of the z-component
     * of the momentum vector. This variable is useful for identifying particles
     * which are produced transversely to the beam.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the polar angle of the particle.
     */
    template<class T>
    double polar_angle(const T & p)
    {
        return std::acos(p.start_dir[2]);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, polar_angle, polar_angle);

    /**
     * @brief Variable for the azimuthal angle (w.r.t the z-axis) of the particle.
     * @details The azimuthal angle is defined as the arccosine of the x-component
     * of the momentum vector divided by the square root of the sum of the squares
     * of the x and y components of the momentum vector.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the azimuthal angle of the particle.
     */
    template<class T>
    double azimuthal_angle(const T & p)
    {
        return std::atan2(p.start_dir[1], p.start_dir[0]);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::BothParticle, azimuthal_angle, azimuthal_angle);

    /**
     * @brief Variable for the photon softmax score of the particle.
     * @details The photon softmax score represents the confidence that the
     * network has in the particle being a photon. The score is between 0 and 1,
     * with 1 being the most confident that the particle is a photon.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the photon softmax score of the particle.
     */
    template<class T>
    double photon_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[0];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, photon_softmax, photon_softmax);

    /**
     * @brief Variable for the electron softmax score of the particle.
     * @details The electron softmax score represents the confidence that the
     * network has in the particle being an electron. The score is between 0 and 1,
     * with 1 being the most confident that the particle is an electron.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the electron softmax score of the particle.
     */
    template<class T>
    double electron_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[1];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, electron_softmax, electron_softmax);
    
    /**
     * @brief Variable for the muon softmax score of the particle.
     * @details The muon softmax score represents the confidence that the
     * network has in the particle being a muon. The score is between 0 and 1,
     * with 1 being the most confident that the particle is a muon.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the muon softmax score of the particle.
     */
    template<class T>
    double muon_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[2];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, muon_softmax, muon_softmax);

    /**
     * @brief Variable for the pion softmax score of the particle.
     * @details The pion softmax score represents the confidence that the
     * network has in the particle being a pion. The score is between 0 and 1,
     * with 1 being the most confident that the particle is a pion.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the pion softmax score of the particle.
     */
    template<class T>
    double pion_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[3];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, pion_softmax, pion_softmax);

    /**
     * @brief Variable for the proton softmax score of the particle.
     * @details The proton softmax score represents the confidence that the
     * network has in the particle being a proton. The score is between 0 and 1,
     * with 1 being the most confident that the particle is a proton.
     * @tparam p the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the proton softmax score of the particle.
     */
    template<class T>
    double proton_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[4];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, proton_softmax, proton_softmax);

    /**
     * @brief Variable for the "MIP" softmax score of the particle.
     * @details The "MIP" softmax score is calculated as the sum of the softmax
     * scores for the muon and pion. The score represents the confidence that
     * the network has in the particle being a minimum ionizing particle.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the "MIP" softmax score of the particle.
     */
    template<class T>
    double mip_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[2] + p.pid_scores[3];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, mip_softmax, mip_softmax);

    /**
     * @brief Variable for the "hadron" softmax score of the particle.
     * @details The "hadron" softmax score is calculated as the sum of the softmax
     * scores for the pion and proton. The score represents the confidence that
     * the network has in the particle being a hadron.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the "hadron" softmax score of the particle.
     */
    template<class T>
    double hadron_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.pid_scores[3] + p.pid_scores[4];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, hadron_softmax, hadron_softmax);

    /**
     * @brief Variable for the primary softmax score of the particle.
     * @details The primary softmax score represents the confidence that the
     * network has in the particle being a primary particle. The score is between
     * 0 and 1, with 1 being the most confident that the particle is a primary
     * particle.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the primary softmax score of the particle.
     */
    template<class T>
    double primary_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.primary_scores[1];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, primary_softmax, primary_softmax);

    /**
     * @brief Variable for the secondary softmax score of the particle.
     * @details The secondary softmax score represents the confidence that the
     * network has in the particle being a secondary particle. The score is between
     * 0 and 1, with 1 being the most confident that the particle is a secondary
     * particle.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the secondary softmax score of the particle.
     */
    template<class T>
    double secondary_softmax(const caf::SRParticleDLPProxy & p)
    {
        return p.primary_scores[0];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, secondary_softmax, secondary_softmax);
}
#endif // PARTICLE_VARIABLES_H