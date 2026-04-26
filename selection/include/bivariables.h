/**
 * @file bivariables.h
 * @brief Header file for definitions of bivariables which act on pairs of
 * particles.
 * @details This file contains definitions of bivariables which act on pairs
 * of particles. Each bivariable is implemented as a function which takes two
 * particle objects as arguments and returns a double. These variables are
 * intended to be used to define observables that depend on the relationship
 * between two particles in an interaction.
 * @author mueller@fnal.gov
 */
#ifndef BIVARIABLES_H
#define BIVARIABLES_H
#include <cmath>
#include <array>

#include "include/particle_variables.h"
#include "framework.h"

/**
 * @namespace bvars
 * @brief Namespace for organizing variables which act on pairs of particles.
 * @details This namespace is intended to be used for organizing variables which
 * act on pairs of particles. Each variable is implemented as a function which
 * takes two particle objects as arguments and returns a double. The function
 * should be templated on the type of particle object if the variable is
 * intended to be used on both true and reconstructed particles.
 */
namespace bvars
{
    /**
     * @brief Opening angle between two particles.
     * @details Computes the arccosine of the dot product of the two
     * particles' start direction vectors.
     * @tparam T the type of particle (true or reco).
     * @param a the first particle.
     * @param b the second particle.
     * @return the opening angle in radians.
     */
    template<class T>
    double opening_angle(const T & a, const T & b)
    {
        return std::acos(
            pvars::start_dir_x(a) * pvars::start_dir_x(b) +
            pvars::start_dir_y(a) * pvars::start_dir_y(b) +
            pvars::start_dir_z(a) * pvars::start_dir_z(b)
        );
    }
    REGISTER_BIVAR_SCOPE(RegistrationScope::BothParticle, opening_angle, opening_angle);

    /**
     * @brief Distance between start points of two particles.
     * @details Computes the Euclidean distance between the start points
     * of two particles.
     * @tparam T the type of particle (true or reco).
     * @param a the first particle.
     * @param b the second particle.
     * @return the Euclidean distance between start points.
     */
    template<class T>
    double start_distance(const T & a, const T & b)
    {
        double dx = pvars::start_x(a) - pvars::start_x(b);
        double dy = pvars::start_y(a) - pvars::start_y(b);
        double dz = pvars::start_z(a) - pvars::start_z(b);
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
    REGISTER_BIVAR_SCOPE(RegistrationScope::BothParticle, start_distance, start_distance);

    // ---- Internal pi0 helpers (non-standard signatures, not registered) ----
    // All KE values are in MeV, distances in cm, momenta in MeV/c.

    /**
     * @brief Kinetic energy of a single pi0 photon daughter in MeV.
     * @details Uses p.ke for true particles (SRParticleTruthDLP stores KE in MeV)
     * and p.calo_ke for reco particles (also MeV).
     */
    template<class T>
    double pi0_shower_ke(const T & p)
    {
        if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
            return p.ke;
        else
            return p.calo_ke;
    }

    /**
     * @brief Distance from the interaction vertex to the shower start point (cm).
     */
    template<class T>
    double pi0_conv_dist(const T & p, double vx, double vy, double vz)
    {
        double dx = p.start_point[0] - vx;
        double dy = p.start_point[1] - vy;
        double dz = p.start_point[2] - vz;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    /**
     * @brief Cosine of the opening angle between two pi0 photon daughters.
     * @details Reco: dot product of vertex-to-shower-start unit vectors.
     * True: dot product of momentum unit vectors.
     */
    template<class T>
    double pi0_opening_costheta_impl(const T & p0, const T & p1,
                                     double vx, double vy, double vz)
    {
        if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
        {
            double px0 = p0.momentum[0], py0 = p0.momentum[1], pz0 = p0.momentum[2];
            double px1 = p1.momentum[0], py1 = p1.momentum[1], pz1 = p1.momentum[2];
            double r0 = std::sqrt(px0*px0 + py0*py0 + pz0*pz0);
            double r1 = std::sqrt(px1*px1 + py1*py1 + pz1*pz1);
            return (px0*px1 + py0*py1 + pz0*pz1) / (r0 * r1);
        }
        else
        {
            double dx0 = p0.start_point[0]-vx, dy0 = p0.start_point[1]-vy, dz0 = p0.start_point[2]-vz;
            double dx1 = p1.start_point[0]-vx, dy1 = p1.start_point[1]-vy, dz1 = p1.start_point[2]-vz;
            double r0 = std::sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);
            double r1 = std::sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
            return (dx0*dx1 + dy0*dy1 + dz0*dz1) / (r0 * r1);
        }
    }

    template<class T>
    double pi0_leading_dirx(const T & p0, const T & /*p1*/)
    {
        if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
            return p0.momentum[0] / std::sqrt(p0.momentum[0]*p0.momentum[0] + p0.momentum[1]*p0.momentum[1] + p0.momentum[2]*p0.momentum[2]);
        else
        {
            double dx = p0.start_point[0]-context::current<T>()->vertex[0];
            double dy = p0.start_point[1]-context::current<T>()->vertex[1];
            double dz = p0.start_point[2]-context::current<T>()->vertex[2];
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            return dx / r;
        }
    }
    REGISTER_BIVAR_SCOPE(RegistrationScope::BothParticle, pi0_leading_dirx, pi0_leading_dirx);

    template<class T>
    double pi0_leading_diry(const T & p0, const T & /*p1*/)
    {
        if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
            return p0.momentum[1] / std::sqrt(p0.momentum[0]*p0.momentum[0] + p0.momentum[1]*p0.momentum[1] + p0.momentum[2]*p0.momentum[2]);
        else
        {
            double dx = p0.start_point[0]-context::current<T>()->vertex[0];
            double dy = p0.start_point[1]-context::current<T>()->vertex[1];
            double dz = p0.start_point[2]-context::current<T>()->vertex[2];
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            return dy / r;
        }
    }
    REGISTER_BIVAR_SCOPE(RegistrationScope::BothParticle, pi0_leading_diry, pi0_leading_diry);

    template<class T>
    double pi0_leading_dirz(const T & p0, const T & /*p1*/)
    {
        if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
            return p0.momentum[2] / std::sqrt(p0.momentum[0]*p0.momentum[0] + p0.momentum[1]*p0.momentum[1] + p0.momentum[2]*p0.momentum[2]);
        else
        {
            double dx = p0.start_point[0]-context::current<T>()->vertex[0];
            double dy = p0.start_point[1]-context::current<T>()->vertex[1];
            double dz = p0.start_point[2]-context::current<T>()->vertex[2];
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            return dz / r;
        }
    }
    REGISTER_BIVAR_SCOPE(RegistrationScope::BothParticle, pi0_leading_dirz, pi0_leading_dirz);

    /**
     * @brief Pi0 3-momentum vector in MeV/c.
     * @details Reco: calo_ke (MeV) * vertex-to-start unit vector, summed over
     * both photons. True: sum of the two daughter momentum vectors
     * (SRParticleTruthDLP stores momentum in MeV/c).
     */
    template<class T>
    std::array<double, 3> pi0_momentum_vec_impl(const T & p0, const T & p1,
                                                double vx, double vy, double vz)
    {
        if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
        {
            return { p0.momentum[0] + p1.momentum[0],
                     p0.momentum[1] + p1.momentum[1],
                     p0.momentum[2] + p1.momentum[2] };
        }
        else
        {
            double ke0 = p0.calo_ke, ke1 = p1.calo_ke;
            double dx0 = p0.start_point[0]-vx, dy0 = p0.start_point[1]-vy, dz0 = p0.start_point[2]-vz;
            double dx1 = p1.start_point[0]-vx, dy1 = p1.start_point[1]-vy, dz1 = p1.start_point[2]-vz;
            double r0 = std::sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);
            double r1 = std::sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
            return { ke0*dx0/r0 + ke1*dx1/r1,
                     ke0*dy0/r0 + ke1*dy1/r1,
                     ke0*dz0/r0 + ke1*dz1/r1 };
        }
    }

    // ---- Registered pi0 kinematic bivariables ----
    // All use the standard (p0, p1) -> double signature and read the parent
    // interaction via context::current<T>() when vertex coordinates are needed.
    // The framework sets context::current_{true,reco} in the bivar dispatch
    // lambda before invoking any bivariable, so the pointer is always valid.

    /**
     * @brief Leading photon KE in GeV: pi0_shower_ke(p0) / 1000.
     */
    template<class T>
    double pi0_leading_ke(const T & p0, const T & /*p1*/)
    {
        return pi0_shower_ke(p0) / 1000.0;
    }
    REGISTER_BIVAR_SCOPE(RegistrationScope::BothParticle, pi0_leading_ke, pi0_leading_ke);

    /**
     * @brief Subleading photon KE in GeV: pi0_shower_ke(p1) / 1000.
     */
    template<class T>
    double pi0_subleading_ke(const T & /*p0*/, const T & p1)
    {
        return pi0_shower_ke(p1) / 1000.0;
    }
    REGISTER_BIVAR_SCOPE(RegistrationScope::BothParticle, pi0_subleading_ke, pi0_subleading_ke);

    /**
     * @brief Total pi0 energy in GeV: (ke0 + ke1) / 1000.
     */
    template<class T>
    double pi0_energy(const T & p0, const T & p1)
    {
        return (pi0_shower_ke(p0) + pi0_shower_ke(p1)) / 1000.0;
    }
    REGISTER_BIVAR_SCOPE(RegistrationScope::BothParticle, pi0_energy, pi0_energy);

    /**
     * @brief Distance from vertex to leading shower start point (cm).
     */
    template<class T>
    double pi0_leading_conv_dist(const T & p0, const T & /*p1*/)
    {
        auto * obj = context::current<T>();
        return pi0_conv_dist(p0, obj->vertex[0], obj->vertex[1], obj->vertex[2]);
    }
    REGISTER_BIVAR_SCOPE(RegistrationScope::BothParticle, pi0_leading_conv_dist, pi0_leading_conv_dist);

    /**
     * @brief Distance from vertex to subleading shower start point (cm).
     */
    template<class T>
    double pi0_subleading_conv_dist(const T & /*p0*/, const T & p1)
    {
        auto * obj = context::current<T>();
        return pi0_conv_dist(p1, obj->vertex[0], obj->vertex[1], obj->vertex[2]);
    }
    REGISTER_BIVAR_SCOPE(RegistrationScope::BothParticle, pi0_subleading_conv_dist, pi0_subleading_conv_dist);

    /**
     * @brief Cosine of the diphoton opening angle.
     * @details Reco: vertex-to-shower-start dot product. True: momentum unit vectors.
     */
    template<class T>
    double pi0_opening_costheta(const T & p0, const T & p1)
    {
        auto * obj = context::current<T>();
        return pi0_opening_costheta_impl(p0, p1, obj->vertex[0], obj->vertex[1], obj->vertex[2]);
    }
    REGISTER_BIVAR_SCOPE(RegistrationScope::BothParticle, pi0_opening_costheta, pi0_opening_costheta);

    /**
     * @brief Diphoton invariant mass in MeV: sqrt(2 ke0 ke1 (1 - cos theta)).
     */
    template<class T>
    double pi0_invariant_mass(const T & p0, const T & p1)
    {
        auto * obj = context::current<T>();
        double vx = obj->vertex[0], vy = obj->vertex[1], vz = obj->vertex[2];
        double ke0 = pi0_shower_ke(p0), ke1 = pi0_shower_ke(p1);
        double ct  = pi0_opening_costheta_impl(p0, p1, vx, vy, vz);
        return std::sqrt(2.0 * ke0 * ke1 * (1.0 - ct));
    }
    REGISTER_BIVAR_SCOPE(RegistrationScope::BothParticle, pi0_invariant_mass, pi0_invariant_mass);

    /**
     * @brief Maximum allowed diphoton mass in MeV: (ke0+ke1) sqrt((1-cos theta)/2).
     */
    template<class T>
    double pi0_mass_max(const T & p0, const T & p1)
    {
        auto * obj = context::current<T>();
        double vx = obj->vertex[0], vy = obj->vertex[1], vz = obj->vertex[2];
        double ke0 = pi0_shower_ke(p0), ke1 = pi0_shower_ke(p1);
        double ct  = pi0_opening_costheta_impl(p0, p1, vx, vy, vz);
        return (ke0 + ke1) * std::sqrt((1.0 - ct) / 2.0);
    }
    REGISTER_BIVAR_SCOPE(RegistrationScope::BothParticle, pi0_mass_max, pi0_mass_max);

    /**
     * @brief Pi0 momentum magnitude in GeV.
     */
    template<class T>
    double pi0_momentum_mag(const T & p0, const T & p1)
    {
        auto * obj = context::current<T>();
        auto [mpx, mpy, mpz] = pi0_momentum_vec_impl(p0, p1, obj->vertex[0], obj->vertex[1], obj->vertex[2]);
        return std::sqrt(mpx*mpx + mpy*mpy + mpz*mpz) / 1000.0;
    }
    REGISTER_BIVAR_SCOPE(RegistrationScope::BothParticle, pi0_momentum_mag, pi0_momentum_mag);

    /**
     * @brief Cosine of the angle between the pi0 momentum and the beam axis (z).
     */
    template<class T>
    double pi0_beam_costheta(const T & p0, const T & p1)
    {
        auto * obj = context::current<T>();
        auto [mpx, mpy, mpz] = pi0_momentum_vec_impl(p0, p1, obj->vertex[0], obj->vertex[1], obj->vertex[2]);
        double mag = std::sqrt(mpx*mpx + mpy*mpy + mpz*mpz);
        return mpz / mag;
    }
    REGISTER_BIVAR_SCOPE(RegistrationScope::BothParticle, pi0_beam_costheta, pi0_beam_costheta);
}
#endif // BIVARIABLES_H
