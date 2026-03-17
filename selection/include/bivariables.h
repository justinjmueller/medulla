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
            a.start_dir[0] * b.start_dir[0] +
            a.start_dir[1] * b.start_dir[1] +
            a.start_dir[2] * b.start_dir[2]
        );
    }
    REGISTER_BIVAR_SCOPE(RegistrationScope::BothParticle, opening_angle, opening_angle);

    /**
     * @brief Invariant mass of two particles.
     * @details Computes the invariant mass of two particles using their
     * energy and momentum components.
     * @tparam T the type of particle (true or reco).
     * @param a the first particle.
     * @param b the second particle.
     * @return the invariant mass in MeV.
     */
    template<class T>
    double invariant_mass(const T & a, const T & b)
    {
        double ea = pvars::energy(a), eb = pvars::energy(b);
        double px = pvars::px(a) + pvars::px(b);
        double py = pvars::py(a) + pvars::py(b);
        double pz = pvars::pz(a) + pvars::pz(b);
        double e  = ea + eb;
        double m2 = e*e - px*px - py*py - pz*pz;
        return (m2 > 0) ? std::sqrt(m2) : 0.0;
    }
    REGISTER_BIVAR_SCOPE(RegistrationScope::BothParticle, invariant_mass, invariant_mass);

    /**
     * @brief Sum of kinetic energies of two particles.
     * @tparam T the type of particle (true or reco).
     * @param a the first particle.
     * @param b the second particle.
     * @return the sum of kinetic energies.
     */
    template<class T>
    double ke_sum(const T & a, const T & b)
    {
        return pvars::ke(a) + pvars::ke(b);
    }
    REGISTER_BIVAR_SCOPE(RegistrationScope::BothParticle, ke_sum, ke_sum);

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
}
#endif // BIVARIABLES_H
