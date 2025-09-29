/**
 * @file scorers.h
 * @brief Header file for functions which are intended to be used by the user
 * at the configuration file level to define PID / scoring / enumeration of
 * categories.
 * @details This file contains functions which are intended to be used by the
 * user at the configuration file level to define PID / scoring / enumeration
 * of categories. Though these may be implemented as functions in the @ref vars
 * or @ref pvars namespace, they are moved here to provide a cleaner interface.
 * @author mueller@fnal.gov
 */
#ifndef SCORERS_H
#define SCORERS_H

namespace pvars
{
    /**
     * @brief Enumeration for particle primary classification.
     * @details This is the order that the particle ids occur in (most) vectors in the framework,
     * with a unknown option to handle exceptions.
     */
    enum Particle_t
    {
        kPhoton   =  0,
        kElectron =  1,
        kMuon     =  2,
        kPion     =  3,
        kProton   =  4,
        kUnknown  = -1
    };

    /**
     * @brief Function pointer for the primary classification function.
     * @details This function pointer is used to allow the user to configure
     * the primary classification function used in the analysis. The primary
     * classification function is used to assign the primary classification of
     * the particle based on the softmax scores of the particle. This method
     * specifically allows the user to override the default primary
     * classification function with a custom one, which can be set in the
     * configuration file.
     */
    extern std::shared_ptr<VarFn<RParticleType>> primfn;

    /**
     * @brief Function pointer for the PID function.
     * @details This function pointer is used to allow the user to configure
     * the PID function used in the analysis. The PID function is used to
     * assign the PID of the particle based on the softmax scores of the
     * particle. This method specifically allows the user to override the
     * default PID function with a custom one, which can be set in the
     * configuration file.
     */
    extern std::shared_ptr<VarFn<RParticleType>> pidfn;

    /**
     * @brief Variable for the particle's primary classification.
     * @details This variable returns the primary classification of the
     * particle. The primary classification is determined upstream in the SPINE
     * reconstruction and is based on the softmax scores of the particle. This
     * function uses the "nominal" primary classification that is made upstream
     * in the SPINE reconstruction.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the primary classification of the particle.
     */
    template<class T>
    double default_primary_classification(const T & p)
    {
        return p.is_primary ? 1 : 0;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, default_primary_classification, default_primary_classification);

    /**
     * @brief Variable for assigning primary classification based on the
     * particle's softmax scores.
     * @details This variable assigns a primary classification based on the
     * softmax scores of the particle. This function places a relaxed threshold
     * on the primary softmax score to reduce observed inefficiencies in the
     * primary classification.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to apply the variable on.
     * @return the primary classification of the particle.
     */
    template<class T>
    double lax_primary_classification(const T & p)
    {
        return p.primary_scores[1] > 0.10 ? 1 : 0;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, lax_primary_classification, lax_primary_classification);

    /**
     * @brief Variable for the particle's PID.
     * @details This variable returns the PID of the particle. The PID is
     * determined by the softmax scores of the particle. This function uses the
     * "nominal" PID decision that is made upstream in the SPINE reconstruction.
     * @param p the particle to apply the variable on.
     * @return the PID of the particle.
     */
    template<class T>
    double default_pid(const T & p)
    {
        return p.pid;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, default_pid, default_pid);

    /**
     * @brief Variable for assigning PID based on the particle's softmax
     * scores.
     * @details This variable assigns a PID based on the softmax scores of the
     * particle. Nominally, the PID is assigned based on the highest softmax
     * score, but the PID can be overridden directly by this function to
     * effectively loosen the requirements to assign a muon type PID to a
     * particle.
     * @param p the particle to apply the variable on.
     * @return the PID of the particle.
     */
    template<class T>
    double lax_muon_pid(const T & p)
    {
        double pid = std::numeric_limits<double>::quiet_NaN();
        if(p.pid_scores[pvars::kMuon] > 0.25)
            pid = 2;
        else
        {
            size_t high_index(0);
            for(size_t i(0); i < 5; ++i)
                if(p.pid_scores[i] > p.pid_scores[high_index]) high_index = i;
            pid = high_index;
        }
        return pid;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, lax_muon_pid, lax_muon_pid);
} // namespace pvars
#endif