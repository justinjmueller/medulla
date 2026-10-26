/**
 * @file bdt_pid.h
 * @brief Particle identification using the XGBoost BDTs.
 * @details This file defines PID functions which can be configured as the
 * global PID function (`pidfn = "bdt_pid"` in the [general] block of the
 * configuration file), so that every PID-dependent cut, selector, and
 * variable defers to the BDTs. Showers (shape 0) are classified by the
 * shower_multi model (photon / electron / other) and tracks (shape 1) by the
 * track_multi model (proton / muon / pion / other). Particles classified as
 * "other" are assigned @ref pvars::kOther, which matches no particle species.
 * The per-class BDT scores are also available as particle variables. The
 * models are read from `general.bdt_model_dir` (see export_bdt_cpp.py).
 * @author bmweiss
 */
#ifndef BDT_PID_H
#define BDT_PID_H

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

#include "framework.h"
#include "scorers.h"
#include "bdt_inference.h"

namespace pvars
{
    /**
     * @brief Directory containing the exported BDT models.
     * @details Set from `general.bdt_model_dir` in the configuration file.
     */
    inline std::string bdt_model_dir = "/exp/sbnd/app/users/bmweiss/workdir/pi0_analysis/BDTs/models/cpp";

    /**
     * @brief Enumeration of the available BDT models.
     */
    enum BDTModel_t
    {
        kShowerMulti,   ///< photon / electron / other
        kTrackMulti,    ///< proton / muon / pion / other
        kShowerGE,      ///< binary: P(photon) vs. electron
        kTrackProton,   ///< binary: P(proton)
        kTrackMuon,     ///< binary: P(muon)
        kNumBDTModels
    };

    /**
     * @brief Retrieve a BDT model, loading it on first use.
     * @param m the model to retrieve.
     * @return the loaded model.
     */
    inline const bdt::XGBForest & bdt_model(BDTModel_t m)
    {
        static std::array<bdt::XGBForest, kNumBDTModels> models;
        static const std::array<std::string, kNumBDTModels> names =
            {"shower_multi", "track_multi", "shower_ge", "track_proton", "track_muon"};
        if(!models[m].loaded())
            models[m].load(bdt_model_dir + "/" + names[m] + ".txt");
        return models[m];
    }

    /**
     * @brief Collect the raw BDT inputs from a reconstructed particle.
     * @tparam T the type of particle (reco).
     * @param p the particle to read.
     * @return the raw inputs to the BDT feature builder.
     */
    template<class T>
    bdt::ParticleInputs bdt_inputs(const T & p)
    {
        bdt::ParticleInputs in;
        for(size_t i(0); i < 6; ++i)
            in.pid_scores[i] = p.pid_scores[i];
        // Only the muon, pion, and proton hypotheses are features.
        in.csda_ke_per_pid.fill(0);
        in.mcs_ke_per_pid.fill(0);
        for(size_t i : {kMuon, kPion, kProton})
        {
            in.csda_ke_per_pid[i] = p.csda_ke_per_pid[i];
            in.mcs_ke_per_pid[i] = p.mcs_ke_per_pid[i];
        }
        in.primary_score = p.primary_scores[1];
        in.ke = p.ke;
        in.calo_ke = p.calo_ke;
        in.csda_ke = p.csda_ke;
        in.mcs_ke = p.mcs_ke;
        in.depositions_sum = p.depositions_sum;
        in.size = p.size;
        in.length = p.length;
        in.start_dedx = p.start_dedx;
        in.is_contained = p.is_contained;
        in.start_dir_z = p.start_dir[2];
        in.vertex_distance = p.vertex_distance;
        in.start_straightness = p.start_straightness;
        in.directional_spread = p.directional_spread;
        in.axial_spread = p.axial_spread;
        return in;
    }

    /**
     * @brief Evaluate a BDT model on a particle.
     * @details The PID function is called many times per particle by the
     * cuts and selectors, so results are cached. The cache key is a handful
     * of fields that identify the particle (cheaper to read than the full
     * feature vector); the features are only built on a cache miss.
     * @tparam T the type of particle (reco).
     * @param p the particle to evaluate.
     * @param m the model to evaluate.
     * @return the class probabilities (one entry for a binary model).
     */
    template<class T>
    std::vector<double> bdt_scores(const T & p, BDTModel_t m)
    {
        static std::array<std::unordered_map<std::string, std::vector<double>>, kNumBDTModels> cache;
        const std::array<double, 8> id = {(double)p.id, (double)p.interaction_id, (double)p.size,
                                          (double)p.depositions_sum, (double)p.length, (double)p.start_point[0],
                                          (double)p.start_point[1], (double)p.start_point[2]};
        std::string key(reinterpret_cast<const char *>(id.data()), sizeof(id));
        auto & c = cache[m];
        auto it = c.find(key);
        if(it != c.end())
            return it->second;
        if(c.size() > 100000)
            c.clear();
        return c.emplace(std::move(key), bdt_model(m).predict(bdt::build_features(bdt_inputs(p)))).first->second;
    }

    /**
     * @brief Map a BDT class name to the framework particle enumeration.
     * @param name the class name stored with the model.
     * @return the corresponding Particle_t (kOther if not a particle species).
     */
    inline double bdt_class_to_pid(const std::string & name)
    {
        if(name == "photon")   return kPhoton;
        if(name == "electron") return kElectron;
        if(name == "muon")     return kMuon;
        if(name == "pion")     return kPion;
        if(name == "proton")   return kProton;
        return kOther;
    }

    /**
     * @brief PID from the argmax of a multi-class BDT model.
     * @tparam T the type of particle (reco).
     * @param p the particle to classify.
     * @param m the multi-class model to use.
     * @return the PID of the particle (kOther for the "other" class).
     */
    template<class T>
    double bdt_argmax_pid(const T & p, BDTModel_t m)
    {
        const std::vector<double> scores = bdt_scores(p, m);
        size_t high_index = std::max_element(scores.begin(), scores.end()) - scores.begin();
        return bdt_class_to_pid(bdt_model(m).classes()[high_index]);
    }

    /**
     * @brief Variable for the BDT PID of a shower (photon / electron / other).
     * @tparam T the type of particle (reco).
     * @param p the particle to apply the variable on.
     * @return the PID of the particle according to the shower BDT.
     */
    template<class T>
    double bdt_shower_pid(const T & p)
    {
        return bdt_argmax_pid(p, kShowerMulti);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_shower_pid, bdt_shower_pid);

    /**
     * @brief Variable for the BDT PID of a track (proton / muon / pion / other).
     * @tparam T the type of particle (reco).
     * @param p the particle to apply the variable on.
     * @return the PID of the particle according to the track BDT.
     */
    template<class T>
    double bdt_track_pid(const T & p)
    {
        return bdt_argmax_pid(p, kTrackMulti);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_track_pid, bdt_track_pid);

    /**
     * @brief Variable for assigning PID with the BDTs.
     * @details Primary showers use @ref bdt_shower_pid and primary tracks
     * use @ref bdt_track_pid. The BDTs were trained only on primaries
     * (primary_scores[1] > 0.10) of these two semantic types, so other
     * particles (non-primaries per the configured primfn, Michel, delta, low
     * energy) keep the SPINE PID with the boosted muon score
     * (@ref boosted_muon_pid). Skipping non-primaries also avoids evaluating
     * the BDTs on the many particles every PID cut ignores.
     * @tparam T the type of particle (reco).
     * @param p the particle to apply the variable on.
     * @return the PID of the particle.
     */
    template<class T>
    double bdt_pid(const T & p)
    {
        if(!(*primfn)(p)) return boosted_muon_pid(p);
        if(p.shape == 0) return bdt_shower_pid(p);
        if(p.shape == 1) return bdt_track_pid(p);
        return boosted_muon_pid(p);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_pid, bdt_pid);

    /**
     * @brief Score of one class of a BDT model.
     * @tparam T the type of particle (reco).
     * @param p the particle to evaluate.
     * @param m the model to evaluate.
     * @param shape the semantic type the model applies to.
     * @param cls the index of the class (ignored for binary models).
     * @return the class score, or PLACEHOLDERVALUE if the particle has the
     * wrong semantic type for the model.
     */
    template<class T>
    double bdt_class_score(const T & p, BDTModel_t m, int shape, size_t cls)
    {
        if(p.shape != shape)
            return PLACEHOLDERVALUE;
        const std::vector<double> scores = bdt_scores(p, m);
        return scores.size() == 1 ? scores[0] : scores[cls];
    }

    /// @brief Shower BDT photon score.
    template<class T> double bdt_photon_score(const T & p) { return bdt_class_score(p, kShowerMulti, 0, 0); }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_photon_score, bdt_photon_score);

    /// @brief Shower BDT electron score.
    template<class T> double bdt_electron_score(const T & p) { return bdt_class_score(p, kShowerMulti, 0, 1); }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_electron_score, bdt_electron_score);

    /// @brief Shower BDT "other" score.
    template<class T> double bdt_shower_other_score(const T & p) { return bdt_class_score(p, kShowerMulti, 0, 2); }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_shower_other_score, bdt_shower_other_score);

    /// @brief Track BDT proton score.
    template<class T> double bdt_proton_score(const T & p) { return bdt_class_score(p, kTrackMulti, 1, 0); }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_proton_score, bdt_proton_score);

    /// @brief Track BDT muon score.
    template<class T> double bdt_muon_score(const T & p) { return bdt_class_score(p, kTrackMulti, 1, 1); }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_muon_score, bdt_muon_score);

    /// @brief Track BDT pion score.
    template<class T> double bdt_pion_score(const T & p) { return bdt_class_score(p, kTrackMulti, 1, 2); }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_pion_score, bdt_pion_score);

    /// @brief Track BDT "other" score.
    template<class T> double bdt_track_other_score(const T & p) { return bdt_class_score(p, kTrackMulti, 1, 3); }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_track_other_score, bdt_track_other_score);

    /// @brief Binary shower BDT score: P(photon) vs. electron.
    template<class T> double bdt_photon_vs_electron_score(const T & p) { return bdt_class_score(p, kShowerGE, 0, 0); }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_photon_vs_electron_score, bdt_photon_vs_electron_score);

    /// @brief Binary track BDT score: P(proton).
    template<class T> double bdt_proton_binary_score(const T & p) { return bdt_class_score(p, kTrackProton, 1, 0); }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_proton_binary_score, bdt_proton_binary_score);

    /// @brief Binary track BDT score: P(muon).
    template<class T> double bdt_muon_binary_score(const T & p) { return bdt_class_score(p, kTrackMuon, 1, 0); }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_muon_binary_score, bdt_muon_binary_score);
} // namespace pvars
#endif // BDT_PID_H
