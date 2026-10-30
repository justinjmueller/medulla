/**
 * @file bdt_pid.h
 * @brief Particle identification using the binary XGBoost BDTs.
 * @details This file defines PID functions which can be configured as the
 * global PID function (`pidfn = "bdt_binary_pid"` or
 * `"bdt_binary_boosted_pid"` in the [general] block of the configuration
 * file), so that every PID-dependent cut, selector, and variable defers to
 * the BDTs. Three single-switch BDTs are used: shower_ge (photon vs.
 * electron) for showers (shape 0), and track_muon and track_proton for tracks
 * (shape 1); a track takes the species with the largest of P(mu) (2 * P(mu)
 * for the boosted function), (1 - P(mu)) * (1 - P(p)) and P(p). The BDT
 * scores are also available as particle variables. The models are read from
 * `general.bdt_model_dir` (see export_bdt_cpp.py).
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
        static const std::array<std::string, kNumBDTModels> names = {"shower_ge", "track_proton", "track_muon"};
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
     * @return the model output (a single probability for the binary models).
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
     * @brief Variable for the binary-BDT PID of a shower (photon / electron).
     * @details Uses the shower_ge model (trained on true photons and
     * electrons only): P(photon) > 0.5 gives a photon, otherwise an electron.
     * @tparam T the type of particle (reco).
     * @param p the particle to apply the variable on.
     * @return the PID of the particle according to the binary shower BDT.
     */
    template<class T>
    double bdt_binary_shower_pid(const T & p)
    {
        return bdt_scores(p, kShowerGE)[0] > 0.5 ? kPhoton : kElectron;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_binary_shower_pid, bdt_binary_shower_pid);

    /**
     * @brief Track PID from the binary BDTs as the largest of three values.
     * @details With P(mu) from track_muon and P(p) from track_proton (there
     * is no pion model), the three values are muon_boost * P(mu),
     * (1 - P(mu)) * (1 - P(p)) (neither switch) and P(p); the track takes the
     * species with the largest value. Ties go to muon, then proton, then pion.
     * Not registered; used by @ref bdt_binary_track_pid and
     * @ref bdt_binary_boosted_track_pid.
     * @tparam T the type of particle (reco).
     * @param p the particle to classify.
     * @param muon_boost the factor applied to the muon value.
     * @return the PID of the particle.
     */
    template<class T>
    double bdt_binary_argmax_track_pid(const T & p, double muon_boost)
    {
        const double p_muon(bdt_scores(p, kTrackMuon)[0]);
        const double p_proton(bdt_scores(p, kTrackProton)[0]);
        const double muon(muon_boost * p_muon);
        const double pion((1 - p_muon) * (1 - p_proton));
        if(muon >= pion && muon >= p_proton) return kMuon;
        if(p_proton >= pion) return kProton;
        return kPion;
    }

    /**
     * @brief Variable for the binary-BDT PID of a track (muon / proton / pion).
     * @details The largest of P(mu), (1 - P(mu)) * (1 - P(p)) and P(p), see
     * @ref bdt_binary_argmax_track_pid.
     * @tparam T the type of particle (reco).
     * @param p the particle to apply the variable on.
     * @return the PID of the particle according to the binary track BDTs.
     */
    template<class T>
    double bdt_binary_track_pid(const T & p)
    {
        return bdt_binary_argmax_track_pid(p, 1.0);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_binary_track_pid, bdt_binary_track_pid);

    /**
     * @brief Variable for assigning PID with the binary BDTs.
     * @details Primary showers use @ref bdt_binary_shower_pid and primary
     * tracks @ref bdt_binary_track_pid. The BDTs were trained only on
     * primaries (primary_scores[1] > 0.10) of these two semantic types, so
     * other particles (non-primaries per the configured primfn, Michel, delta,
     * low energy) keep the SPINE PID with the boosted muon score
     * (@ref boosted_muon_pid). Skipping non-primaries also avoids evaluating
     * the BDTs on the many particles every PID cut ignores.
     * @tparam T the type of particle (reco).
     * @param p the particle to apply the variable on.
     * @return the PID of the particle.
     */
    template<class T>
    double bdt_binary_pid(const T & p)
    {
        if(!(*primfn)(p)) return boosted_muon_pid(p);
        if(p.shape == 0) return bdt_binary_shower_pid(p);
        if(p.shape == 1) return bdt_binary_track_pid(p);
        return boosted_muon_pid(p);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_binary_pid, bdt_binary_pid);

    /**
     * @brief Variable for the binary-BDT PID of a track with the muon score boosted.
     * @details As @ref bdt_binary_track_pid, but the muon value is doubled
     * (cf. @ref boosted_muon_pid): the largest of 2 * P(mu),
     * (1 - P(mu)) * (1 - P(p)) and P(p), see @ref bdt_binary_argmax_track_pid.
     * @tparam T the type of particle (reco).
     * @param p the particle to apply the variable on.
     * @return the PID of the particle according to the boosted binary track BDTs.
     */
    template<class T>
    double bdt_binary_boosted_track_pid(const T & p)
    {
        static constexpr double kMuonBoost = 2.0;
        return bdt_binary_argmax_track_pid(p, kMuonBoost);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_binary_boosted_track_pid, bdt_binary_boosted_track_pid);

    /**
     * @brief Variable for assigning PID with the binary BDTs, muon score boosted.
     * @details As @ref bdt_binary_pid, but primary tracks use
     * @ref bdt_binary_boosted_track_pid.
     * @tparam T the type of particle (reco).
     * @param p the particle to apply the variable on.
     * @return the PID of the particle.
     */
    template<class T>
    double bdt_binary_boosted_pid(const T & p)
    {
        if(!(*primfn)(p)) return boosted_muon_pid(p);
        if(p.shape == 0) return bdt_binary_shower_pid(p);
        if(p.shape == 1) return bdt_binary_boosted_track_pid(p);
        return boosted_muon_pid(p);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_binary_boosted_pid, bdt_binary_boosted_pid);

    /**
     * @brief Output of a binary BDT model.
     * @tparam T the type of particle (reco).
     * @param p the particle to evaluate.
     * @param m the model to evaluate.
     * @param shape the semantic type the model applies to.
     * @return the model output, or PLACEHOLDERVALUE if the particle has the
     * wrong semantic type for the model.
     */
    template<class T>
    double bdt_binary_score(const T & p, BDTModel_t m, int shape)
    {
        if(p.shape != shape)
            return PLACEHOLDERVALUE;
        return bdt_scores(p, m)[0];
    }

    /// @brief Binary shower BDT score: P(photon) vs. electron.
    template<class T> double bdt_photon_vs_electron_score(const T & p) { return bdt_binary_score(p, kShowerGE, 0); }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_photon_vs_electron_score, bdt_photon_vs_electron_score);

    /// @brief Binary track BDT score: P(proton).
    template<class T> double bdt_proton_binary_score(const T & p) { return bdt_binary_score(p, kTrackProton, 1); }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_proton_binary_score, bdt_proton_binary_score);

    /// @brief Binary track BDT score: P(muon).
    template<class T> double bdt_muon_binary_score(const T & p) { return bdt_binary_score(p, kTrackMuon, 1); }
    REGISTER_VAR_SCOPE(RegistrationScope::RecoParticle, bdt_muon_binary_score, bdt_muon_binary_score);
} // namespace pvars
#endif // BDT_PID_H
