/**
 * @file bdt_inference.h
 * @brief Self-contained evaluator for the XGBoost particle-ID BDTs.
 * @details The BDTs are trained in python (pi0_analysis/BDTs) and exported
 * by export_bdt_cpp.py to a flat text format. This header reads that format
 * and evaluates the forest exactly as XGBoost does: float32 comparisons,
 * "x < threshold" goes left, missing (NaN) values follow the learned
 * default direction. It has no dependence on ROOT or the CAF classes so the
 * same code can be validated standalone against the python predictions.
 * @author bmweiss
 */
#ifndef BDT_INFERENCE_H
#define BDT_INFERENCE_H

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * @namespace bdt
 * @brief Namespace for the BDT particle-ID inference.
 */
namespace bdt
{
    /**
     * @brief Raw per-particle quantities needed to build the BDT features.
     * @details These mirror the rec.dlp.particles.* branches read by
     * extract_particles.py before any sanitization or derivation.
     */
    struct ParticleInputs
    {
        std::array<double, 6> pid_scores;
        double primary_score;
        double ke, calo_ke, csda_ke, mcs_ke;
        double depositions_sum;
        double size;
        double length;
        double start_dedx;
        double is_contained;
        double start_dir_z;
        double vertex_distance, start_straightness, directional_spread, axial_spread;
        std::array<double, 6> csda_ke_per_pid;
        std::array<double, 6> mcs_ke_per_pid;
    };

    /**
     * @brief Identifiers for every feature known to the builder.
     */
    enum Feature_t
    {
        kPidScore0, kPidScore1, kPidScore2, kPidScore3, kPidScore4, kPidScore5,
        kPrimaryScore, kKE, kCaloKE, kDepositionsSum, kSize, kLength, kStartDedx,
        kIsContained, kStartDirZ, kQPerVoxel,
        kVertexDistance, kStartStraightness, kDirectionalSpread, kAxialSpread,
        kCsdaKE2, kCsdaKE3, kCsdaKE4, kMcsKE2, kMcsKE3, kMcsKE4,
        kCsdaKE, kMcsKE, kMcsCsdaRatioMu, kMcsCsdaRatioP, kDqdxMean, kVoxelsPerCm,
        kNumFeatures
    };

    using Features = std::array<float, kNumFeatures>;

    /**
     * @brief Map from the python feature names to the feature identifiers.
     */
    inline const std::map<std::string, Feature_t> & feature_names()
    {
        static const std::map<std::string, Feature_t> names = {
            {"pid_scores_0", kPidScore0}, {"pid_scores_1", kPidScore1}, {"pid_scores_2", kPidScore2},
            {"pid_scores_3", kPidScore3}, {"pid_scores_4", kPidScore4}, {"pid_scores_5", kPidScore5},
            {"primary_score", kPrimaryScore}, {"ke", kKE}, {"calo_ke", kCaloKE},
            {"depositions_sum", kDepositionsSum}, {"size", kSize}, {"length", kLength},
            {"start_dedx", kStartDedx}, {"is_contained", kIsContained}, {"start_dir_z", kStartDirZ},
            {"q_per_voxel", kQPerVoxel}, {"vertex_distance", kVertexDistance},
            {"start_straightness", kStartStraightness}, {"directional_spread", kDirectionalSpread},
            {"axial_spread", kAxialSpread},
            {"csda_ke_per_pid_2", kCsdaKE2}, {"csda_ke_per_pid_3", kCsdaKE3}, {"csda_ke_per_pid_4", kCsdaKE4},
            {"mcs_ke_per_pid_2", kMcsKE2}, {"mcs_ke_per_pid_3", kMcsKE3}, {"mcs_ke_per_pid_4", kMcsKE4},
            {"csda_ke", kCsdaKE}, {"mcs_ke", kMcsKE}, {"mcs_csda_ratio_mu", kMcsCsdaRatioMu},
            {"mcs_csda_ratio_p", kMcsCsdaRatioP}, {"dqdx_mean", kDqdxMean}, {"voxels_per_cm", kVoxelsPerCm}
        };
        return names;
    }

    /**
     * @brief Build every feature from the raw inputs.
     * @details Reproduces extract_particles.py:sanitize() exactly: negative
     * sentinels -> NaN, MCS failures (> 5e4) -> NaN, non-positive CSDA KE per
     * PID -> NaN, then the derived ratios, then +-inf -> NaN. Arithmetic is
     * done in double (as pandas does) and cast to float32 (as XGBoost does).
     * @param in the raw particle inputs.
     * @return the full feature array indexed by Feature_t.
     */
    inline Features build_features(const ParticleInputs & in)
    {
        constexpr double nan = std::numeric_limits<double>::quiet_NaN();
        constexpr double kMcsFail = 5e4;
        auto neg_to_nan  = [&](double v) { return v < 0 ? nan : v; };
        auto mcs_to_nan  = [&](double v) { return v > kMcsFail ? nan : v; };
        auto csda_to_nan = [&](double v) { return v <= 0 ? nan : v; };

        std::array<double, kNumFeatures> f;
        for(size_t i(0); i < 6; ++i)
            f[kPidScore0 + i] = in.pid_scores[i];
        f[kPrimaryScore]      = in.primary_score;
        f[kKE]                = neg_to_nan(in.ke);
        f[kCaloKE]            = neg_to_nan(in.calo_ke);
        f[kDepositionsSum]    = in.depositions_sum;
        f[kSize]              = in.size;
        f[kLength]            = in.length;
        f[kStartDedx]         = neg_to_nan(in.start_dedx);
        f[kIsContained]       = in.is_contained;
        f[kStartDirZ]         = in.start_dir_z;
        f[kVertexDistance]    = neg_to_nan(in.vertex_distance);
        f[kStartStraightness] = neg_to_nan(in.start_straightness);
        f[kDirectionalSpread] = neg_to_nan(in.directional_spread);
        f[kAxialSpread]       = neg_to_nan(in.axial_spread);
        f[kCsdaKE2]           = csda_to_nan(in.csda_ke_per_pid[2]);
        f[kCsdaKE3]           = csda_to_nan(in.csda_ke_per_pid[3]);
        f[kCsdaKE4]           = csda_to_nan(in.csda_ke_per_pid[4]);
        f[kMcsKE2]            = mcs_to_nan(in.mcs_ke_per_pid[2]);
        f[kMcsKE3]            = mcs_to_nan(in.mcs_ke_per_pid[3]);
        f[kMcsKE4]            = mcs_to_nan(in.mcs_ke_per_pid[4]);
        f[kCsdaKE]            = neg_to_nan(in.csda_ke);
        f[kMcsKE]             = mcs_to_nan(neg_to_nan(in.mcs_ke));

        f[kQPerVoxel]         = f[kDepositionsSum] / f[kSize];
        f[kDqdxMean]          = f[kDepositionsSum] / f[kLength];
        f[kVoxelsPerCm]       = f[kSize] / f[kLength];
        f[kMcsCsdaRatioMu]    = f[kMcsKE2] / f[kCsdaKE2];
        f[kMcsCsdaRatioP]     = f[kMcsKE4] / f[kCsdaKE4];

        Features out;
        for(size_t i(0); i < kNumFeatures; ++i)
            out[i] = std::isinf(f[i]) ? std::numeric_limits<float>::quiet_NaN() : static_cast<float>(f[i]);
        return out;
    }

    /**
     * @brief A gradient-boosted forest exported from XGBoost.
     */
    class XGBForest
    {
    public:
        /**
         * @brief Load the forest from the text format written by export_bdt_cpp.py.
         * @param path the path to the model file.
         */
        void load(const std::string & path)
        {
            std::ifstream in(path);
            if(!in)
                throw std::runtime_error("bdt::XGBForest: cannot open model file " + path);

            std::string key;
            int version(0);
            in >> key >> version;
            if(key != "xgbforest" || version != 1)
                throw std::runtime_error("bdt::XGBForest: bad header in " + path);

            size_t n(0), num_trees(0);
            while(in >> key)
            {
                if(key == "task") in >> task_;
                else if(key == "objective") in >> objective_;
                else if(key == "num_class") in >> num_class_;
                else if(key == "base_margin") in >> base_margin_;
                else if(key == "features")
                {
                    in >> n;
                    feature_index_.resize(n);
                    for(size_t i(0); i < n; ++i)
                    {
                        std::string name;
                        in >> name;
                        auto it = feature_names().find(name);
                        if(it == feature_names().end())
                            throw std::runtime_error("bdt::XGBForest: unknown feature " + name + " in " + path);
                        feature_index_[i] = it->second;
                    }
                }
                else if(key == "classes")
                {
                    in >> n;
                    classes_.resize(n);
                    for(auto & c : classes_) in >> c;
                }
                else if(key == "num_trees")
                {
                    in >> num_trees;
                    break;
                }
                else throw std::runtime_error("bdt::XGBForest: unknown key " + key + " in " + path);
            }
            if(objective_ != "binary:logistic" && objective_ != "multi:softprob")
                throw std::runtime_error("bdt::XGBForest: unsupported objective " + objective_);

            trees_.resize(num_trees);
            for(auto & t : trees_)
            {
                size_t num_nodes(0);
                in >> key >> num_nodes >> t.cls;
                if(key != "tree" || !in || t.cls >= num_class_)
                    throw std::runtime_error("bdt::XGBForest: malformed tree block in " + path);
                t.nodes.resize(num_nodes);
                for(auto & node : t.nodes)
                {
                    double value;
                    int default_left;
                    in >> node.left >> node.right >> node.feature >> value >> default_left;
                    node.value = static_cast<float>(value);
                    node.default_left = default_left != 0;
                }
            }
            if(!in)
                throw std::runtime_error("bdt::XGBForest: truncated model file " + path);
            loaded_ = true;
        }

        /**
         * @brief Evaluate the forest on a full feature array.
         * @param features the feature array from @ref build_features.
         * @return class probabilities (one entry for a binary model: P(class 1)).
         */
        std::vector<double> predict(const Features & features) const
        {
            std::vector<float> x(feature_index_.size());
            for(size_t i(0); i < x.size(); ++i)
                x[i] = features[feature_index_[i]];

            std::vector<double> margin(num_class_, base_margin_);
            for(const auto & t : trees_)
            {
                int n(0);
                while(t.nodes[n].left != -1)
                {
                    const Node & node = t.nodes[n];
                    const float v = x[node.feature];
                    if(std::isnan(v)) n = node.default_left ? node.left : node.right;
                    else n = (v < node.value) ? node.left : node.right;
                }
                margin[t.cls] += t.nodes[n].value;
            }

            if(objective_ == "binary:logistic")
                return {1.0 / (1.0 + std::exp(-margin[0]))};

            double max_margin = *std::max_element(margin.begin(), margin.end());
            double sum(0);
            for(double & m : margin) { m = std::exp(m - max_margin); sum += m; }
            for(double & m : margin) m /= sum;
            return margin;
        }

        bool loaded() const { return loaded_; }
        const std::string & task() const { return task_; }
        const std::vector<std::string> & classes() const { return classes_; }

    private:
        struct Node
        {
            int left, right, feature;
            float value;  // split threshold, or leaf weight when left == -1
            bool default_left;
        };
        struct Tree
        {
            size_t cls;
            std::vector<Node> nodes;
        };

        bool loaded_ = false;
        std::string task_, objective_;
        size_t num_class_ = 1;
        double base_margin_ = 0;
        std::vector<Feature_t> feature_index_;
        std::vector<std::string> classes_;
        std::vector<Tree> trees_;
    };
} // namespace bdt
#endif // BDT_INFERENCE_H
