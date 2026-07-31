/**
 * @file systematic.h
 * @brief Header file for the Systematic class.
 * @details This file contains the header for the Systematic class. The
 * Systematic class is used to encapsulate the different types of
 * systematics that can be applied to the analysis.
 * @author mueller@fnal.gov
 */
#ifndef SYSTEMATIC_H
#define SYSTEMATIC_H
#include <algorithm>
#include <limits>

#include "configuration.h"

#include "TTree.h"

namespace sys
{
    /**
     * @brief Enumeration for the types of systematics.
     * @details The enumeration defines the types of systematics that can be
     * applied to the analysis. The types are:
     * - kMULTISIM: "universe-style" systematics that can be used directly to
     *   evaluate a covariance matrix.
     * - kMULTISIGMA: "spline-based" systematics that can be used to cast the
     *   systematic as a pull-term.
     * - kVARIATION: "variation-based" systematics, which reflect detector
     *   model variations.
     */
    enum class Type { kMULTISIM, kMULTISIGMA, kVARIATION };

    /**
     * @brief Class representing a systematic.
     * @details The Systematic class encapsulates the properties of a
     * systematic parameter such as its name, index, type, ordinate
     * (variation-specific), and other configuration information. It also
     * contains the TTrees associated with the systematic, which holds the
     * universe weights and z-scores for the systematic.
     */
    class Systematic
    {
    private:
        std::string name;
        size_t index;
        Type type;
        std::string ordinate;
        std::vector<std::string> points;
        std::vector<double> scale;
        TTree * tree;
        std::vector<double> * weights;
        std::vector<double> * nsigma;
        double clip_min;
        double clip_max;

    public:
        /**
         * @brief Construct a new Systematic object
         * @details This constructor initializes the Systematic object with the
         * configuration table and the TTree associated with the systematic.
         * @param table the configuration table containing the systematic
         * configuration details
         * @param t the TTree associated with the systematic
         * @param default_clip the [min, max] weight clip range to use if the
         * systematic's own configuration table does not provide a "clip"
         * field. An empty vector means no clipping is applied.
         */
        Systematic(cfg::ConfigurationTable & table, TTree * t,
                   std::vector<double> default_clip = {});

        /**
         * @brief Get the index of the systematic parameter.
         * @details This function returns the index of the systematic parameter
         * in the TTree. The index is used to identify the systematic in the
         * CAF input files.
         * @return size_t the index of the systematic parameter.
         */
        size_t get_index();

        /**
         * @brief Get the type of the systematic.
         * @details This function returns the type of the systematic, which can
         * be one of the types defined in the Type enumeration.
         * @return Type the type of the systematic.
         */
        Type get_type();

        /**
         * @brief Get the TTree associated with the systematic.
         * @details This function returns a pointer to the TTree associated
         * with the systematic. The TTree contains the universe weights and
         * z-scores for the systematic.
         * @return TTree* pointer to the TTree associated with the systematic.
         */
        TTree * get_tree();

        /**
         * @brief Get a reference to the weights vector.
         * @details This function returns a reference to the weights vector
         * associated with the systematic. The weights vector contains the
         * universe weights for the systematic.
         * @return std::vector<double>*& reference to the weights vector.
         */
        std::vector<double> * & get_weights();

        /**
         * @brief Get a reference to the nsigma vector.
         * @details This function returns a reference to the nsigma vector
         * associated with the systematic. The nsigma vector contains the
         * z-scores for the systematic.
         * @return std::vector<double>*& reference to the nsigma vector.
         */
        std::vector<double> * & get_nsigma();

        /**
         * @brief Clip a weight to the configured [min, max] range.
         * @details This function clamps the input weight to the range
         * resolved for this systematic (either its own "clip" field, the
         * general default passed in at construction, or, if neither is
         * present, no clipping at all).
         * @param weight the weight to clip.
         * @return double the clipped weight.
         */
        double clip(double weight) const;
    };
} // namespace sys
#endif // SYSTEMATIC_H