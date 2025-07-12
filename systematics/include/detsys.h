/**
 * @file detsys.h
 * @brief Header file for the DetsysCalculator class.
 * @details This file contains the header for the DetsysCalculator class. The
 * DetsysCalculator class is used to calculate the weights for the detector
 * systematics using a spline interpolation of the ratio of the nominal and
 * sample histograms. Configuration of the class is done using a TOML-based
 * configuration file.
 * @see sys::cfg::ConfigurationTable
 * @author mueller@fnal.gov
 */
#ifndef DETSYS_H
#define DETSYS_H
#include <map>
#include <vector>

#include "configuration.h"
#include "utilities.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TSpline.h"
#include "TFile.h"

/**
 * @namespace sys::detsys
 * @brief Namespace for functions that calculate the detector systematics
 * weights.
 * @details This namespace contains functionality for calculating the detector
 * systematics weights. The main class in this namespace that provides this
 * functionality is the DetsysCalculator class.
 */
namespace sys::detsys
{
    /**
     * @class DetsysCalculator
     * @brief Class for calculating the detector systematics weights.
     * @details This class is used to calculate the weights for the detector
     * systematics using a spline interpolation of the ratio of the nominal and
     * sample histograms. The class is configured using a TOML-based configuration
     * file. There is a main configuration block "variations" that specifies the
     * list of all variations to be used in the construction of the splines. Each
     * detector systematic result is separately configured in the list of subtables
     * "detsys". The configuration for each detector systematic includes the name
     * of the systematic, the list of points to be used in the spline construction,
     * and the scale factors for each point.
     * @see sys::cfg::ConfigurationTable
     */
    class DetsysCalculator
    {
    private:
        bool initialized;
        TDirectory * histogram_directory;
        TDirectory * result_directory;
        std::map<std::string, TH1D *> histograms;
        std::map<std::string, std::vector<double>> zscores;
        std::map<std::string, std::vector<TSpline3 *>> splines;
        std::map<std::string, TH1D *> hdummies;
        std::string variable;
        std::vector<std::string> variations;
        std::map<std::string, TH1D *> detsys_results1D;
        std::map<std::string, TH2D *> detsys_results2D;
        size_t nuniverses;
        double nominal_count;
        std::vector<double> random_zscores;

    public:
        /**
         * @brief Constructor for the DetsysCalculator class.
         * @details This constructor initializes the DetsysCalculator class
         * using the configuration table, the output file, and the input file.
         * The configuration table is used to configure the class. The output
         * file is used to write the histograms and splines. The input file is
         * used to read the histograms that define the variations.
         * @param table The configuration table.
         * @param output The output file.
         * @param input The input file.
         * @see sys::cfg::ConfigurationTable
         */
        DetsysCalculator(sys::cfg::ConfigurationTable & table, TFile * output, TFile * input);

        /**
         * @brief Default constructor for the DetsysCalculator class.
         * @details This constructor initializes the DetsysCalculator class when
         * no arguments are provided. The class is not initialized and is marked
         * as such.
         */
        DetsysCalculator();

        /**
         * @brief Add a variable to the list of result histograms.
         * @details This function adds a variable to the list of result
         * histograms by creating a new TH1D and TH2D for each configured
         * detector variation. The name of each histogram follows the naming
         * scheme "<variable>_<detsys>_1D" and "<variable>_<detsys>_2D".
         * @param variable The SysVariable object to be added to the list of
         * result histograms.
         * @return void
         */
        void add_variable(SysVariable & variable);

        /**
         * @brief Accessor method for the initialized flag.
         * @details This function returns the initialized flag.
         * @return The initialized flag.
         */
        bool is_initialized();

        /**
         * @brief Accessor method for the variable name.
         * @details This function returns the variable name.
         * @return The variable name.
         */
        std::string get_variable();

        /**
         * @brief Get the number of universes.
         * @details This function returns the number of universes.
         * @return The number of universes.
         */
        size_t get_nuniverses();

        /**
         * @brief Increment the nominal count by the specified value.
         * @details This function increments the nominal count by the specified
         * value.
         * @param value The value by which the nominal count is to be incremented.
         * @return void
         */
        void increment_nominal_count(double value);

        /**
         * @brief Get the z-scores for a specified detector systematic.
         * @details This function returns the z-scores for a specified detector
         * systematic.
         * @param name The name of the detector systematic.
         * @return The z-scores for the specified detector systematic.
         */
        std::vector<double> get_zscores(std::string name);

        /**
         * @brief Write the variation histograms to the output file.
         * @details This function writes the histograms to the output file. The
         * directory used for writing the histograms is set by the configuration
         * parameter "variations.histogram_destination" and the name of each
         * histogram is set by the configuration parameter "variations.keys".
         * @return void
         */
        void write();

        /**
         * @brief Write the result detector systematic variation histograms to
         * the output file.
         * @details This function writes the result detector systematic variation
         * histograms to the output file. The directory used for writing the
         * histograms is set by the configuration parameter
         * "variations.result_destination". The name of each histogram is set by
         * the configuration parameter "variations.keys". The histogram consists
         * of the binned variable (X) and the weighted entries across the
         * universes (Y).
         * @return void
         */
        void write_results();

        /**
         * @brief Access the histograms by name.
         * @details This function allows the histograms to be accessed by name.
         * @param key The name of the histogram to be accessed.
         * @return A pointer to the histogram with the specified name.
         */
        TH1D * operator[](std::string key);

        /**
         * @brief Get the weight for a given value and z-score for a specified
         * detector systematic.
         * @details This function calculates the weight for a given value and
         * z-score for a specified detector systematic. The weight is calculated
         * using a spline interpolation of the ratio of the ordinate and
         * variation histograms.
         * @param name The name of the detector systematic.
         * @param value The value for which the weight is to be calculated.
         * @param zscore The z-score for which the weight is to be calculated.
         * @return The weight for the specified value and z-score.
         */
        double get_weight(std::string name, double value, double zscore);

        /**
         * @brief Add a value to the histogram for a specified detector
         * systematic.
         * @details This function adds a value to the histogram for a specified
         * detector systematic while respecting the pre-rolled z-scores of
         * the random universes.
         * @param varname The name of the variable.
         * @param binvar The value of the binning variable.
         * @param detsysname The name of the detector systematic.
         * @param value The value to be added to the histogram.
         * @return void
         */
        void add_value(std::string varname, double binvar, std::string detsysname, double value);
    };
} // namespace sys::detsys
#endif // DETSYS_H