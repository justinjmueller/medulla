/**
 * @file utilities.h
 * @brief Header file for definitions of utility functions acting on
 * interactions.
 * @details This file contains definitions of utility functions which are used
 * to support the implementation of analysis variables and cuts. These functions
 * are intended to be used to simplify the implementation of variables and cuts
 * by providing common functionality which can be reused across multiple
 * variables and cuts.
 * @author mueller@fnal.gov
 */
#ifndef UTILITIES_H
#define UTILITIES_H
#include <vector>
#include <array>

#include "framework.h"
#include "include/particle_variables.h"
#include "include/particle_cuts.h"

/**
 * @namespace utilities
 * @brief Namespace for organizing utility functions for supporting analysis
 * variables and cuts.
 * @details This namespace is intended to be used for organizing utility
 * functions which are used to support the implementation of analysis variables
 * and cuts. These functions are intended to be used to simplify the
 * implementation of variables and cuts by providing common functionality which
 * can be reused across multiple variables and cuts.
 * @note The namespace is intended to be used in conjunction with the
 * vars and cuts namespaces, which are used for organizing variables and cuts
 * which act on interactions.
 */
namespace utilities
{
    /**
     * @brief Count the primaries of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to find the topology of.
     * @return the count of primaries of each particle type within the
     * interaction.
     */
    template<class T>
        std::vector<uint32_t> count_primaries(const T & obj)
        {
            std::vector<uint32_t> counts(5, 0);
            for(auto &p : obj.particles)
            {
                if(pcuts::final_state_signal(p))
                    ++counts[pvars::pid(p)];
            }
            return counts;
        }

    /**
     * @brief List of good runs for ICARUS Run 2.
     * @details This list contains the run numbers that are considered good
     * for analysis in ICARUS Run 2. The list is used to filter out runs that
     * are known to have issues or are not suitable for analysis.
     * @note This is intended to be removed once we have a proper version-
     * controlled run list committed to sbnana.
     */
    constexpr std::array<unsigned int, 229>
    icarus_good_runs_run2{  9301,  9302,  9303,  9307,  9308,  9309,  9310,  9311,  9312,  9313,
                            9314,  9316,  9317,  9318,  9327,  9328,  9329,  9330,  9332,  9333,
                            9335,  9337,  9338,  9339,  9340,  9341,  9342,  9343,  9344,  9346,
                            9347,  9353,  9354,  9356,  9357,  9358,  9359,  9360,  9361,  9362,
                            9363,  9364,  9365,  9366,  9367,  9380,  9383,  9384,  9385,  9386,
                            9387,  9388,  9389,  9390,  9391,  9392,  9393,  9394,  9409,  9412,
                            9415,/*9435,*/9436,  9437,  9438,  9439,  9441,  9445,  9448,  9450,
                            9458,  9460,  9472,  9473,  9474,  9477,  9478,  9481,  9482,  9499,
                            9504,  9518,  9563,  9565,  9569,  9570,  9583,  9587,  9588,  9589,
                            9590,  9593,  9594,  9595,  9597,  9598,  9599,  9602,  9610,  9626,
                            9627,  9631,  9647,  9648,  9649,  9658,  9672,  9675,  9688,  9689,
                            9690,  9691,  9692,  9693,  9694,  9695,  9696,  9697,  9698,  9699,
                            9700,  9703,  9704,  9705,  9714,  9715,  9716,  9717,  9721,  9723,
                            9724,  9725,  9726,  9728,  9729,  9730,  9731,  9732,  9733,  9734,
                            9735,  9743,  9744,  9745,  9746,  9747,  9750,  9752,  9753,  9755,
                            9758,  9762,  9763,  9764,  9765,  9783,  9788,  9791,  9792,  9793,
                            9794,  9795,  9796,  9807,  9834,  9835,  9837,  9838,  9840,  9841,
                            9844,  9847,  9849,  9851,  9854,  9855,  9860,  9862,  9868,  9870,
                            9892,  9894,  9896,  9897,  9914,  9917,  9919,  9921,  9922,  9924,
                            9925,  9926,  9929,  9932,  9935,  9940,  9941,  9942,  9944,  9945,
                            9946,  9949,  9950,  9951,  9953,  9954,  9956,  9959,  9960,  9961,
                            9970,  9971,  9974,  9977,  9979,  9981,  9982,  9986, 10054, 10059,
                           10061, 10062, 10064, 10065, 10066, 10067, 10084, 10085, 10096, 10097};
    
    /**
     * @brief Boolean function to check if a run is in the list of good runs
     * for ICARUS Run 2.
     * @param run the run number to check.
     * @return true if the run is in the list of good runs, false otherwise.
     */
    bool is_icarus_good_run(unsigned int run)
    {
        return std::find(icarus_good_runs_run2.begin(), icarus_good_runs_run2.end(), run) != icarus_good_runs_run2.end();
    }

    /**
     * @brief Utility function to find the index of the first optical flash
     * @details This function returns the index of the first optical flash in
     * the event, shifted by a given amount. It is intended to be used to find
     * the flash nearest to the trigger of the event. The optional shift is
     * used to allow for tuning to correct for the natural offset of the
     * reconstructed flash time from the trigger time, which is not zero
     * despite all systems being referenced to the trigger. This function uses
     * the `firsttime` field of the optical flash.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @param shift the amount to shift the time of the flash by, default is 0.
     * @return size_t the index of the first optical flash in the event, or
     * kNoMatch if no flash is found.
     */
    template<typename T>
    size_t first_opflash_firsttime(const T & sr, double shift=0.0)
    {
        // This variable returns the index of the first optical flash in the
        // event, shifted by a given amount.
        if(sr.opflashes.empty()) return kNoMatch;
        size_t first_flash_index = 0;
        double min_time_diff = std::abs(sr.opflashes[0].firsttime - shift);
        for(size_t i = 1; i < sr.opflashes.size(); ++i)
        {
            double time_diff = std::abs(sr.opflashes[i].firsttime - shift);
            if(time_diff < min_time_diff)
            {
                min_time_diff = time_diff;
                first_flash_index = i;
            }
        }
        return first_flash_index;
    }

    /**
     * @brief Utility function to find the index of the first optical flash
     * @details This function returns the index of the first optical flash in
     * the event, shifted by a given amount. It is intended to be used to find
     * the flash nearest to the trigger of the event. The optional shift is
     * used to allow for tuning to correct for the natural offset of the
     * reconstructed flash time from the trigger time, which is not zero
     * despite all systems being referenced to the trigger. This function uses
     * the `time` field of the optical flash.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @param shift the amount to shift the time of the flash by, default is 0.
     * @return size_t the index of the first optical flash in the event, or
     * kNoMatch if no flash is found.
     */
    template<typename T>
    size_t first_opflash_rawtime(const T & sr, double shift=0.0)
    {
        // This variable returns the index of the first optical flash in the
        // event, shifted by a given amount.
        if(sr.opflashes.empty()) return kNoMatch;
        size_t first_flash_index = 0;
        double min_time_diff = std::abs(sr.opflashes[0].time - shift);
        for(size_t i = 1; i < sr.opflashes.size(); ++i)
        {
            double time_diff = std::abs(sr.opflashes[i].time - shift);
            if(time_diff < min_time_diff)
            {
                min_time_diff = time_diff;
                first_flash_index = i;
            }
        }
        return first_flash_index;
    }
}
#endif // UTILITIES_H
