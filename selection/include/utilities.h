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
                    ++counts[PIDFUNC(p)];
            }
            return counts;
        }
}
#endif // UTILITIES_H