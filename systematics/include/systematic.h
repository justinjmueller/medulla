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
#include "configuration.h"

#include "TTree.h"

namespace sys
{
    enum class Type { kMULTISIM, kMULTISIGMA, kVARIATION };

    class Systematic
    {
    private:
        std::string name;
        size_t index;
        Type type;
        std::string ordinate;
        std::vector<std::string> points;
        std::vector<double> scale;
        std::vector<double> nsigma;
        TTree * tree;
        std::vector<double> * weights;
        std::vector<double> * zscores;
        
    public:
        Systematic(sys::cfg::ConfigurationTable & table, TTree * t);
        size_t get_index();
        Type get_type();
        TTree * get_tree();
        std::vector<double> * & get_weights();
    };
} // namespace sys
#endif // SYSTEMATIC_H