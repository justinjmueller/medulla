/**
 * @file systematic.cc
 * @brief Source file for the Systematic class.
 * @details This file contains the source for the Systematic class. The
 * Systematic class is used to encapsulate the different types of
 * systematics that can be applied to the analysis.
 * @author mueller@fnal.gov
 */
#include "systematic.h"
#include "configuration.h"

#include "TTree.h"

// Construct a new Systematic object.
sys::Systematic::Systematic(cfg::ConfigurationTable & table, TTree * t)
    : name(table.get_string_field("name")),
      index(table.get_int_field("index")),
      type(table.get_string_field("type") == "multisim" ? Type::kMULTISIM : table.get_string_field("type") == "multisigma" ? Type::kMULTISIGMA : Type::kVARIATION),
      tree(t),
      weights(new std::vector<double>()),
      nsigma(new std::vector<double>())
{
    if(table.has_field("nsigma"))
        *nsigma = table.get_double_vector("nsigma");

    if(type == Type::kVARIATION)
    {
        ordinate = table.get_string_field("ordinate");
        points = table.get_string_vector("points");
    }
    else
    {
        if(table.has_field("scale"))
            scale = table.get_double_vector("scale");
        else
            scale = std::vector<double>(points.size(), 1);
    }
}

// Get the index of the systematic parameter.
size_t sys::Systematic::get_index()
{
    return index;
}

// Get the type of the systematic.
sys::Type sys::Systematic::get_type()
{
    return type;
}

// Get the TTree associated with the systematic.
TTree * sys::Systematic::get_tree()
{
    return tree;
}

// Get a reference to the weights vector.
std::vector<double> * & sys::Systematic::get_weights()
{
    return weights;
}

// Get a reference to the nsigma vector.
std::vector<double> * & sys::Systematic::get_nsigma()
{
    return nsigma;
}