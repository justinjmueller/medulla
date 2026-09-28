/**
 * @file trees.h
 * @brief Header file for the trees namespace.
 * @details This file contains the header for the trees namespace. The trees
 * namespace contains functions that read and interface with the TTrees
 * produced by the CAFAna analysis framework. Different "copying" actions can
 * be performed on the TTrees, such as a simple copy or adding systematics to
 * the output file based on the selected signal candidates and the configured
 * systematics.
 * @author mueller@fnal.gov
 */
#ifndef TREES_H
#define TREES_H
#include <iostream>

#include "detsys.h"
#include "configuration.h"

#include "TFile.h"

/**
 * @namespace sys::trees
 * @brief Namespace for functions that read and interface with the TTrees
 * produced by the CAFAna analysis framework.
 * @details This namespace contains functions that read and interface with the
 * TTrees produced by the CAFAna analysis framework. Different "copying" 
 * actions can be performed on the TTrees, such as a simple copy or the
 * addition of systematics to the output file based on the selected signal
 * candidates and the configured systematics.
 */
namespace sys::trees
{
    /**
     * @brief Type definitions for systematic indexing (variable name and
     * index).
     */
    typedef std::pair<std::string, int64_t> syst_t;
    typedef std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, double> index_t;

    /**
     * @brief Copy the input TTree to the output TTree.
     * @details This function copies the input TTree to the output TTree. The
     * function loops over the input TTree and copies the values of the branches
     * to the output TTree. The output TTree is created with the same branches
     * as the input TTree.
     * @param table The table that contains the configuration for the tree.
     * @param output The output TFile.
     * @param input The input TFile.
     * @return void
     */
    void copy_tree(cfg::ConfigurationTable & table, TFile * output, TFile * input);

    /**
     * @brief Add reweightable systematics to the output TTree.
     * @details This function adds reweightable systematics to the output
     * TTree. The function loops over the input TTree to build a map for the
     * selected signal candidates to their index in the input TTree. The
     * function then loops over the neutrinos in the CAF input files and
     * populates the output TTree with the selected signal candidates and the
     * universe weights for matched neutrinos.
     * @param table The table that contains the configuration for the tree.
     * @param output The output TFile.
     * @param input The input TFile.
     * @return void
     */
    void copy_with_weight_systematics(cfg::ConfigurationTable & config, cfg::ConfigurationTable & table, TFile * output, TFile * input, sys::detsys::DetsysCalculator & calc, sys::detsys::DetsysCalculator * calc_nue = nullptr, const std::vector<int> & nue_categories = {});

    /**
     * @brief Apply pre-loaded detector variation weights without reading CAF files.
     * @details Phase-2 counterpart to copy_with_weight_systematics. Loops
     * directly over the input selection tree and evaluates pre-built splines
     * from a DetsysCalculator loaded via the spline-loading constructor. No
     * WeightReader is needed because detector variation weights depend only on
     * the reconstructed variable value and the pre-rolled z-scores, both of
     * which are available without re-reading the original CAF files.
     *
     * If @p calc_nue is non-null and @p nue_categories is non-empty, events
     * whose true_category falls in @p nue_categories use @p calc_nue to
     * evaluate splines; all other events use @p calc. This allows nue-enhanced
     * splines (built from a nue-enriched sample) to be applied to nue-flavor
     * events while background events receive weights from the standard splines.
     *
     * @param config The global configuration table.
     * @param table The per-tree configuration sub-table.
     * @param output The output TFile.
     * @param input The input TFile (individual job selection output).
     * @param calc A DetsysCalculator initialised with standard pre-built splines.
     * @param calc_nue Optional DetsysCalculator with nue-enhanced splines.
     *                 Pass nullptr to use @p calc for all events (default behaviour).
     * @param nue_categories true_category values that should use @p calc_nue.
     * @return void
     */
    void copy_with_detsys_weights(cfg::ConfigurationTable & config, cfg::ConfigurationTable & table, TFile * output, TFile * input, sys::detsys::DetsysCalculator & calc, sys::detsys::DetsysCalculator * calc_nue = nullptr, const std::vector<int> & nue_categories = {});
}
#endif