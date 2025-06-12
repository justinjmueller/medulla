/**
 * @file analysis.h
 * @brief Header file for the Analysis class, a class designed to streamline
 * the running of multiple samples through CAFAna.
 * @details The Analysis class operates under the principle that the full
 * dataset in an analysis is comprised of multiple samples, each of which
 * represents a different component of the signal or background. Under this
 * paradigm, an analysis consists of running the same set of variables/cuts
 * on each sample, and then combining the results at the end. The Analysis
 * class is designed to facilitate this process by allowing the user to
 * specify the variables and cuts to be applied to each sample, along with the
 * full list of samples to be run. The class then handles the running of the
 * samples and the storing of the results.
 * @author mueller@fnal.gov
 */
#ifndef ANALYSIS_H
#define ANALYSIS_H
#include <vector>
#include <string>

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

#include "TDirectory.h"
#include "TFile.h"

/**
 * @namespace ana
 * @brief Namespace for the Analysis class and related functions.
 * @details The ana namespace contains the Analysis class, which is designed
 * to streamline the running of multiple samples through CAFAna. The namespace
 * is also used to organize analysis-related functions and variables within
 * CAFAna.
 */
namespace ana
{
    /**
     * @struct Sample
     * @brief Struct to store information about a sample in the analysis.
     * @details This struct is used to store information about a sample in the
     * analysis, including the name of the sample, the SpectrumLoader object
     * representing the sample, and a boolean indicating whether the sample is
     * a simulation sample. The simulation flag is used to determine if truth
     * information is available for the sample.
     */
    struct Sample
    {
        std::string name;
        ana::SpectrumLoader * loader;
        bool is_sim;
    };

    /**
     * @struct TreeSet
     * @brief Struct to store information about a set of variables that comprise
     * a Tree in the analysis.
     * @details This struct is used to store information about a set of variables
     * that comprise a Tree in the analysis. The Tree is used to store the results
     * of the analysis in a TTree in the output ROOT file. The struct contains the
     * name of the Tree, the names of the variables, the SpillMultiVars that
     * implement the variables, and a boolean indicating whether the Tree
     * represents a simulation sample. The simulation flag is used to determine
     * if truth information is present in the Tree.
     */
    struct TreeSet
    {
        std::string name;
        std::vector<std::string> names;
        std::vector<ana::SpillMultiVar> vars;
        bool is_sim;
    };

    /**
     * @class Analysis
     * @brief Class designed to streamline the running of multiple samples
     * through CAFAna.
     * @details The Analysis class operates under the principle that the full
     * dataset in an analysis is comprised of multiple samples, each of which
     * represents a different component of the signal or background. Under this
     * paradigm, an analysis consists of running the same set of variables/cuts
     * on each sample, and then combining the results at the end. The Analysis
     * class is designed to facilitate this process by allowing the user to
     * specify the variables and cuts to be applied to each sample, along with
     * the full list of samples to be run. The class then handles the running of
     * the samples and the storing of the results.
     */
    class Analysis
    {
        public:
            Analysis(std::string name);
            void AddLoader(std::string name, ana::SpectrumLoader * loader, bool is_sim);
            void AddTree(std::string name, std::map<std::string, ana::SpillMultiVar> & vars, bool is_sim);
            void AddTreeForSample(std::string sname, std::string name, std::map<std::string, ana::SpillMultiVar> & vars, bool is_sim);
            void Go();
        private:
            std::string name;
            std::vector<Sample> samples;
            std::vector<TreeSet> trees;
            std::map<std::pair<std::string, std::string>, TreeSet> trees_map;
    };

    /**
     * @brief Default constructor for the Analysis class.
     * @details This constructor initializes an instance of the Analysis class
     * with a default name. 
     * @param name The name of the Analysis class. This name is used in the
     * creation of the output ROOT file.
     * @return A new instance of the Analysis class.
     */
    Analysis::Analysis(std::string name)
    {
        this->name = name;
    }

    /**
     * @brief Add a SpectrumLoader to the Analysis class.
     * @details This function allows the user to add a SpectrumLoader to the
     * Analysis class. The SpectrumLoader represents a sample in the analysis,
     * and is used to load the data from the ROOT file and apply the cuts and
     * variables to the data.
     * @param name The name of the SpectrumLoader to be added to the Analysis
     * class.
     * @param loader The SpectrumLoader to be added to the Analysis class.
     * @param is_sim A boolean indicating whether the SpectrumLoader represents
     * a simulation sample, which is principally used to determine if truth
     * information is available.
     * @return void
     */
    void Analysis::AddLoader(std::string name, ana::SpectrumLoader * loader, bool is_sim)
    {
        samples.push_back({name, loader, is_sim});
    }

    /**
     * @brief Add a Tree to the Analysis class (set of variables and names).
     * @details This function allows the user to add a new Tree to the Analysis
     * class. A Tree is a set of variables with associated names that are
     * applied to the data in the SpectrumLoader. The Tree is used to store the
     * results of the analysis in a TTree in the output ROOT file.
     * @param names A vector of strings containing the names of the variables
     * to use in the Tree.
     * @param vars A vector of SpillMultiVar objects implementing the variables
     * to use in the Tree.
     * @param is_sim A boolean indicating whether the Tree represents a
     * simulation sample, which is principally used to determine if truth
     * information is available.
     * @return void
     */
    void Analysis::AddTree(std::string name, std::map<std::string, ana::SpillMultiVar> & vars, bool is_sim)
    {
        std::vector<std::string> n;
        std::vector<ana::SpillMultiVar> v;
        for(const auto & [name, var] : vars)
        {
            n.push_back(name);
            v.push_back(var);
        }
        trees.push_back({name, n, v, is_sim});
    }

    /**
     * @brief Add a Tree to the Analysis class for a specific sample.
     * @details This function allows the user to add a new Tree to the Analysis
     * class for a specific sample. A Tree is a set of variables with associated
     * names that are applied to the data in the SpectrumLoader. The Tree is used
     * to store the results of the analysis in a TTree in the output ROOT file.
     * @param sname The name of the sample to which the Tree belongs.
     * @param name The name of the Tree to be added to the Analysis class.
     * @param vars A map of variable names to SpillMultiVar objects implementing
     * the variables to use in the Tree.
     * @param is_sim A boolean indicating whether the Tree represents a simulation
     * sample, which is principally used to determine if truth information is
     * available.
     * @return void
     */
    void Analysis::AddTreeForSample(std::string sname, std::string name, std::map<std::string, ana::SpillMultiVar> & vars, bool is_sim)
    {
        std::vector<std::string> n;
        std::vector<ana::SpillMultiVar> v;
        for(const auto & [name, var] : vars)
        {
            n.push_back(name);
            v.push_back(var);
        }
        trees_map[std::make_pair(sname, name)] = {name, n, v, is_sim};
    }

    /**
     * @brief Run the analysis on the specified samples.
     * @details This function runs the analysis on the configured samples by
     * looping over each sample, creating the Trees for each sample, then
     * running the analysis on the sample to populate the Trees with the
     * results of the analysis. The results are stored in a TFile in the output
     * ROOT file in a parent directory named "events" and a subdirectory for
     * each sample.
     * @return void
     */
    void Analysis::Go()
    {
        TFile * f = new TFile(std::string(name + ".root").c_str(), "RECREATE");
        TDirectory * dir = f->mkdir("events");
        dir->cd();

        for(const Sample & s : samples)
        {
            TDirectory * subdir = dir->mkdir(s.name.c_str());
            subdir->cd();
            std::vector<ana::Tree*> sbruce_trees;
            for(const TreeSet & t : trees)
            {
                if(t.is_sim && !s.is_sim)
                    continue;
                sbruce_trees.push_back(new ana::Tree(t.name, t.names, *s.loader, t.vars, ana::kNoSpillCut, true));
            }
            for(const auto & [name, t] : trees_map)
            {
                if((t.is_sim && !s.is_sim) || name.first != s.name)
                    continue;
                sbruce_trees.push_back(new ana::Tree(t.name, t.names, *s.loader, t.vars, ana::kNoSpillCut, true));
            }

            s.loader->Go();
            for(const ana::Tree * t : sbruce_trees)
            {
                t->SaveTo(subdir);
                delete t;
            }
            dir->cd();
        }
        f->Close();
    }
}
#endif // ANALYSIS_H