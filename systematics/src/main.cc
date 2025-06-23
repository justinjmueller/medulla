/**
 * @file main.cc
 * @brief Main function for the code that adds TTrees with systematics to the
 * output ROOT file from CAFAna analyses.
 * @details This code is designed to add TTrees with systematics to the output
 * ROOT file. The code reads in the input ROOT file, which contains the TTrees
 * produced by the CAFAna analysis framework with the sBruce tree format. The 
 * code then matches the selected interactions in the sBruce trees with the 
 * corresponding universe weights from the original CAF files. The final output
 * is a ROOT file with TTrees that contain the selected interactions and some
 * additional TTrees containing the universe weights for the configured
 * systematics.
 * @author mueller@fnal.gov
 */
#include <iostream>
#include <string>

#include "configuration.h"
#include "trees.h"
#include "detsys.h"

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TInterpreter.h"

int main(int argc, char * argv[])
{
    /**
     * @brief Ignore ROOT warnings.
     * @details This block ignores ROOT warnings. Sometimes there are
     * mismatches in the dictionary for the StandardRecord object, which
     * causes ROOT to print warnings. These warnings can be safely ignored.
     */
    gErrorIgnoreLevel = kError;

    /**
     * @brief Check the number of arguments. The code expects the configuration
     * file as the only argument.
     * @details This block checks the number of arguments. The code expects the
     * configuration file as the only argument. If the number of arguments is not
     * correct, the code prints the usage and exits with an error code.
     */
    if(argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <configuration.toml>" << std::endl;
        return 1;
    }

    /**
     * @brief Read the configuration file.
     * @details This block reads the configuration file by using the
     * @ref cfg::Configuration class. The configuration
     * file is read using the toml++ library. If the configuration file is not
     * found or if there is an error in parsing the configuration file, the
     * code prints an error message and exits with an error code. Each requisite
     * field in the configuration file is checked for validity by the 
     * @ref cfg::validate() function
     * @see cfg::Configuration
     * @see cfg::ConfigurationError
     * @see cfg::validate()
     */
    cfg::ConfigurationTable config;
    try
    {
        config.set_config(argv[1]);
    }
    catch(const cfg::ConfigurationError & e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    /**
     * @brief Open the input and output ROOT files.
     * @details This block opens the input and output ROOT files. The input
     * ROOT file is the file that contains the TTrees produced by the CAFAna
     * analysis framework. The output ROOT file is the file that will contain
     * the TTrees that are produced by this code.
     */
    TFile * input = TFile::Open(config.get_string_field("input.path").c_str(), "READ");
    TFile * output = TFile::Open(config.get_string_field("output.path").c_str(), "RECREATE");

    /**
     * @brief Load the DetsysCalculator, if configured.
     * @details This block loads the DetsysCalculator if it is configured in
     * the configuration file. The DetsysCalculator is used to calculate the
     * detector systematics weights using a spline interpolation of the ratio
     * of the nominal and sample histograms.
     * @see sys::detsys::DetsysCalculator
     */
    sys::detsys::DetsysCalculator calc;
    if(config.has_field("variations"))
    {
        calc = sys::detsys::DetsysCalculator(config, output, input);
        calc.write();
    }

    /**
     * @brief Begin main loop over trees in the configuration file.
     * @details This block begins the main loop over the trees in the
     * configuration file. Each tree is a sub-table in the configuration file,
     * and can be fetched as a vector of sub-tables using the function
     * @ref cfg::get_subtables(). The main body of the loop then reads
     * delegates the handling of the tree to the appropriate function.
     * @see cfg::get_subtables()
     * @see cfg::ConfigurationTable
     */
    std::vector<cfg::ConfigurationTable> tables;
    try
    {
       tables = config.get_subtables("tree");
    }
    catch(const cfg::ConfigurationError & e)
    {
        /**
         * @TODO reconsider if this should cause the code to exit.
         */
        std::cout << "No trees found in the configuration file." << std::endl;
    }

    for(cfg::ConfigurationTable & table : tables)
    {
        std::cout << "Processing tree: " << table.get_string_field("origin") << std::endl;
        std::string type(table.get_string_field("action"));
        if(type == "copy")
            sys::trees::copy_tree(table, output, input);
        else if(type == "add_weights")
            sys::trees::copy_with_weight_systematics(config, table, output, input, calc);
    }

    input->Close();
    output->Close();

    return 0;
}