/**
 * @file utilities.h
 * @brief Header and implementation of functions generally useful for the
 * systematics framework within SPINE.
 * @details This file contains the implementation of functions that are
 * generally useful for the systematics framework within SPINE. These functions
 * are used to read in the selected signal candidates.
 * @author mueller@fnal.gov
 */
#ifndef UTILITIES_H
#define UTILITIES_H
#include <iostream>
#include <map>
#include <string>
#include <tuple>

#include "TTreeReader.h"
#include "TTreeReaderValue.h"

/**
 * @class SysVariable
 * @brief Simple struct to hold a variable's definition.
 * @details This struct holds a variable's definition, which includes the name
 * of the variable, the number of bins, the minimum value, and the maximum
 * value. This is used to contain the configuration of the variables that are
 * used in the systematics framework.
 * @see cfg::ConfigurationTable
 */
struct SysVariable
{
    SysVariable(cfg::ConfigurationTable & table)
        : name(table.get_string_field("name"))
    {
        std::vector<double> bins = table.get_double_vector("bins");
        nbins = bins[0];
        min = bins[1];
        max = bins[2];
    }
    std::string name;
    size_t nbins;
    double min;
    double max;
};

/**
 * @brief Create a nested directory structure in the output ROOT file.
 * @details This function creates a nested directory structure in the output
 * ROOT file. The directory structure is parsed from the desired output
 * destination in the configuration file.
 * @param parent The parent directory in the output ROOT file.
 * @param path The desired output destination.
 * @return directory The final directory in the output ROOT file.
 */
TDirectory * create_directory(TDirectory * parent, std::string path)
{
    TDirectory * directory = parent;
    size_t pos(path.find_last_of("/"));
    if(pos != std::string::npos)
    {
        std::string subdirs(path.substr(0, pos));
        std::string subdir;
        while((pos = subdirs.find_first_of("/")) != std::string::npos)
        {
            subdir = subdirs.substr(0, pos);
            if(directory->Get(subdir.c_str()) == 0)
                directory = directory->mkdir(subdir.c_str());
            else
                directory = directory->GetDirectory(subdir.c_str());
            subdirs = subdirs.substr(pos+1);
        }
        if(directory->Get(subdirs.c_str()) == 0)
            directory = directory->mkdir(subdirs.c_str());
        else
            directory = directory->GetDirectory(subdirs.c_str());
    }
    return directory;
}

/**
 * @brief Retrieve the parent directory of the object specified in the path.
 * @details This function retrieves the parent directory of the object specified
 * in the path by recursively crawling down the directory structure until the
 * parent directory is found.
 * @param parent The parent directory in the output ROOT file.
 * @param path The desired output destination.
 * @return directory The parent directory of the object specified in the path.
 */
TDirectory * get_parent_directory(TDirectory * parent, std::string path)
{
    TDirectory * directory = parent;
    size_t pos(path.find_last_of("/"));
    if(pos != std::string::npos)
    {
        std::string subdirs(path.substr(0, pos));
        std::string subdir;
        while((pos = subdirs.find_first_of("/")) != std::string::npos)
        {
            subdir = subdirs.substr(0, pos);
            directory = directory->GetDirectory(subdir.c_str());
            subdirs = subdirs.substr(pos+1);
        }
        directory = directory->GetDirectory(subdirs.c_str());
    }
    else
        directory = parent;
    return directory;
}

#endif // UTILITIES_H