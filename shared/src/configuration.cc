/**
 * @file configuration.cc
 * @brief Implementation of the ConfigurationTable class and related functions.
 * @details This file contains the implementation of the ConfigurationTable
 * class and related functions. The ConfigurationTable class is used to read
 * and interface with the TOML-based configuration file. The configuration file
 * consolidates all the parameters that may be analysis-specific or that may
 * change frequently. The main backend employed in this project is the toml++
 * library.
 * @author mueller@fnal.gov
 */
#include <string>

#include "configuration.h"
#include "toml++/toml.h"

namespace cfg
{
    // Set the configuration file by reading the TOML configuration table and
    // validating the configuration file.
    void ConfigurationTable::set_config(const std::string & path)
    {
        // Parse the TOML file.
        try
        {
            doc = toml::parse_file(path);
        }
        catch(const std::exception& e)
        {
            throw ConfigurationError("Failed to parse configuration file: " + std::string(e.what()));
        }

        // We need to store the parsed table as the root, and set the scope to
        // the root initially as well.
        root = &doc;
        scope = toml::node_view<const toml::node>(*root);
    }

    // Check that the requested field is present in the configuration file.
    void ConfigurationTable::check_field(const std::string & field) const
    {
        std::optional<std::string> value(scope.at_path(field).value<std::string>());
        if(!value)
            throw ConfigurationError("Field " + field + " not found in the configuration file.");
    }

    // Check that the requested field is present in the configuration file.
    bool ConfigurationTable::has_field(const std::string & field) const
    {
        // Check if the field is present in the scope.
        return static_cast<bool>(scope.at_path(field));
    }

    // Retrieve the requested boolean field from the configuration table.
    bool ConfigurationTable::get_bool_field(const std::string & field) const
    {
        std::optional<bool> value(scope.at_path(field).value<bool>());
        if(!value)
            throw ConfigurationError("Field " + field + " (bool) not found in the configuration file.");
        return *value;
    }

    // Retrieve the requested boolean field from the configuration table.
    bool ConfigurationTable::get_bool_field(const std::string & field, bool default_field) const
    {
        std::optional<bool> value(scope.at_path(field).value<bool>());
        if(!value)
            return default_field;
        return *value;
    }

    // Retrieve the requested string field from the configuration table.
    std::string ConfigurationTable::get_string_field(const std::string & field) const
    {
        std::optional<std::string> value(scope.at_path(field).value<std::string>());
        if(!value)
            throw ConfigurationError("Field " + field + " (string) not found in the configuration file.");
        return *value;
    }

    // Retrieve the requested string field from the configuration table.
    std::string ConfigurationTable::get_string_field(const std::string & field, const std::string & default_value) const
    {
        std::optional<std::string> value(scope.at_path(field).value<std::string>());
        if(!value)
            return default_value;
        return *value;
    }

    // Retrieve the requested vector of strings from the configuration table.
    std::vector<std::string> ConfigurationTable::get_string_vector(const std::string & field) const
    {
        std::vector<std::string> values;
        const toml::array * elements = scope.at_path(field).as_array();
        for(auto & e : *elements)
            values.push_back(*e.value<std::string>());
        return values;
    }

    // Resolve a scalar that may be given as a literal or a string reference.
    double cfg::ConfigurationTable::resolve_numeric_node(toml::node_view<const toml::node> nv) const
    {
        // Check if the node_view is a valid numeric type. If so, we simply
        // return the value directly.
        if(const auto * v = nv.as_floating_point())
            return v->get();
        if(const auto * v = nv.as_integer())
            return static_cast<double>(v->get());

        // Else, the node_view is not a numeric type and may be a reference to
        // a parameter defined elsewhere in the configuration file.
        if(const auto * sv = nv.as_string())
        {
            // Get the string value and prepare the path.
            std::string ref  = sv->get();
            std::string path = ref;

            // This is a shorthand for parameters defined in the "parameters"
            // table of the configuration file. It does a subsititution of the
            // '@' character with "parameters." to form the full path.
            // E.g., "@value" -> "parameters.value"
            if (!ref.empty() && ref.front() == '@')
                path = "parameters." + ref.substr(1);

            // Look up the path in the configuration table and check that it is
            // a valid numeric type. If so, return the value.
            auto ref_nv = lookup(path);
            if(ref_nv)
            {
                if(const auto * f = ref_nv.as_floating_point())
                    return f->get();
                if(const auto * i = ref_nv.as_integer())
                    return static_cast<double>(i->get());
            }

            // If we reach here, the reference is missing or not a numeric
            // type. Throw an exception to alert the user.
            throw ConfigurationError(
                "Reference '" + ref + "' (resolved path '" + path +
                "') is missing or not a numeric type.");
        }

        // If we reach here, the node_view is neither a numeric type nor a
        // string. Throw an exception to alert the user.
        throw ConfigurationError("Field is neither a numeric type nor a string.");
    }

    // Retrieve the requested integer field from the configuration table.
    int64_t ConfigurationTable::get_int_field(const std::string & field) const
    {
        // Get the node_view for the requested field.
        const toml::node_view<const toml::node> nv = scope.at_path(field);

        // Check that the node_view is valid.
        if(!nv)
            throw ConfigurationError("Field " + field + " (int) not found in the configuration file.");

        // Call the resolver on the node_view.
        return (int)resolve_numeric_node(nv);
    }

    // Retrieve the requested double field from the configuration table.
    double cfg::ConfigurationTable::get_double_field(const std::string& field) const
    {
        // Get the node_view for the requested field.
        const toml::node_view<const toml::node> nv = scope.at_path(field);

        // Check that the node_view is valid.
        if(!nv)
            throw ConfigurationError("Field " + field + " (double) not found in the configuration file.");

        // Call the resolver on the node_view.
        return resolve_numeric_node(nv);
    }

    std::vector<double> cfg::ConfigurationTable::get_double_vector(const std::string & field) const
    {
        // Get the node_view for the requested field. If the field is not found,
        // an exception is thrown.
        const toml::node_view<const toml::node> nv = scope.at_path(field);
        if(!nv)
            throw ConfigurationError("Field " + field + " (double[]) not found in the configuration file.");

        // Check that the node_view is an array. If not, throw an exception.
        const toml::array * elements = nv.as_array();
        if(!elements)
            throw ConfigurationError("Field " + field + " is not an array.");

        // Iterate over the elements of the array and resolve each element to a
        // double value. Store the values in a vector and return it.
        std::vector<double> values;
        values.reserve(elements->size());
        for (const auto & e : *elements)
        {
            toml::node_view<const toml::node> e_nv(e);
            values.push_back(resolve_numeric_node(e_nv));
        }

        return values;
    }

    // Get a list of all subtables matching the requested table name.
    std::vector<ConfigurationTable> ConfigurationTable::get_subtables(const std::string & table) const
    {
        // Prepare a vector to hold the subtables.
        std::vector<ConfigurationTable> tables;

        // Look up the requested table in the configuration table. It may be
        // either an array of tables or a single table.
        const auto * arr = scope.at_path(table).as_array();
        const auto * tbl = scope.at_path(table).as_table();
        
        // Check if we found an array of tables. If so, iterate over the
        // elements and construct a ConfigurationTable for each element, using
        // the root table as the root and the element as the scope.
        if(arr)
        {
            tables.reserve(arr->size());
            for(std::size_t i(0); i < arr->size(); ++i)
            {
                toml::node_view<const toml::node> nv(&(*arr)[i]);
                tables.emplace_back(root, nv);
            }
        }
        // Check if we found a single table. If so, construct a
        // ConfigurationTable using the root table as the root and the table
        // as the scope.
        else if(tbl)
        {
            tables.emplace_back(root, scope.at_path(table));
        }
        // If we found neither an array of tables nor a single table, throw an
        // exception to alert the user.
        else
        {
            throw ConfigurationError("Table " + table + " not found in the configuration file.");
        }
        return tables;
    }

    toml::node_view<const toml::node> ConfigurationTable::lookup(const std::string& path) const
    {
        // There are two places to look for the path: first in the local scope,
        // then in the root table. First, we check the local scope.
        auto sv  = scope.at_path(path);
        if(sv)
            return sv;

        // If the path is not found in the local scope, we check the root
        // table.
        if(root)
        {
            // Use node_view to avoid copying the root table.
            auto rv = toml::node_view<const toml::node>(*root).at_path(path);
            if(rv)
                return rv;
        }

        // If the path is not found in either place, return an empty node_view.
        return {};
    }

} // namespace systematics::configuration