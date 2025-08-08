/**
 * @file configuration.h
 * @brief Header of the ConfigurationTable class and related functions.
 * @details This file contains the header for the ConfigurationTable class and
 * related functions. The ConfigurationTable class is used to read and
 * interface with the TOML-based configuration file. The configuration file
 * consolidates all the parameters that may be analysis-specific or that may
 * change frequently. The main backend employed in this project is the toml++
 * library.
 * @author mueller@fnal.gov
 */
#ifndef CONFIGURATION_H
#define CONFIGURATION_H
#include <toml++/toml.h>

/**
 * @namespace cfg
 * @brief Namespace for functions that read and interface with the TOML-based
 * configuration file.
 * @details This namespace contains functions that read and interface with the
 * TOML-based configuration file. The configuration file consolidates all the
 * parameters that may be analysis-specific or that may change frequently. The
 * main backend employed in this project is toml++.
 */
namespace cfg
{
    /**
     * @class A class inheriting from std::exception that is used as an error
     * type for configuration-related errors.
     * @details This class is used as an error type for configuration-related
     * errors. The class inherits from std::exception and is used to throw
     * exceptions when the configuration file is not found or when a requested
     * field is not present in the configuration file.
     * @see std::exception
     */
    class ConfigurationError : public std::exception
    {
    public:
        /**
         * @brief Constructor for the ConfigurationError class.
         * @details This constructor takes a string as an argument and stores it
         * as the error message.
         * @param message The error message.
         */
        ConfigurationError(const std::string & message) : message(message) {}

        /**
         * @brief Get the error message.
         * @details This function returns the error message.
         * @return The error message.
         */
        const char * what() const noexcept override { return message.c_str(); }
    
    private:
        std::string message; ///< The error message.
    };

    class ConfigurationTable
    {
    public:
        /**
         * @brief Standard constructor for the ConfigurationTable class.
         * @details This constructor initializes an empty TOML configuration
         * table.
         * @return void
         * @throw None
         */
        ConfigurationTable() = default;

        /**
         * @brief Constructor for the Configuration class.
         * @details This constructor takes a TOML configuration table as an
         * argument and stores it as a class member. This constructor is meant
         * to be used in instances where this class is encapsulating a
         * sub-table of the original configuration file.
         * @param config The TOML configuration table.
         */
        ConfigurationTable(const toml::table & config);

        /**
         * @brief Set the configuration table to a table loaded from a file.
         * @details This function sets the configuration table to a table
         * loaded from a file. The function reads the configuration file using
         * the toml++ library and stores the configuration table as a class
         * member. The function also validates the configuration file by
         * calling the function @ref validate().
         * @param path The path to the configuration file.
         * @return void
         * @throw ConfigurationError
         */
        void set_config(const std::string & path);

        /**
         * @brief Check that the requested field is present in the configuration
         * file.
         * @details This function checks that the requested field is present in the
         * configuration file. If the field is not present, the function throws an
         * exception.
         * @param config The TOML configuration table.
         * @param field The field that is requested.
         * @return void
         * @throw ConfigurationError
         */
        void check_field(const std::string & field) const;

        /**
         * @brief Check that the requested field is present in the configuration
         * file.
         * @details This function checks that the requested field is present in the
         * configuration file. It returns a boolean value indicating whether the
         * field is present and does not throw an exception.
         * @param field The field that is requested.
         * @return A boolean value indicating whether the field is present.
         */
        bool has_field(const std::string & field) const;

        /**
         * @brief Get the requested boolean field from the ConfigurationTable.
         * @details This function gets the requested boolean field from the
         * ConfigurationTable. If the field is not present, the function throws
         * an exception.
         * @param field The name of the field that is requested.
         * @return The value of the requested boolean field.
         * @throw ConfigurationError
         */
        bool get_bool_field(const std::string & field) const;

        /**
         * @brief Get the requested boolean field from the ConfigurationTable.
         * @details This function gets the requested boolean field from the
         * ConfigurationTable. If the field is not present, the provided
         * default value is returned instead of throwing an exception.
         * @param field The name of the field that is requested.
         * @param default_value The default value to return if the field is not
         * present.
         * @return The value of the requested boolean field.
         */
        bool get_bool_field(const std::string & field, bool default_value) const;

        /**
         * @brief Get the requested string field from the ConfigurationTable.
         * @details This function gets the requested string field from the
         * ConfigurationTable. If the field is not present, the function throws
         * an exception.
         * @param field The name of the field that is requested.
         * @return The value of the requested string field.
         * @throw ConfigurationError
         */
        std::string get_string_field(const std::string & field) const;

        /**
         * @brief Get the requested string field from the ConfigurationTable.
         * @details This function gets the requested string field from the
         * ConfigurationTable. If the field is not present, the provided
         * default value is returned instead of throwing an exception.
         * @param field The name of the field that is requested.
         * @param default_value The default value to return if the field is not
         * present.
         * @return The value of the requested string field.
         */
        std::string get_string_field(const std::string & field, const std::string & default_value) const;

        /**
         * @brief Get a list of all strings matching the requested field name.
         * @details This function gets a list of all strings matching the requested
         * field name. The function returns a vector of strings.
         * @param field The field that is requested.
         * @return A vector of strings.
         * @throw ConfigurationError
         */
        std::vector<std::string> get_string_vector(const std::string & field) const;

        /**
         * @brief Get the requested integer field from the ConfigurationTable.
         * @details This function gets the requested integer field from the
         * ConfigurationTable. If the field is not present, the function throws
         * an exception.
         * @param field The name of the field that is requested.
         * @return The value of the requested integer field.
         * @throw ConfigurationError
         */
        int64_t get_int_field(const std::string & field) const;

        /**
         * @brief Get the requested double field from the ConfigurationTable.
         * @details This function gets the requested double field from the
         * ConfigurationTable. If the field is not present, the function throws
         * an exception.
         * @param field The name of the field that is requested.
         * @return The value of the requested double field.
         * @throw ConfigurationError
         */
        double get_double_field(const std::string & field) const;

        /**
         * @brief Get a list of all doubles matching the requested field name.
         * @details This function gets a list of all doubles matching the
         * requested field name. The function returns a vector of doubles.
         * @param field The field that is requested.
         * @return A vector of doubles.
         * @throw ConfigurationError
         */
        std::vector<double> get_double_vector(const std::string & field) const;

        /**
         * @brief Get a list of all subtables matching the requested table name.
         * @details This function gets a list of all subtables matching the
         * requested table name. The function returns a vector of ConfigurationTable
         * objects.
         * @param table The table that is requested.
         * @return A vector of ConfigurationTable objects.
         * @throw ConfigurationError
         * @see ConfigurationTable
         */
        std::vector<ConfigurationTable> get_subtables(const std::string & table) const;
    
    private:
        toml::table config; ///< The TOML configuration table.
    };
}
#endif // CONFIGURATION_H