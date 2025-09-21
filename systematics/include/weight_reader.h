/**
 * @file weight_reader.h
 * @brief Header file for the WeightReader class.
 * @details This class provides methods to read and access weight information
 * from structured and flat CAF files. It can accept a file pattern, a list of
 * files, or a single file as input. The class provides methods to get the
 * number of neutrinos, weight groups, and universes, as well as the weight
 * values themselves.
 * @author mueller@fnal.gov
 */
#ifndef WEIGHT_READER_H
#define WEIGHT_READER_H
#include <chrono>

#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"

namespace sys
{
    /**
     * @class WeightReader
     * @brief A class to read and access weight information from CAF files.
     * @details This class provides methods to read and access weight information
     * from structured and flat CAF files. It can accept a file pattern, a list of
     * files, or a single file as input. The class provides methods to advance to
     * the next entry, set the weight group index, and access metadata such as
     * run, subrun, and event numbers. It also provides methods to get the number
     * of neutrinos, weight groups, and universes, as well as the weight values
     * themselves.
     */
    class WeightReader
    {
        public:
        
        /**
         * @brief Constructor for the WeightReader class.
         * @details This constructor initializes the WeightReader object
         * with the input file or file pattern. It sets up the TChain using the
         * relevant methods and initializes the machinery to read the weight
         * information.
         * @param input The input file or file pattern to read.
         */
        explicit WeightReader(const std::string & input);
        
        /**
         * @brief Destructor for the WeightReader class.
         * @details This destructor cleans up the resources used by the
         * WeightReader object.
         */
        ~WeightReader() = default;

        /**
         * @brief Advance to the next entry in the TChain.
         * @details This method advances the TChain to the next entry and
         * updates the entry index. It returns true if successful, false
         * otherwise. Internally, it will correctly handle differences between
         * structured and flat CAF files.
         * @return True if successful, false otherwise.
         */
        bool next();

        /**
         * @brief Set the weight group index.
         * @details This method sets the weight group index for the current
         * entry. This determines the systematic parameter to be accessed
         * for the current event.
         * @param index The index of the weight group to set.
         * @return void
         */
        void set(size_t index);

        /**
         * @brief Get the run number.
         * @details This method returns the run number for the current
         * entry.
         * @return The run number.
         */
        uint32_t get_run() const { return **run; }

        /**
         * @brief Get the subrun number.
         * @details This method returns the subrun number for the current
         * entry.
         * @return The subrun number.
         */
        uint32_t get_subrun() const { return **subrun; }

        /**
         * @brief Get the event number.
         * @details This method returns the event number for the current
         * entry.
         * @return The event number.
         */
        uint32_t get_event() const { return **event; }

        /**
         * @brief Accessor method for the number of neutrinos.
         * @details This method returns the number of neutrinos for the
         * current entry. It handles both structured and flat CAF files.
         * @return The number of neutrinos.
         */
        uint32_t get_nnu() const;

        /**
         * @brief Accessor method for the number of weight groups.
         * @details This method returns the number of weight groups for the
         * current entry. It handles both structured and flat CAF files.
         * @param i The index of the neutrino.
         * @return The number of weight groups for the specified neutrino.
         */
        uint32_t get_nwgt(Int_t i) const;

        /**
         * @brief Accessor method for the number of universes.
         * @details This method returns the number of universes for the
         * specified neutrino and weight group. It handles both structured
         * and flat CAF files.
         * @param idn The index of the neutrino.
         * @return The number of universes for the specified neutrino and
         * weight group.
         */
        uint32_t get_nuniv(size_t idn) const;

        /**
         * @brief Accessor method for the weight value.
         * @details This method returns the weight value for the specified
         * neutrino, weight group, and universe. It handles both structured
         * and flat CAF files.
         * @param idn The index of the neutrino.
         * @param idu The index of the universe.
         * @return The weight value for the specified neutrino, weight group,
         * and universe.
         * @note This method assumes that the weight group index has been set
         * using the @ref set method.
         */
        float get_weight(size_t idn, size_t idu) const;

        /**
         * @brief Accessor method for the neutrino energy.
         * @details This method returns the neutrino energy for the
         * specified neutrino. It handles both structured and flat CAF files.
         * @param idn The index of the neutrino.
         * @return The neutrino energy for the specified neutrino.
         */
        float get_energy(size_t idn) const;

        private:

        /**
         * @brief A simple progress bar for the TChain.
         * @details This method provides a simple progress bar for the TChain
         * to indicate the progress of the reading process.
         * @param entry The current entry index.
         * @param total The total number of entries.
         * @return void
         */
        void progress_bar(size_t entry, size_t total) const;

        bool isflat; // Flag to indicate if the input file is flat or structured
        TChain chain; // TChain to hold the input files
        size_t entry; // Current entry index in the TChain

        std::unique_ptr<TTreeReader> reader; // TTreeReader for structured CAF files
        
        // Metadata
        std::unique_ptr<TTreeReaderValue<uint32_t>> run; // Run number
        std::unique_ptr<TTreeReaderValue<uint32_t>> subrun; // Subrun number
        std::unique_ptr<TTreeReaderValue<uint32_t>> event; // Event number

        // Event-level indexing
        Int_t nnu; // Number of neutrinos
        std::unique_ptr<TTreeReaderValue<uint64_t>> nnu_structured; // Number of neutrinos for structured CAF files

        // Neutrino-level indexing
        Int_t nwgt[10]; // Number of weight groups for each neutrino
        Int_t iwgt[10]; // Index of the weight group for each neutrino
        size_t idx; // Index of the current weight group
        Float_t nu_energy[10]; // Neutrino energy for each neutrino
        std::unique_ptr<TTreeReaderArray<Float_t>> nu_energy_structured; // Neutrino energy for structured CAF files

        // Systematic-level indexing
        Int_t nuniv[15000]; // Number of universes for each weight group
        Int_t iuniv[15000]; // Index of the universe for each weight group
        Float_t wgts[500000]; // Weight values for each universe

        // MC-truth branch
        std::unique_ptr<TTreeReaderArray<caf::SRTrueInteraction>> mc; // MC-truth data for structured CAF files

        // Progress bar timestamp
        mutable std::chrono::steady_clock::time_point progress_start_time; // Start time for the progress bar
        mutable bool progress_started = false; // Flag to indicate if the progress bar has started
        mutable int last_printed_percent = -1; // Last printed percent for the progress bar
    };
} // namespace sys
#endif  // WEIGHT_READER_H