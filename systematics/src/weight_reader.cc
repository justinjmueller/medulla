/**
 * @file
 * @brief A simple weight reader for CAF / flat CAF files.
 * @details This class reads the weight information from a TTree and provides
 * accessor methods for retrieving metadata (such as run, subrun, and event
 * numbers) and the weights themselves by index.
 * @author mueller@fnal.gov
 */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <sstream>

#include "weight_reader.h"

#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"

// Constructor for the WeightReader class.
sys::WeightReader::WeightReader(const std::string & input)
: entry(0),
  chain("recTree"),
  idx(0),
  progress_started(false)
{
    if(input.find("*") != std::string::npos)
    {
        // Input is a pattern for a set of files
        isflat = input.find("flat") != std::string::npos;
        chain.Add(input.c_str());
    }
    else if(input.find(".txt") != std::string::npos)
    {
        // Input is a .txt file containing a list of files
        std::ifstream infile(input);
        std::string line;
        while(std::getline(infile, line))
            chain.Add(line.c_str());
        isflat = line.find("flat") != std::string::npos;
    }
    else
    {
        // Input is a single .root file
        isflat = input.find("flat") != std::string::npos;
        chain.Add(input.c_str());
    }
    
    if(isflat)
    {
        // Set the TTreeReader to nullptr for flat CAF files
        reader = nullptr;

        // Metadata branches
        chain.SetBranchAddress("rec.hdr.run", &run_value);
        chain.SetBranchAddress("rec.hdr.subrun", &subrun_value);
        chain.SetBranchAddress("rec.hdr.evt", &event_value);

        // Event-level indexing
        chain.SetBranchAddress("rec.mc.nu..length", &nnu);

        // Neutrino-level indexing
        chain.SetBranchAddress("rec.mc.nu.wgt..length", nwgt);
        chain.SetBranchAddress("rec.mc.nu.wgt..idx", &iwgt);
        chain.SetBranchAddress("rec.mc.nu.E", &nu_energy);

        // Systematic-level indexing
        chain.SetBranchAddress("rec.mc.nu.wgt.univ..length", &nuniv);
        chain.SetBranchAddress("rec.mc.nu.wgt.univ..idx", &iuniv);
        chain.SetBranchAddress("rec.mc.nu.wgt.univ", wgts);
        chain.GetEntry(0);
    }
    else
    {
        // Create the TTreeReader
        reader = std::make_unique<TTreeReader>(&chain);

        // Metadata branches
        run = std::make_unique<TTreeReaderValue<uint32_t>>(*reader, "rec.hdr.run");
        subrun = std::make_unique<TTreeReaderValue<uint32_t>>(*reader, "rec.hdr.subrun");
        event = std::make_unique<TTreeReaderValue<uint32_t>>(*reader, "rec.hdr.evt");

        // MC-truth branches
        nnu_structured = std::make_unique<TTreeReaderValue<uint64_t>>(*reader, "rec.mc.nnu");
        mc = std::make_unique<TTreeReaderArray<caf::SRTrueInteraction>>(*reader, "rec.mc.nu");
        nu_energy_structured = std::make_unique<TTreeReaderArray<Float_t>>(*reader, "rec.mc.nu.E");
    }
    //reader->Next();
}

// Advance to the next entry in the TChain.
bool sys::WeightReader::next()
{
    // Update the progress bar
    this->progress_bar(entry+1, chain.GetEntries());

    // Handler for flat CAF files
    if(isflat)
    {
        if(!chain.GetTree() || chain.GetEntry(entry) == 0)
            return false;
    }
    // Handler for structured CAF files
    else
    {
        if(!reader || !reader->Next())
            return false;
        run_value = **run;
        subrun_value = **subrun;
        event_value = **event;
    }
    ++entry;

    return true;
}

// Set the weight group index.
void sys::WeightReader::set(size_t index)
{
    idx = index;
}

// Accessor method for the number of neutrinos.
uint32_t sys::WeightReader::get_nnu() const
{
    return isflat ? nnu : **nnu_structured;
}

// Accessor method for the number of weight groups.
uint32_t sys::WeightReader::get_nwgt(Int_t i) const
{
    if(i < 0 || i >= (Int_t)this->get_nnu())
        throw std::out_of_range("WeightReader: Index out of range in 'get_nwgt()'");
    
    return isflat ? nwgt[i] : (*mc)[i].wgt.size();
}

// Accessor method for the number of universes.
uint32_t sys::WeightReader::get_nuniv(size_t idn) const
{
    if(idn >= get_nnu() || idx >= get_nwgt(idn))
        throw std::out_of_range("WeightReader: Index out of range in 'get_nuniv()'");

    if(isflat)
    {
        size_t n = iwgt[idn] + idx;
        if(n >= 10000)
            throw std::out_of_range("WeightReader: Weight group index out of range in 'get_nuniv()'");
        return nuniv[n];
    }
    else
    {
        return (*mc)[idn].wgt[idx].univ.size();
    }
}

// Accessor method for the weight value.
float sys::WeightReader::get_weight(size_t idn, size_t idu) const
{
    if(isflat)
    {
        if(idn >= get_nnu())
            throw std::out_of_range("WeightReader: Neutrino index out of range in 'get_weight()'");

        size_t n = iwgt[idn] + idx;
        if(n >= 10000)        
            throw std::out_of_range("WeightReader: Weight group index out of range in 'get_weight()'");
        
        size_t univ_offset = iuniv[n] + idu;
        if(univ_offset >= 150000 || iuniv[n] < 0)
            throw std::out_of_range("WeightReader: Universe index out of range in 'get_weight()'");
        
        return wgts[univ_offset + idu];
    }
    else
        return (*mc)[idn].wgt[idx].univ[idu];
}

// Accessor method for the neutrino energy.
float sys::WeightReader::get_energy(size_t idn) const
{
    return isflat ? nu_energy[idn] : (*nu_energy_structured)[idn];
}

// Simple progress bar for the TChain.
void sys::WeightReader::progress_bar(size_t entry, size_t total) const
{
    // Start the clock if it hasn't been started yet.
    if(!progress_started)
    {
        progress_start_time = std::chrono::steady_clock::now();
        progress_started = true;
    }

    // Calculate and display fractional progress.
    float percent = (float)entry / total;
    int percent_int = static_cast<int>(percent*1000.0);
    if(percent_int == last_printed_percent && entry != total)
        return;
    last_printed_percent = percent_int;

    // Clear the line and print progress bar
    std::cout << "\r\033[K[";  // \r = carriage return, \033[K = clear to end of line

    int pos = static_cast<int>(50 * percent);
    for(int i = 0; i < 50; ++i)
    {
        if(i < pos) std::cout << "=";
        else if(i == pos) std::cout << ">";
        else std::cout << " ";
    }

    std::cout << "] " << std::fixed << std::setprecision(2) << percent * 100.0 << "%  ";

    // Calculate time elapsed and estimated time remaining.
    auto now = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(now - progress_start_time).count();
    double eta = percent > 0.0 ? elapsed / percent - elapsed : 0.0;

    auto format_time = [](double seconds) -> std::string
    {
        int h = static_cast<int>(seconds) / 3600;
        int m = (static_cast<int>(seconds) % 3600) / 60;
        double s = seconds - h * 3600 - m * 60;

        std::ostringstream oss;
        oss << std::setfill('0') << std::setw(2) << h << ":"
            << std::setw(2) << m << ":"
            << std::setw(4) << std::fixed << std::setprecision(1) << s;
        return oss.str();
    };

    std::cout << "Elapsed: " << format_time(elapsed)
              << ", ETA: " << format_time(eta) << std::flush;

    if(entry == total)
    {
        std::cout << std::endl;
        progress_started = false;
    }
}