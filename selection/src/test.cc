/**
 * @file test.cc
 * @brief Implementation of the test functions for the SPINE analysis
 * framework.
 * @details This file contains the implementation of the test functions
 * for the SPINE analysis framework. The tests are designed to
 * ensure that the framework's logic and functionality work as expected,
 * including interaction- and particle-level variables, matching between
 * truth and reconstructed objects, and all manner of cases surrounding
 * these functionalities.
 * @author mueller@fnal.gov
 */
#include <iostream>

#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"
#include "sbnanaobj/StandardRecord/SRParticleDLP.h"
#include "sbnanaobj/StandardRecord/SRParticleTruthDLP.h"

#include "test.h"

// Collect the offset for a given object.
template<typename T>
size_t collect_offset(const T & obj)
{
    if constexpr (std::is_same_v<T, caf::SRInteractionDLP> || std::is_same_v<T, caf::SRInteractionTruthDLP>)
    {
        return obj.particles.size();
    }
    else if constexpr (std::is_same_v<T, caf::SRParticleDLP> || std::is_same_v<T, caf::SRParticleTruthDLP>)
    {
        return 1;
    }
    else if constexpr (std::is_same_v<T, std::vector<caf::SRInteractionDLP>>)
    {
        size_t total_offset = 0;
        for(const auto & interaction : obj)
        {
            total_offset += interaction.particles.size();
        }
    }
    else if constexpr (std::is_same_v<T, std::vector<caf::SRInteractionTruthDLP>>)
    {
        size_t total_offset = 0;
        for(const auto & interaction : obj)
        {
            total_offset += interaction.particles.size();
        }
        return total_offset;
    }

    return 0;
}

// Generate a particle of the specified type.
template<typename T>
T generate_particle(int64_t id, int64_t pid)
{
    T particle;
    particle.id = id;
    particle.interaction_id = 0;
    particle.is_primary = true;
    particle.is_contained = false;
    particle.pid = pid;
    particle.ke = ENERGY_SCALE;
    particle.csda_ke = ENERGY_SCALE;
    particle.mcs_ke = ENERGY_SCALE;
    particle.calo_ke = ENERGY_SCALE;
    particle.mass = 0.0;

    if constexpr (std::is_same_v<T, caf::SRParticleTruthDLP>)
        particle.energy_init = ENERGY_SCALE;

    return particle;
}

// Generate an interaction of the specified type with the requested particle
// multiplicity.
template<typename T>
T generate_interaction(int64_t id, int64_t poffset, multiplicity_t mult, bool assign_fm)
{
    T interaction;
    interaction.id = id;

    for(size_t midx(0); midx < mult.size(); ++midx)
    {
        int64_t m = mult[midx];
        for(int64_t i(0); i < m; ++i)
        {
            if constexpr (std::is_same_v<T, caf::SRInteractionDLP>)
                interaction.particles.push_back(generate_particle<caf::SRParticleDLP>(poffset + collect_offset(interaction), midx));
            else if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLP>)
                interaction.particles.push_back(generate_particle<caf::SRParticleTruthDLP>(id * 100 + i, midx));
        }
    }

    interaction.vertex[0] = -210.0;
    interaction.vertex[1] = 0.0;
    interaction.vertex[2] = 0.0;
    if(assign_fm)
    {
        interaction.flash_times.push_back(1.0);
        interaction.is_flash_matched = 1;
    }

    return interaction;
}

// Pair two interactions together, setting the match IDs.
template<typename T, typename U>
void pair(T & left, U & right)
{
    left.match_ids.push_back(right.id);
    right.match_ids.push_back(left.id);
}

// Mark the particles as contained.
void mark_contained(caf::SRInteractionDLP * reco_interaction,
                    caf::SRInteractionTruthDLP * true_interaction)
{
    for(auto & p : reco_interaction->particles)
    {
        p.is_contained = true;
    }
    if(true_interaction)
    {
        for(auto & p : true_interaction->particles)
        {
            p.is_contained = true;
        }
    }
}

// Set the event metadata.
void write_event(caf::StandardRecord * rec, int64_t run, int64_t subrun,
                 int64_t event_num, TH1F * pot, TH1F * nevt, TTree * tree, int32_t trigger_time)
{
    // Set the size of each interaction vector.
    rec->ndlp = rec->dlp.size();
    rec->ndlp_true = rec->dlp_true.size();

    // Set the exposure information.
    rec->hdr.pot = 1.0;

    // Set the trigger time.
    rec->hdr.triggerinfo.global_trigger_time = trigger_time;

    // Set the run, subrun, and event numbers in the header.
    rec->hdr.run = run;
    rec->hdr.subrun = subrun;
    rec->hdr.evt = event_num;

    // Fill the exposure histograms.
    pot->Fill(rec->hdr.pot);
    nevt->Fill(1);

    // Fill the tree with the event data.
    tree->Fill();

    // Clear the StandardRecord interaction vectors.
    rec->dlp.clear();
    rec->dlp_true.clear();
    rec->ndlp = 0;
    rec->ndlp_true = 0;
}

// Read the event data from the TTree at the specified path.
std::vector<row_t> read_event_data(const std::string & name)
{
    TTree * t = dynamic_cast<TTree *>(gDirectory->Get(name.c_str()));
    if(!t)
    {
        throw std::runtime_error("Could not find TTree " + name + " in the current directory.");
    }

    // Triplet metadata.
    Int_t run, subrun, event;

    // Everything else is a double. Retrieve the list of branches.
    auto branches = t->GetListOfBranches();
    std::vector<std::string> branch_names;
    for(auto const & branch : *branches)
    {
        auto b = dynamic_cast<TBranch *>(branch);
        if(!b || std::string(b->GetName()) == "Run" || std::string(b->GetName()) == "Subrun" || std::string(b->GetName()) == "Evt")
            continue;
        branch_names.push_back(b->GetName());
    }

    // Double array to hold the values for each branch.
    double values[branch_names.size()];

    // Set the branch addresses.
    t->SetBranchAddress("Run", &run);
    t->SetBranchAddress("Subrun", &subrun);
    t->SetBranchAddress("Evt", &event);
    for(size_t i = 0; i < branch_names.size(); ++i)
    {
        t->SetBranchAddress(branch_names[i].c_str(), &values[i]);
    }

    // Read the entries from the TTree.
    std::vector<row_t> rows;
    for(int i = 0; i < t->GetEntries(); ++i)
    {
        t->GetEntry(i);
        row_t row;
        row["Run"] = run;
        row["Subrun"] = subrun;
        row["Evt"] = event;

        // Fill the row with the values from the branches.
        for(size_t j = 0; j < branch_names.size(); ++j)
        {
            row[branch_names[j]] = values[j];
        }
        rows.push_back(row);
    }

    return rows;
}

// Match the (Run, Subrun, Evt) metadata between a row_t object and a
// condition_t object.
bool match_metadata(const row_t & row, const condition_t & condition)
{
    for(const auto & field : condition.second)
    {
        if(field.first == "Run" || field.first == "Subrun" || field.first == "Evt")
        {
            if(row.at(field.first) != field.second)
                return false; // Metadata mismatch.
        }
    }
    return true;
}

// Match a set of row_t objects against a set of condition_t objects.
void match_conditions(const std::vector<row_t> & rows, const std::vector<condition_t> & conditions)
{
    for(const auto & condition : conditions)
    {
        bool found = false;
        for(const auto & row : rows)
        {
            if(match_metadata(row, condition))
            {
                found = true;
                // Check if all fields in the condition match the corresponding
                // fields in the row.
                bool match = true;
                for(const auto & field : condition.second)
                {
                    if(row.find(field.first) == row.end()
                       || (!std::isnan(field.second) && row.at(field.first) != field.second)
                       || (std::isnan(field.second) && !std::isnan(row.at(field.first))))
                    {
                        match = false;
                        break; // Mismatch found, no need to check further.
                    }
                }
                if(match)
                {
                    if(condition.first.find('!') == std::string::npos)
                    {
                        std::cout << "\033[32mValidation passed:\033[0m   " << condition.first << "." << std::endl;
                        found = true;
                        break;
                    }
                    else
                    {
                        // If the condition starts with '!', it means we expect it to not match.
                        std::cout << "\033[31mValidation failed:\033[0m   " << condition.first.substr(1) << "." << std::endl;
                        found = true;
                        break;
                    }
                }
                else if(found && !match)                
                {
                    std::cerr << "\033[33mValidation mismatch:\033[0m " << condition.first << "." << std::endl;
                    
                    // Print the fields that are mismatched.
                    for(const auto & field : condition.second)
                    {
                        if(row.find(field.first) == row.end() || row.at(field.first) != field.second)
                        {
                            std::cout << "    " << field.first
                                        << " - expected: " << field.second
                                        << ", got: " << (row.find(field.first) != row.end() ? std::to_string(row.at(field.first)) : "N/A") << std::endl;
                        }
                    }
                }
            }
        }
        if(!found && condition.first.find('!') == std::string::npos)
            std::cout << "\033[31mValidation failed:\033[0m   " << condition.first << "." << std::endl;
        else if(!found && condition.first.find('!') != std::string::npos)
            std::cout << "\033[32mValidation passed:\033[0m   " << condition.first.substr(1) << "." << std::endl;
    }
}

// Explicit template instantiations for the functions defined above.
template caf::SRInteractionDLP generate_interaction<caf::SRInteractionDLP>(int64_t, int64_t, multiplicity_t, bool);
template caf::SRInteractionTruthDLP generate_interaction<caf::SRInteractionTruthDLP>(int64_t, int64_t, multiplicity_t, bool);
template void pair<caf::SRInteractionDLP, caf::SRInteractionTruthDLP>(caf::SRInteractionDLP &, caf::SRInteractionTruthDLP &);
template void pair<caf::SRInteractionTruthDLP, caf::SRInteractionDLP>(caf::SRInteractionTruthDLP &, caf::SRInteractionDLP &);
template void pair<caf::SRParticleDLP, caf::SRParticleTruthDLP>(caf::SRParticleDLP &, caf::SRParticleTruthDLP &);
template void pair<caf::SRParticleTruthDLP, caf::SRParticleDLP>(caf::SRParticleTruthDLP &, caf::SRParticleDLP &);