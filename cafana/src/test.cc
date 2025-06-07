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
    particle.pid = pid;
    particle.ke = ENERGY_SCALE;
    particle.csda_ke = ENERGY_SCALE + 1.0;
    particle.mcs_ke = ENERGY_SCALE + 2.0;
    particle.calo_ke = ENERGY_SCALE + 3.0;

    return particle;
}

// Generate an interaction of the specified type with the requested particle
// multiplicity.
template<typename T>
T generate_interaction(int64_t id, int64_t poffset, multiplicity_t mult)
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

    return interaction;
}

// Pair two interactions together, setting the match IDs.
template<typename T, typename U>
void pair(T & left, U & right)
{
    left.match_ids.push_back(right.id);
}

// Set the event metadata.
void write_event(caf::StandardRecord * rec, int64_t run, int64_t subrun,
                 int64_t event_num, TH1F * pot, TH1F * nevt, TTree * tree)
{
    // Set the size of each interaction vector.
    rec->ndlp = rec->dlp.size();
    rec->ndlp_true = rec->dlp_true.size();

    // Set the exposure information.
    rec->hdr.pot = 1.0;

    // Set the run, subrun, and event numbers in the header.
    rec->hdr.run = run;
    rec->hdr.subrun = subrun;
    rec->hdr.evt = event_num;

    // Fill the exposure histograms.
    pot->Fill(rec->hdr.pot);
    nevt->Fill(1);

    // Fill the tree with the event data.
    tree->Fill();
}

// Explicit template instantiations for the functions defined above.
template caf::SRInteractionDLP generate_interaction<caf::SRInteractionDLP>(int64_t, int64_t, multiplicity_t);
template caf::SRInteractionTruthDLP generate_interaction<caf::SRInteractionTruthDLP>(int64_t, int64_t, multiplicity_t);
template void pair<caf::SRInteractionDLP, caf::SRInteractionTruthDLP>(caf::SRInteractionDLP &, caf::SRInteractionTruthDLP &);
template void pair<caf::SRInteractionTruthDLP, caf::SRInteractionDLP>(caf::SRInteractionTruthDLP &, caf::SRInteractionDLP &);