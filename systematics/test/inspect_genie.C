/**
 * @file inspect_genie.C
 * @brief Inspect the GENIE event records written alongside selected candidates.
 * @details This macro prints one line per entry of the record tree so that the
 * records can be eyeballed against the selected candidate tree they are meant
 * to be aligned with. It relies on the GENIE dictionaries being loadable, which
 * is the same requirement the writing code has, so it doubles as a check that
 * the environment is set up correctly.
 *
 * A default-initialized record is recognisable by its null 'event' pointer:
 * genie::NtpMCEventRecord's constructor leaves it null, so an entry that shows
 * "(default-initialized)" is one where no GENIE record was found for the
 * corresponding selected candidate.
 *
 * Usage (with the environment from setup.sh):
 *   root -l -b -q 'inspect_genie.C("output.root", "events/full", "selected")'
 *
 * The second argument is the directory holding the trees and the third is the
 * selected candidate tree's name; the record tree is that name + "_genieTree".
 *
 * @author mueller@fnal.gov
 */
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"

#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Interaction/Interaction.h"

void inspect_genie(const char * path,
                   const char * directory = "",
                   const char * selected = "selected",
                   Long64_t max_entries = 20)
{
    TFile * file = TFile::Open(path, "READ");
    if(file == nullptr || file->IsZombie())
    {
        printf("Could not open %s\n", path);
        return;
    }

    TDirectory * dir = std::string(directory).empty()
        ? (TDirectory *) file
        : file->GetDirectory(directory);
    if(dir == nullptr)
    {
        printf("No such directory: %s\n", directory);
        return;
    }

    const std::string record_name = std::string(selected) + "_genieTree";
    TTree * records = (TTree *) dir->Get(record_name.c_str());
    TTree * candidates = (TTree *) dir->Get(selected);

    if(records == nullptr)
    {
        printf("No record tree '%s' in %s. Either 'general.store_genie_evt_rec' was\n"
               "not set, or one of the checks disabled it -- look for the warning in\n"
               "the run_systematics output.\n", record_name.c_str(), directory);
        return;
    }

    /**
     * @brief The alignment check that matters.
     * @details Entry N of the record tree is meant to belong to entry N of the
     * selected candidate tree, so the two must have the same length. A mismatch
     * means records are being attributed to the wrong candidates.
     */
    printf("Record tree '%s': %lld entries\n", record_name.c_str(), records->GetEntries());
    if(candidates != nullptr)
    {
        printf("Selected tree '%s': %lld entries\n", selected, candidates->GetEntries());
        printf("ALIGNMENT: %s\n\n",
               candidates->GetEntries() == records->GetEntries() ? "OK" : "*** MISMATCH ***");
    }
    else
        printf("(selected tree '%s' not found; skipping the alignment check)\n\n", selected);

    genie::NtpMCEventRecord * record = nullptr;
    records->SetBranchAddress("GenieEvtRec", &record);

    // If the candidate tree carries the neutrino energy, print it alongside the
    // probe energy from the record: the two should agree for a correct match.
    double true_energy = 0;
    const bool have_energy = candidates != nullptr
        && candidates->GetBranch("true_neutrino_energy") != nullptr;
    if(have_energy)
        candidates->SetBranchAddress("true_neutrino_energy", &true_energy);

    const Long64_t n = std::min(records->GetEntries(), max_entries);
    for(Long64_t i(0); i < n; ++i)
    {
        records->GetEntry(i);
        if(have_energy)
            candidates->GetEntry(i);

        if(record == nullptr || record->event == nullptr)
        {
            printf("  [%4lld] (default-initialized -- no GENIE record for this candidate)\n", i);
            continue;
        }

        genie::GHepParticle * probe = record->event->Probe();
        printf("  [%4lld] probe pdg=%6d  E=%8.4f GeV  particles=%3d",
               i,
               probe != nullptr ? probe->Pdg() : 0,
               probe != nullptr ? probe->E() : 0.0,
               record->event->GetEntries());

        if(have_energy)
            printf("   | candidate true_neutrino_energy=%8.4f  %s",
                   true_energy,
                   (probe != nullptr && std::abs(probe->E() - true_energy) < 1e-3)
                       ? "MATCH" : "<-- CHECK");
        printf("\n");
    }

    if(records->GetEntries() > n)
        printf("  ... %lld more entries\n", records->GetEntries() - n);

    file->Close();
}
