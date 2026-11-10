// Batched integrity check of a job's staged input files.
//
// submit.sh used to run one ROOT process per input file, which cost about a
// second each in startup alone: ~25 s of a ~75 s job, more than the copy that
// fetched the files and nearly as much as the selection itself. Checking every
// file in a single process removes that, in the same way campaign.py's scan
// already does with scan_check.C.
//
// Reads a list of staged file paths (one per line) and writes one report line
// per file, in the same order:
//
//     <entries>|<path>     the file opened and recTree has <entries> entries
//     -1|<path>            the file opened but its recTree could not be read
//     BAD|<path>           the file is missing, unopenable, or a zombie
//
// The entry count is what lets the job record INPUT_EVENTS, the per-job
// normalization for run-time metrics and a completeness cross-check against the
// events actually written. A file reported BAD is dropped from the job.
//
// Usage:
//   root -l -b -q 'check_inputs.C("filelist.txt","report.txt")'

#include <fstream>
#include <string>

#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

void check_inputs(const char * filelist_path, const char * report_path)
{
    std::ifstream in(filelist_path);
    std::ofstream out(report_path);
    if(!in.is_open())
    {
        // No list means nothing can be checked; leave the report empty and let
        // the caller treat every file as unverified rather than as good.
        out.flush();
        gSystem->Exit(1);
    }

    std::string line;
    while(std::getline(in, line))
    {
        if(line.empty()) continue;

        TFile * f = TFile::Open(line.c_str());
        if(f == nullptr || f->IsZombie())
        {
            out << "BAD|" << line << "\n";
            if(f) delete f;
            continue;
        }

        TObject * obj = f->Get("recTree");
        if(obj != nullptr && obj->InheritsFrom(TTree::Class()))
            out << static_cast<TTree *>(obj)->GetEntries() << "|" << line << "\n";
        else
            out << "-1|" << line << "\n";

        f->Close();
        delete f;
    }

    out.flush();
    out.close();
    gSystem->Exit(0);
}
