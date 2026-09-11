// Build a synthetic medulla output file for testing validate_pair.C.
//
// Fixtures are written with ROOT rather than uproot so that the tests
// depend only on the toolchain that actually runs on the grid. (uproot
// pulls in pyarrow, which in the SL7 analysis container is linked against
// a newer glibc than the one present, so importing it fails outright.)
//
// One function covers both halves of a pair: a selection output is just
// the case with no weight tables and no non-matched tree. Every knob
// corresponds to a way real output has been seen to go wrong -- a missing
// tree, a short weight table, a file that opens cleanly with nothing in
// it -- so the tests can construct those precisely.
//
// Usage:
//   make_fixture(path, sample, tree, n_main, n_non, n_table,
//                pot, livetime, tables_csv, include_tree, include_dir)

#include <string>
#include <vector>

#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TTree.h"

namespace
{
    // A tree with the Run/Subrun/Evt branches every medulla output tree
    // carries, filled with n rows.
    void fixture_tree(TDirectory * dir, const std::string & name, Long64_t n)
    {
        dir->cd();
        TTree * t = new TTree(name.c_str(), name.c_str());
        int run = 0, subrun = 0, evt = 0;
        t->Branch("Run", &run);
        t->Branch("Subrun", &subrun);
        t->Branch("Evt", &evt);
        for(Long64_t i = 0; i < n; ++i)
        {
            run = 1;
            subrun = 1;
            evt = static_cast<int>(i);
            t->Fill();
        }
        t->Write();
    }

    // Single-bin exposure histogram, matching what medulla writes and what
    // validate_pair.C reads out of bin 1.
    void fixture_exposure(TDirectory * dir, const std::string & name, double value)
    {
        dir->cd();
        TH1D * h = new TH1D(name.c_str(), name.c_str(), 1, 0, 1);
        h->SetBinContent(1, value);
        h->Write(name.c_str());
    }

    std::vector<std::string> csv_split(const std::string & s)
    {
        std::vector<std::string> out;
        std::string cur;
        for(char c : s)
        {
            if(c == ',') { if(!cur.empty()) out.push_back(cur); cur.clear(); }
            else cur += c;
        }
        if(!cur.empty()) out.push_back(cur);
        return out;
    }
}

void make_fixture(const char * path,
                  const char * sample,
                  const char * tree,
                  Long64_t n_main,
                  Long64_t n_non,
                  Long64_t n_table,
                  double pot,
                  double livetime,
                  const char * tables,
                  bool include_tree,
                  bool include_dir)
{
    TFile * f = TFile::Open(path, "RECREATE");

    if(!include_dir)
    {
        // A file that opens fine and holds nothing of interest -- what a
        // systematics run that skipped every tree leaves behind.
        fixture_exposure(f, "placeholder", 0.0);
        f->Close();
        return;
    }

    // mkdir("a/b") creates the whole hierarchy but returns "a", not "b", so
    // look the leaf up explicitly -- otherwise everything below would land
    // one level too high, in events/ rather than events/<sample>/.
    const std::string dirname = std::string("events/") + sample;
    f->mkdir(dirname.c_str());
    TDirectory * d = f->GetDirectory(dirname.c_str());
    if(d == nullptr)
    {
        f->Close();
        gSystem->Exit(1);
    }
    fixture_exposure(d, "POT", pot);
    fixture_exposure(d, "Livetime", livetime);

    const std::string tname(tree);
    if(include_tree)
        fixture_tree(d, tname, n_main);
    if(n_non > 0)
        fixture_tree(d, tname + "_nonmatched", n_non);
    for(const std::string & t : csv_split(tables))
        fixture_tree(d, tname + "_" + t + "Tree", n_table);

    f->Close();
}
