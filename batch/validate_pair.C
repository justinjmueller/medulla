// Worker-side validation of a job's (selection, systematics) output pair,
// run by submit.sh on the grid node *before* either file is transferred
// back. The point is to never put a bad or half-complete result on dCache:
// the two artifacts are copied back only if this macro exits 0.
//
// The failure this primarily exists to catch is asymmetric output --
// output.root fine, output_sys.root missing or silently empty. That mode is
// structural rather than accidental: systematics/src/main.cc skips any
// configured tree that is absent from its input and still returns 0, so a
// run that produced nothing at all reports success. Checking exit codes
// alone therefore cannot detect it; the output has to be inspected.
//
// Checks are driven by validation_manifest.txt (written at project creation
// by utilities.write_validation_manifest from the same object that becomes
// systematics.toml), which says which trees must exist and, for each,
// whether it was a plain copy or had weight tables added.
//
// A machine-readable KEY=VALUE report is written to report_path for the
// job's validation record, and is produced even when a check fails -- the
// record is the only trustworthy per-job signal, since a batched grid
// process exits with the status of whatever ran last.
//
// Usage:
//   root -l -b -q 'validate_pair.C("output.root","output_sys.root",
//                                  "sbnd","validation_manifest.txt",
//                                  "report.txt")'
//
// Exit codes (recorded as VALIDATE_EXIT_CODE, and mapped to a failure
// category by the caller):
//   0  pass
//   1  usage or manifest error
//   2  selection output bad (open/zombie/missing sample dir/exposure/tree)
//   3  systematics output bad (open/zombie/missing sample dir)
//   4  exposure (POT/Livetime) mismatch between the two files
//   5  tree-level failure (missing tree, entry-count mismatch, empty match)

#include <cmath>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TKey.h"
#include "TList.h"
#include "TSystem.h"
#include "TTree.h"

namespace
{
    // One manifest line: what the systematics step owed us for one tree.
    struct ManifestEntry
    {
        std::string origin;                    // events/<sample>/<tree>
        std::string name;                      // <tree>
        std::string action;                    // "copy" or "add_weights"
        std::vector<std::string> table_types;  // e.g. multisim, multisigma
    };

    std::vector<std::string> split_on(const std::string & s, char sep)
    {
        std::vector<std::string> out;
        std::string cur;
        for(char c : s)
        {
            if(c == sep) { out.push_back(cur); cur.clear(); }
            else cur += c;
        }
        out.push_back(cur);
        return out;
    }

    // Fetch an object only if it is really of the expected class. A
    // straight C-style cast would happily reinterpret, say, the POT
    // histogram as a TTree and then crash on first use.
    TObject * get_checked(TDirectory * dir, const std::string & name, TClass * cls)
    {
        if(dir == nullptr) return nullptr;
        TObject * obj = dir->Get(name.c_str());
        if(obj == nullptr || !obj->InheritsFrom(cls)) return nullptr;
        return obj;
    }

    // Entries in a tree, or -1 when the tree is absent. -1 is used rather
    // than 0 so that "missing" stays distinguishable from "present but
    // empty" -- the two mean very different things for a _nonmatched tree.
    Long64_t tree_entries(TDirectory * dir, const std::string & name)
    {
        TObject * obj = get_checked(dir, name, TTree::Class());
        if(obj == nullptr) return -1;
        return static_cast<TTree *>(obj)->GetEntries();
    }

    // Bin-1 content of a 1-D exposure histogram, or NaN when absent.
    double exposure_value(TDirectory * dir, const std::string & name)
    {
        TObject * obj = get_checked(dir, name, TH1::Class());
        if(obj == nullptr) return std::nan("");
        return static_cast<TH1 *>(obj)->GetBinContent(1);
    }

    // Stamp the verdict and leave. Every exit goes through here because
    // gSystem->Exit() calls exit(), which does not run the destructor of a
    // local ofstream -- writing the report and exiting separately would
    // truncate it on precisely the failure paths the record exists for.
    void finish(std::ofstream & rep, const std::string & status, int code)
    {
        rep << "STATUS=" << status << "\n";
        rep.flush();
        rep.close();
        gSystem->Exit(code);
    }
}

void validate_pair(const char * nosyst_path,
                   const char * syst_path,
                   const char * sample,
                   const char * manifest_path,
                   const char * report_path,
                   double nonmatched_warn_frac = 0.01)
{
    std::ofstream rep(report_path);
    rep << "SAMPLE=" << sample << "\n";

    // -----------------------------------------------------------------
    // Manifest
    // -----------------------------------------------------------------
    // The manifest describes every sample in the project; this job ran
    // exactly one, so keep only the entries under this sample's directory.
    const std::string prefix = std::string("events/") + sample + "/";
    std::vector<ManifestEntry> expected;
    {
        std::ifstream man(manifest_path);
        if(!man.is_open())
        {
            rep << "ERROR=manifest_unreadable:" << manifest_path << "\n";
            finish(rep, "manifest_error", 1);
        }
        std::string line;
        while(std::getline(man, line))
        {
            if(line.empty()) continue;
            std::vector<std::string> f = split_on(line, '|');
            if(f.size() < 4)
            {
                rep << "ERROR=manifest_malformed_line:" << line << "\n";
                finish(rep, "manifest_error", 1);
            }
            if(f[0].rfind(prefix, 0) != 0) continue;  // another sample
            ManifestEntry e;
            e.origin = f[0];
            e.name   = f[1];
            e.action = f[2];
            if(!f[3].empty()) e.table_types = split_on(f[3], ',');
            expected.push_back(e);
        }
    }
    rep << "N_EXPECTED_TREES=" << expected.size() << "\n";
    if(expected.empty())
    {
        // No manifest entry for this sample means the project was built
        // without this sample, or the sample name the job was given does
        // not match the one the configuration was generated from. Either
        // way nothing downstream can be trusted.
        rep << "ERROR=no_manifest_entries_for_sample\n";
        finish(rep, "manifest_error", 1);
    }

    // -----------------------------------------------------------------
    // Selection output
    // -----------------------------------------------------------------
    TFile * fno = TFile::Open(nosyst_path, "READ");
    if(fno == nullptr || fno->IsZombie())
    {
        rep << "ERROR=nosyst_unopenable\n";
        finish(rep, "output_empty", 2);
    }

    // The expected sample name is known here (it comes from the job's DB
    // row), so require that exact directory rather than inferring it from
    // "the one subdirectory that happens to be present".
    TDirectory * dno = fno->GetDirectory(("events/" + std::string(sample)).c_str());
    if(dno == nullptr)
    {
        rep << "ERROR=nosyst_missing_sample_dir\n";
        finish(rep, "output_empty", 2);
    }

    const double pot_no  = exposure_value(dno, "POT");
    const double live_no = exposure_value(dno, "Livetime");
    rep << "POT=" << pot_no << "\n";
    rep << "LIVETIME=" << live_no << "\n";
    if(std::isnan(pot_no) || std::isnan(live_no))
    {
        rep << "ERROR=nosyst_missing_exposure\n";
        finish(rep, "output_empty", 2);
    }
    // Offbeam/intime samples legitimately carry zero POT and are bookkept
    // by livetime instead, so require only that *some* exposure was
    // recorded rather than POT > 0 specifically.
    if(!(pot_no > 0.0 || live_no > 0.0))
    {
        rep << "ERROR=nosyst_zero_exposure\n";
        finish(rep, "output_empty", 2);
    }

    // Every tree the manifest names must exist in the selection output.
    std::map<std::string, Long64_t> n_nosyst;
    for(const ManifestEntry & e : expected)
    {
        const Long64_t n = tree_entries(dno, e.name);
        if(n < 0)
        {
            rep << "ERROR=nosyst_missing_tree:" << e.name << "\n";
            finish(rep, "output_empty", 2);
        }
        n_nosyst[e.name] = n;
    }

    // -----------------------------------------------------------------
    // Systematics output
    // -----------------------------------------------------------------
    TFile * fsy = TFile::Open(syst_path, "READ");
    if(fsy == nullptr || fsy->IsZombie())
    {
        rep << "ERROR=syst_unopenable\n";
        finish(rep, "sys_error", 3);
    }

    TDirectory * dsy = fsy->GetDirectory(("events/" + std::string(sample)).c_str());
    if(dsy == nullptr)
    {
        // This is what a run that skipped every tree looks like: a valid,
        // openable file with no sample directory in it at all.
        rep << "ERROR=syst_missing_sample_dir\n";
        finish(rep, "match_empty", 3);
    }

    // Exposure is copied verbatim from the selection file by
    // systematics/src/trees.cc, so an exact comparison is the right one --
    // any difference means the two files did not come from the same run.
    const double pot_sy  = exposure_value(dsy, "POT");
    const double live_sy = exposure_value(dsy, "Livetime");
    if(std::isnan(pot_sy) || std::isnan(live_sy))
    {
        rep << "ERROR=syst_missing_exposure\n";
        finish(rep, "sys_error", 4);
    }
    if(pot_sy != pot_no || live_sy != live_no)
    {
        rep << "ERROR=exposure_mismatch\n";
        rep << "SYST_POT=" << pot_sy << "\n";
        rep << "SYST_LIVETIME=" << live_sy << "\n";
        finish(rep, "sys_error", 4);
    }

    // -----------------------------------------------------------------
    // Per-tree accounting
    // -----------------------------------------------------------------
    // fail_status keeps the *first* hard failure's category, which is the
    // most specific description of what went wrong; later errors are still
    // listed individually as ERROR= lines so nothing is hidden.
    int         rc = 0;
    std::string fail_status;
    bool        warn_match = false;
    for(const ManifestEntry & e : expected)
    {
        const Long64_t n_in  = n_nosyst[e.name];
        const Long64_t n_sel = tree_entries(dsy, e.name);
        const Long64_t n_non = tree_entries(dsy, e.name + "_nonmatched");

        // <nosyst>,<selected>,<nonmatched>; -1 means the tree was absent.
        rep << "TREE_" << e.name << "=" << n_in << "," << n_sel << "," << n_non << "\n";

        if(n_sel < 0)
        {
            // The silent-skip case: the tree was owed and is simply not there.
            rep << "ERROR=syst_missing_tree:" << e.name << "\n";
            rc = 5;
            if(fail_status.empty()) fail_status = "match_empty";
            continue;
        }

        if(e.action == "add_weights")
        {
            // Matched and non-matched records partition the input, so they
            // must add up. A missing _nonmatched tree means none were
            // written, i.e. zero unmatched.
            const Long64_t n_non_eff = (n_non < 0) ? 0 : n_non;
            if(n_sel + n_non_eff != n_in)
            {
                rep << "ERROR=match_accounting:" << e.name << "\n";
                rc = 5;
                if(fail_status.empty()) fail_status = "match_unknown";
                continue;
            }
            if(n_in > 0 && n_sel == 0)
            {
                // Everything fell through to _nonmatched: the weights were
                // never actually applied to anything.
                rep << "ERROR=match_empty:" << e.name << "\n";
                rc = 5;
                if(fail_status.empty()) fail_status = "match_empty";
                continue;
            }
            // Each weight table is filled once per matched record.
            for(const std::string & t : e.table_types)
            {
                const std::string tname = e.name + "_" + t + "Tree";
                const Long64_t n_tab = tree_entries(dsy, tname);
                if(n_tab < 0)
                {
                    rep << "ERROR=missing_weight_table:" << tname << "\n";
                    rc = 5;
                    if(fail_status.empty()) fail_status = "match_empty";
                    continue;
                }
                if(n_tab != n_sel)
                {
                    rep << "ERROR=weight_table_entries:" << tname
                        << ":" << n_tab << "!=" << n_sel << "\n";
                    rc = 5;
                    if(fail_status.empty()) fail_status = "match_unknown";
                }
            }
            // Soft signal only: a partial match is often a real physics
            // effect rather than a broken job, so it is reported for later
            // triage but does not block the transfer.
            if(n_in > 0 && n_non_eff > nonmatched_warn_frac * n_in)
            {
                rep << "WARN=match_partial:" << e.name
                    << ":" << n_non_eff << "/" << n_in << "\n";
                warn_match = true;
            }
        }
        else
        {
            // A plain copy must be exactly that.
            if(n_sel != n_in)
            {
                rep << "ERROR=copy_entries:" << e.name
                    << ":" << n_sel << "!=" << n_in << "\n";
                rc = 5;
                if(fail_status.empty()) fail_status = "match_partial";
            }
        }
    }

    fno->Close();
    fsy->Close();

    if(rc == 0) finish(rep, warn_match ? "match_partial" : "ok", 0);
    finish(rep, fail_status.empty() ? "sys_error" : fail_status, rc);
}
