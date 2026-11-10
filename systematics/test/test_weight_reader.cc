/**
 * @file test_weight_reader.cc
 * @brief Standalone check that sys::WeightReader visits every entry of its
 * chain, in order, exactly once.
 * @details This is a regression test for a silent data-loss bug. The
 * constructor leaves the TTreeReader and the TChain positioned on entry 0, but
 * next() used to advance before returning, so the caller's first iteration saw
 * entry 1 and entry 0 was never processed. The selected candidates belonging to
 * that event therefore matched no weights and were written to the
 * "_nonmatched" tree instead -- about 0.26% of selected neutrino candidates,
 * one event per job. Nothing caught it: the systematics executable still
 * returned 0, and the matched/non-matched counts still added up, because the
 * lost candidates were accounted for in "_nonmatched".
 *
 * The test drives WeightReader over real CAF files and compares the sequence of
 * (run, subrun, event) it yields against an independent pass over the same
 * chain made with a plain TTreeReader. Using a separate reader as the oracle,
 * rather than re-deriving the expectation from WeightReader itself, is what
 * makes the off-by-one detectable.
 *
 * It also checks the per-neutrino content of each entry, which is a distinct
 * failure mode: the event identity comes from the TTreeReader, while the
 * per-neutrino arrays come from chain.GetEntry() via SetBranchAddress. If those
 * buffers do not hold the entry the reader claims to be on -- for instance if
 * they are unfilled on the first iteration -- then get_nnu() returns 0 there,
 * the caller's `for(idn = 0; idn < get_nnu(); ++idn)` loop body never runs, and
 * that entry's candidates are silently left unmatched. The output is then
 * indistinguishable from the entry having been skipped outright, so checking
 * the event identity alone cannot tell the two apart.
 *
 * Build (with the environment from setup.sh):
 *   g++ -std=c++17 -Wall -O2 \
 *       -I../include -I$SBNANAOBJ_INC -I$SRPROXY_INC -I$(root-config --incdir) \
 *       test_weight_reader.cc ../src/weight_reader.cc \
 *       -o test_weight_reader \
 *       $(root-config --libs) -L$SBNANAOBJ_LIB \
 *       -lsbnanaobj_StandardRecord -lsbnanaobj_StandardRecordProxy
 *
 * Run (the argument takes the same forms WeightReader accepts: one file, a
 * quoted wildcard, or a .txt list of paths):
 *   ./test_weight_reader /pnfs/.../caf-....flat.caf.root
 *   ./test_weight_reader "*flat*.root"        (quote it: the shell must not expand it)
 *
 * Exits 0 if every check passes, 1 otherwise.
 *
 * @author mueller@fnal.gov
 */
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "weight_reader.h"

namespace
{
    using event_t = std::tuple<uint32_t, uint32_t, uint32_t>;

    int failures = 0;

    void check(bool ok, const std::string & what)
    {
        std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << std::endl;
        if(!ok) ++failures;
    }

    /**
     * @brief Populate a TChain the same way WeightReader's constructor does.
     * @details Kept deliberately parallel to that logic so the oracle reads
     * exactly the same files in the same order; everything after this point is
     * independent of WeightReader.
     */
    void add_input(TChain & chain, const std::string & input)
    {
        if(input.find("*") != std::string::npos)
        {
            chain.Add(input.c_str());
        }
        else if(input.find(".txt") != std::string::npos)
        {
            std::ifstream infile(input);
            std::string line;
            while(std::getline(infile, line))
                if(!line.empty()) chain.Add(line.c_str());
        }
        else
        {
            chain.Add(input.c_str());
        }
    }

    /// The (run, subrun, event) of every entry, read independently.
    std::vector<event_t> expected_events(const std::string & input)
    {
        TChain chain("recTree");
        add_input(chain, input);
        TTreeReader reader(&chain);
        TTreeReaderValue<uint32_t> run(reader, "rec.hdr.run");
        TTreeReaderValue<uint32_t> subrun(reader, "rec.hdr.subrun");
        TTreeReaderValue<uint32_t> event(reader, "rec.hdr.evt");
        std::vector<event_t> out;
        while(reader.Next())
            out.emplace_back(*run, *subrun, *event);
        return out;
    }

    /// Number of files the input expands to.
    int n_files(const std::string & input)
    {
        TChain chain("recTree");
        add_input(chain, input);
        return chain.GetNtrees();
    }

    /**
     * @brief The neutrino count of every entry, read independently.
     * @details Flat CAFs only ("rec.mc.nu..length"); returns an empty vector for
     * a structured CAF, and the content checks are then skipped. Read with
     * SetBranchAddress on its own chain, so nothing here shares state with the
     * WeightReader under test.
     */
    std::vector<int> expected_nnu(const std::string & input)
    {
        TChain chain("recTree");
        add_input(chain, input);
        if(chain.GetBranch("rec.mc.nu..length") == nullptr) return {};
        Int_t nnu(0);
        chain.SetBranchAddress("rec.mc.nu..length", &nnu);
        std::vector<int> out;
        const Long64_t n = chain.GetEntries();
        for(Long64_t i(0); i < n; ++i)
        {
            nnu = -1;
            chain.GetEntry(i);
            out.push_back(nnu);
        }
        return out;
    }
}

int main(int argc, char * argv[])
{
    if(argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <caf_file | \"pattern\" | list.txt>" << std::endl;
        return 1;
    }
    const std::string input(argv[1]);

    const std::vector<event_t> expected = expected_events(input);
    std::cout << "Oracle: " << expected.size() << " entries in the chain." << std::endl;
    if(expected.empty())
    {
        std::cerr << "Error: the input chain has no entries; nothing to test." << std::endl;
        return 1;
    }

    // A single-file chain does not exercise the case this test exists to
    // cover. TChain::GetEntries() on a multi-file chain built from a wildcard
    // walks to the end, loading the last file and rebinding branch addresses;
    // with one file there is nothing to walk, so a reader that is disturbed by
    // that walk still passes. Always run this against a realistic file set.
    const int nfiles = n_files(input);
    std::cout << "Chain holds " << nfiles << " file(s)." << std::endl;
    if(nfiles < 2)
        std::cout << "  NOTE  single-file chain: the multi-file path is NOT exercised."
                     " Re-run with a wildcard over a realistic set (e.g. a job's 25 files)."
                  << std::endl;

    const std::vector<int> nnu_expected = expected_nnu(input);
    const bool check_content = !nnu_expected.empty();
    std::cout << (check_content ? "Oracle: per-entry neutrino counts available (flat CAF)."
                                : "Oracle: no 'rec.mc.nu..length' branch; content checks skipped.")
              << std::endl;

    std::vector<event_t> seen;
    std::vector<int> nnu_seen;
    {
        sys::WeightReader reader(input);
        while(reader.next())
        {
            seen.emplace_back(reader.get_run(), reader.get_subrun(), reader.get_event());
            nnu_seen.push_back((int)reader.get_nnu());
        }
    }
    std::cout << "WeightReader visited " << seen.size() << " entries." << std::endl;

    check(seen.size() == expected.size(),
          "visits every entry (" + std::to_string(seen.size()) + " of "
              + std::to_string(expected.size()) + ")");

    // The regression itself: the first entry of the chain must be processed.
    if(!seen.empty())
    {
        const bool first_ok = seen.front() == expected.front();
        check(first_ok, "first entry visited is entry 0 of the chain");
        if(!first_ok)
        {
            const auto & e = expected.front();
            const auto & s = seen.front();
            std::cout << "        expected (" << std::get<0>(e) << "," << std::get<1>(e) << ","
                      << std::get<2>(e) << ") but got (" << std::get<0>(s) << "," << std::get<1>(s)
                      << "," << std::get<2>(s) << ")" << std::endl;
            if(seen.size() + 1 == expected.size() && seen.front() == expected[1])
                std::cout << "        -- the whole sequence is shifted by one: entry 0 was skipped."
                          << std::endl;
        }
        check(seen.back() == expected.back(), "last entry visited is the last entry of the chain");
    }

    // Order matters: the weight lookup keys on (run, subrun, event) read from
    // the same entry the per-neutrino buffers were filled from, so a reader
    // that visited the right set of entries in the wrong order would still
    // mismatch every candidate.
    size_t first_diff = seen.size();
    for(size_t i(0); i < std::min(seen.size(), expected.size()); ++i)
        if(seen[i] != expected[i]) { first_diff = i; break; }
    check(first_diff >= std::min(seen.size(), expected.size()),
          "sequence matches the chain entry by entry"
              + (first_diff < std::min(seen.size(), expected.size())
                     ? " (first difference at " + std::to_string(first_diff) + ")"
                     : ""));

    // Content: the per-neutrino buffers must belong to the entry the reader
    // says it is on. get_nnu() == 0 where the file says otherwise means the
    // caller's per-neutrino loop never runs for that entry, losing its
    // candidates exactly as if the entry had been skipped.
    if(check_content)
    {
        size_t wrong = 0, zero_where_nonzero = 0, first_bad = nnu_seen.size();
        for(size_t i(0); i < std::min(nnu_seen.size(), nnu_expected.size()); ++i)
        {
            if(nnu_seen[i] == nnu_expected[i]) continue;
            ++wrong;
            if(nnu_seen[i] == 0 && nnu_expected[i] > 0) ++zero_where_nonzero;
            if(first_bad == nnu_seen.size()) first_bad = i;
        }
        check(wrong == 0, "get_nnu() matches the file for every entry ("
                              + std::to_string(wrong) + " wrong)");
        if(wrong)
        {
            std::cout << "        first mismatch at iteration " << first_bad << ": reader says "
                      << nnu_seen[first_bad] << ", file says " << nnu_expected[first_bad] << std::endl;
            if(zero_where_nonzero)
                std::cout << "        " << zero_where_nonzero
                          << " entr(y|ies) where the reader reports 0 neutrinos but the file has some:"
                             " those entries' candidates cannot be matched." << std::endl;
        }
    }

    std::set<event_t> unique(seen.begin(), seen.end());
    if(unique.size() != seen.size())
        std::cout << "  NOTE  " << seen.size() - unique.size()
                  << " repeated (run,subrun,event) among the visited entries; this is a property"
                     " of the input files, not of the reader." << std::endl;

    std::cout << (failures == 0 ? "ALL CHECKS PASSED" : std::to_string(failures) + " CHECK(S) FAILED")
              << std::endl;
    return failures == 0 ? 0 : 1;
}
