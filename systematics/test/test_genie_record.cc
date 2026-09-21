/**
 * @file test_genie_record.cc
 * @brief Standalone exercise of the GENIE event record copying.
 * @details This test drives sys::GenieRecordWriter directly against real CAF
 * files, without requiring a selected candidate tree or a full run of the
 * systematics executable. It exists because the path that actually writes a
 * record can only be exercised where the GENIE dictionaries are loadable, which
 * is a property of the environment rather than of the build.
 *
 * The test simulates what the main loop does over a TChain: it copies records
 * out of one file, rolls over to a second file (closing the first, to prove the
 * output tree does not depend on the input staying open), and mixes in the two
 * cases that must produce a default-initialized entry rather than a skipped one
 * -- an out-of-range index, and a file with no record tree at all.
 *
 * Build (with the environment from setup.sh):
 *   g++ -std=c++17 -Wall -O2 \
 *       -I../include -I$(root-config --incdir) \
 *       test_genie_record.cc ../src/genie_record.cc \
 *       -o test_genie_record $(root-config --libs)
 *
 * Run:
 *   ./test_genie_record <caf_file> [second_caf_file]
 *
 * @author mueller@fnal.gov
 */
#include <iostream>
#include <string>
#include <vector>

#include "genie_record.h"

#include "TFile.h"
#include "TTree.h"

namespace
{
    int failures = 0;

    void check(bool condition, const std::string & what)
    {
        std::cout << (condition ? "  [ OK ] " : "  [FAIL] ") << what << std::endl;
        if(!condition)
            ++failures;
    }
}

int main(int argc, char ** argv)
{
    if(argc < 2)
    {
        std::cerr << "usage: " << argv[0] << " <caf_file> [second_caf_file]" << std::endl;
        return 2;
    }

    const std::string output_path = "test_genie_record_out.root";

    /**
     * @brief Report whether the records can be written at all.
     * @details This is the same check the systematics code makes before copying
     * anything. If it fails, the rest of the test cannot run, but that is a
     * statement about the environment rather than a defect: without GENIE set
     * up the records are meant to be skipped with a warning.
     */
    TFile * first = TFile::Open(argv[1], "READ");
    if(first == nullptr || first->IsZombie())
    {
        std::cerr << "Could not open " << argv[1] << std::endl;
        return 2;
    }
    TTree * first_tree = (TTree *) first->Get("GenieEvtRecTree");
    std::cout << "Input: " << argv[1] << std::endl;
    std::cout << "  GenieEvtRecTree entries: "
              << (first_tree != nullptr ? first_tree->GetEntries() : -1) << std::endl;

    std::string reason;
    bool writable = sys::records_are_writable(first_tree, reason);
    std::cout << "  records_are_writable: " << (writable ? "true" : "false") << std::endl;
    if(!writable)
    {
        std::cout << "  reason: " << reason << std::endl;
        std::cout << std::endl
                  << "The GENIE dictionaries are not available, so the write path cannot be"
                  << std::endl
                  << "exercised here. This is the guarded path: the systematics code will skip"
                  << std::endl
                  << "the records with a warning rather than crash. Set up GENIE (see setup.sh)"
                  << std::endl
                  << "and run this again to exercise the write path." << std::endl;
        return 1;
    }

    /**
     * @brief Copy a mixture of present and absent records.
     * @details The expected number of entries is tracked alongside, since the
     * whole point of the default-initialized entries is that the record tree
     * stays aligned with the selected candidate tree.
     */
    size_t expected = 0;
    size_t expected_defaults = 0;
    {
        TFile * out = TFile::Open(output_path.c_str(), "RECREATE");
        sys::GenieRecordWriter writer("selected_genieTree", out);

        std::cout << std::endl << "Copying from the first file..." << std::endl;
        const Long64_t n1 = first_tree->GetEntries();
        for(Long64_t idx : {(Long64_t)0, n1 / 2, n1 - 1})
        {
            writer.fill(first_tree, 0, idx);
            ++expected;
        }

        // An index past the end of the record tree must produce a default entry.
        writer.fill(first_tree, 0, n1 + 1000);
        ++expected;
        ++expected_defaults;

        // So must a negative index.
        writer.fill(first_tree, 0, -1);
        ++expected;
        ++expected_defaults;

        // So must a file that carries no record tree at all.
        writer.fill(nullptr, 0, 0);
        ++expected;
        ++expected_defaults;

        /**
         * @brief Roll over to a second file.
         * @details The first file is closed before any record is read from the
         * second, which is what a TChain does. If the writer had kept sharing
         * the input file's buffers this is where it would fall over.
         */
        if(argc > 2)
        {
            TFile * second = TFile::Open(argv[2], "READ");
            TTree * second_tree = (TTree *) second->Get("GenieEvtRecTree");
            std::cout << "Rolling over to " << argv[2] << " ("
                      << (second_tree != nullptr ? second_tree->GetEntries() : -1)
                      << " entries)..." << std::endl;

            first->Close();
            delete first;
            first = nullptr;

            if(second_tree != nullptr)
            {
                const Long64_t n2 = second_tree->GetEntries();
                for(Long64_t idx : {(Long64_t)0, n2 / 2})
                {
                    writer.fill(second_tree, 1, idx);
                    ++expected;
                }
            }

            writer.write();
            check(writer.get_ndefaulted() == expected_defaults,
                  "default-initialized count is " + std::to_string(expected_defaults));
            second->Close();
            delete second;
        }
        else
        {
            writer.write();
            check(writer.get_ndefaulted() == expected_defaults,
                  "default-initialized count is " + std::to_string(expected_defaults));
        }

        out->Close();
        delete out;
    }
    if(first != nullptr)
    {
        first->Close();
        delete first;
    }

    /**
     * @brief Read the output back and confirm it is intact.
     */
    std::cout << std::endl << "Reading back " << output_path << "..." << std::endl;
    TFile * check_file = TFile::Open(output_path.c_str(), "READ");
    check(check_file != nullptr && !check_file->IsZombie(), "output file opens");
    if(check_file == nullptr || check_file->IsZombie())
        return 1;

    TTree * result = (TTree *) check_file->Get("selected_genieTree");
    check(result != nullptr, "output tree 'selected_genieTree' is present");
    if(result != nullptr)
    {
        check((size_t) result->GetEntries() == expected,
              "output has " + std::to_string(expected) + " entries (one per candidate), got "
                  + std::to_string(result->GetEntries()));
        check(result->GetBranch("GenieEvtRec") != nullptr,
              "output carries the 'GenieEvtRec' branch");

        // Reading every entry back exercises the streaming in the other
        // direction, which is where a half-written record would show up.
        Long64_t nread = 0;
        for(Long64_t i(0); i < result->GetEntries(); ++i)
            nread += result->GetEntry(i) > 0 ? 1 : 0;
        check(nread == result->GetEntries(), "every entry reads back");
    }
    check_file->Close();
    delete check_file;

    std::cout << std::endl
              << (failures == 0 ? "ALL CHECKS PASSED" : "FAILURES: " + std::to_string(failures))
              << std::endl;
    return failures == 0 ? 0 : 1;
}
