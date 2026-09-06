/**
 * @file genie_record.cc
 * @brief Implementation file for the GenieRecordWriter class.
 * @details This file contains the implementation of the GenieRecordWriter
 * class, which copies the GENIE event records belonging to the selected signal
 * candidates into a tree in the output file that is entry-aligned with the
 * selected candidate tree.
 * @author mueller@fnal.gov
 */
#include <cstring>
#include <iostream>

#include "genie_record.h"

#include "TBranch.h"
#include "TBranchElement.h"
#include "TClass.h"
#include "TDirectory.h"
#include "THashTable.h"
#include "TLeaf.h"
#include "TObjArray.h"
#include "TSystem.h"
#include "TTree.h"

namespace
{
    /**
     * @brief The outcome of trying to load the GENIE dictionaries.
     */
    enum class GenieLoad { kNotFound, kFailed, kOk };

    std::string failed_library; ///< The library that could not be loaded, if any.

    /**
     * @brief Try to make the GENIE dictionaries available, once per job.
     * @details ROOT autoloads the record classes through the rootmap files that
     * GENIE ships, but autoloading pulls in only the one library that declares
     * them, and that is not enough to use them.
     *
     * GENIE's framework libraries are mutually under-linked: every one of them
     * references symbols defined in the others without naming those libraries in
     * its own DT_NEEDED, on the assumption that the application links the whole
     * set together. Loading only the library that declares the record classes
     * therefore leaves symbols unresolved, and because binding is lazy that does
     * not fail at load time -- it fails at the first call, which is partway
     * through a job, on the first candidate that actually has a record. The
     * whole framework set is loaded here so that every symbol is resident before
     * any record is touched. Loading them by soname rather than by path keeps
     * this independent of the GENIE version that happens to be set up.
     * @return Whether the libraries were absent, present but unloadable, or
     * loaded.
     */
    GenieLoad try_load_genie()
    {
        static bool attempted = false;
        static GenieLoad status = GenieLoad::kNotFound;
        if(attempted)
            return status;
        attempted = true;

        /**
         * @note Ordered so that the more depended-upon libraries come first.
         * The order is not strictly required, since ROOT loads with RTLD_GLOBAL
         * and the dependencies between these libraries are circular in any case,
         * but it keeps the sequence readable.
         */
        static const char * const libraries[] = {
            "libGFwUtl", "libGFwMsg", "libGFwReg", "libGFwNum", "libGFwAlg",
            "libGFwParDat", "libGFwInt", "libGFwGHEP", "libGFwEG", "libGFwNtp"
        };

        // The record classes live in libGFwNtp. If that is not on the path then
        // GENIE is simply not set up, which is reported by the dictionary check
        // rather than as a load failure. Probing first also keeps ROOT from
        // printing errors of its own in that entirely expected situation.
        char * probe = gSystem->DynamicPathName("libGFwNtp", kTRUE);
        if(probe == nullptr)
            return status;
        delete [] probe;

        /**
         * @note A failed load has to be treated as fatal to the feature in its
         * own right, rather than relying on the class-level checks below. When a
         * dictionary library is found but cannot be loaded -- in practice
         * because it was built against a different ROOT -- the interpreter has
         * still parsed the dictionary's payload headers by that point, and the
         * class then reports a dictionary, reports itself as loaded, and reports
         * nothing missing, none of which is backed by compiled code.
         */
        status = GenieLoad::kOk;
        for(const char * library : libraries)
        {
            if(gSystem->Load(library) < 0)
            {
                status = GenieLoad::kFailed;
                if(failed_library.empty())
                    failed_library = library;
            }
        }
        return status;
    }

    /**
     * @brief Collect the classes in a hierarchy that have no dictionary.
     * @details The check is recursive because it is not the record class itself
     * that trips the fault so much as the classes nested inside it (the GENIE
     * particle records held in a TClonesArray, in particular).
     * @param cls The class to check.
     * @param first Set to the name of one of the classes that is missing.
     * @param count Set to the number of classes that are missing.
     * @return True if any dictionary is missing, false otherwise.
     */
    bool has_missing_dictionaries(TClass * cls, std::string & first, int & count)
    {
        THashTable missing;
        cls->GetMissingDictionaries(missing, true);
        count = missing.GetEntries();
        if(count == 0)
            return false;
        TIter next(&missing);
        TObject * entry = next();
        first = entry != nullptr ? entry->GetName() : "?";
        return true;
    }
}

// Check that the records held by a tree can actually be written out.
bool sys::records_are_writable(TTree * source, std::string & reason)
{
    if(source == nullptr)
    {
        reason = "no GenieEvtRecTree is present in the input file";
        return false;
    }

    // Attempt this before any class is looked up, so that the dictionaries are
    // in place before ROOT has a chance to build emulated classes instead.
    if(try_load_genie() == GenieLoad::kFailed)
    {
        reason = "the GENIE library '" + failed_library + "' was found but could not be loaded, "
            "which normally means it was built against a different ROOT version than the one "
            "in use. Set up a GENIE built against the same ROOT as sbnana";
        return false;
    }

    TObjArray * branches = source->GetListOfBranches();
    for(Int_t i(0); i < branches->GetEntriesFast(); ++i)
    {
        TBranchElement * element = dynamic_cast<TBranchElement *>(branches->UncheckedAt(i));
        if(element == nullptr)
            continue;

        TClass * cls = TClass::GetClass(element->GetClassName());
        if(cls == nullptr)
        {
            reason = std::string("no class is known for branch '") + element->GetName() + "'";
            return false;
        }

        std::string first;
        int count(0);
        if(has_missing_dictionaries(cls, first, count))
        {
            reason = std::string("no dictionary is loaded for ") + cls->GetName()
                + " (" + std::to_string(count) + " class(es) missing, e.g. " + first
                + "). Set up GENIE in the environment so that ROOT can autoload its "
                "dictionaries, matching the GENIE_VERSION recorded in the CAF file's "
                "'env/envtree' (e.g. 'setup genie v3_04_02a -q e26:prof')";
            return false;
        }
    }

    return true;
}

// Constructor for the GenieRecordWriter class.
sys::GenieRecordWriter::GenieRecordWriter(const std::string & n, TDirectory * d)
: name(n),
  directory(d)
{
}

// Destructor for the GenieRecordWriter class.
sys::GenieRecordWriter::~GenieRecordWriter()
{
    for(Slot & slot : slots)
    {
        if(slot.cls == nullptr)
            continue;
        if(slot.object != nullptr)
            slot.cls->Destructor(slot.object);
        if(slot.blank != nullptr)
            slot.cls->Destructor(slot.blank);
    }
}

// Append the GENIE event record at the requested index.
void sys::GenieRecordWriter::fill(TTree * source, int file_index, int64_t idx)
{
    ++nfilled;

    /**
     * @brief Handle the case where the currently loaded file has no GENIE
     * event record tree at all.
     * @details If the output tree does not exist yet there is no structure to
     * clone, so the default entry is deferred until a file that does have the
     * tree is reached. This keeps the output tree aligned with the selected
     * candidate tree even if the very first files of the input are missing the
     * GENIE event records.
     */
    if(source == nullptr)
    {
        if(output == nullptr)
            ++pending;
        else
            fill_default();
        ++ndefaulted;
        return;
    }

    connect(source, file_index);

    /**
     * @brief Copy the requested record, or fall back to a default entry.
     * @details The index carried by the CAF record is relative to the file that
     * the record lives in, so it is validated against the source tree belonging
     * to the currently loaded file. Note that the CAF record initializes the
     * index to zero rather than to a distinguishable sentinel, so a neutrino
     * with no GENIE event record cannot be told apart from one that legitimately
     * points at the first entry. Only an out-of-range index can be caught here.
     */
    if(idx >= 0 && idx < source->GetEntries())
    {
        source->GetEntry(idx);
        output->Fill();
    }
    else
    {
        fill_default();
        ++ndefaulted;
    }
}

// Create the output TTree and/or re-point the input branches.
void sys::GenieRecordWriter::connect(TTree * source, int file_index)
{
    if(output == nullptr)
    {
        /**
         * @brief Clone the structure, but none of the entries, of the input
         * tree.
         * @details Cloning the structure rather than declaring the branches by
         * hand means that this code never has to name the GENIE classes, and
         * therefore never has to link against GENIE: ROOT builds emulated
         * classes from the StreamerInfo carried by the CAF files. It also means
         * that any auxiliary branches sitting alongside the record itself (the
         * "GENIEEntry" and "SourceFileHash" branches of structured CAF files,
         * for example) are carried over without further effort.
         */
        /**
         * @note TTree::CloneTree() attaches the clone to the current ROOT
         * directory, which at this point is one of the input CAF files rather
         * than the output file. The clone is therefore built while the output
         * directory is made current, and the previous directory is restored
         * afterwards so that nothing else in the caller is disturbed.
         */
        TDirectory * previous = gDirectory;
        directory->cd();
        output = source->CloneTree(0);
        output->SetName(name.c_str());
        output->SetTitle(name.c_str());
        output->SetDirectory(directory);
        output->SetAutoFlush(1000);
        take_ownership();
        if(previous != nullptr)
            previous->cd();

        // Emit the default entries that were owed from before the tree existed.
        for(size_t i(0); i < pending; ++i)
            fill_default();
        pending = 0;

        // Force the input branches to be pointed at our storage below.
        connected = -1;
    }

    /**
     * @brief Point the branches of the input tree at our own storage.
     * @details This has to be redone every time the TChain moves on to a new
     * file, since each file carries its own instance of the GENIE event record
     * tree. The comparison is made on the file index rather than on the tree
     * pointer because a newly opened file may well reuse the address of the tree
     * that was just deleted.
     */
    if(file_index != connected)
    {
        output->CopyAddresses(source);
        connected = file_index;
    }
}

// Take ownership of the storage backing the output branches.
void sys::GenieRecordWriter::take_ownership()
{
    /**
     * @brief Collect the top-level branches of the output tree.
     * @details TTree::CloneTree() leaves the cloned branches sharing the buffers
     * of the input tree. Those buffers belong to a tree in the currently open
     * input file and are destroyed as soon as the TChain moves on, so the output
     * tree would be left holding dangling pointers. The storage is therefore
     * replaced with storage owned by this class.
     */
    TObjArray * branches = output->GetListOfBranches();
    for(Int_t i(0); i < branches->GetEntriesFast(); ++i)
    {
        Slot slot;
        slot.branch = (TBranch *) branches->UncheckedAt(i);
        slot.element = dynamic_cast<TBranchElement *>(slot.branch);
        if(slot.element != nullptr)
        {
            slot.cls = slot.element->GetTargetClass();
            if(slot.cls == nullptr)
                slot.cls = TClass::GetClass(slot.element->GetClassName());
        }

        if(slot.cls == nullptr)
        {
            /**
             * @brief Mirror the storage of a branch of fundamental type.
             * @details Branches that do not hold an object (the "GENIEEntry" and
             * "SourceFileHash" branches, for example) are backed by a local byte
             * buffer sized from the leaves of the branch.
             */
            size_t nbytes(0);
            TObjArray * leaves = slot.branch->GetListOfLeaves();
            for(Int_t j(0); j < leaves->GetEntriesFast(); ++j)
            {
                TLeaf * leaf = (TLeaf *) leaves->UncheckedAt(j);
                nbytes += (size_t) leaf->GetLenType() * (size_t) leaf->GetLen();
            }
            slot.buffer.assign(nbytes > 0 ? nbytes : sizeof(Long64_t), 0);
        }

        slots.push_back(std::move(slot));
    }

    /**
     * @brief Attach the storage to the output branches.
     * @details This is done in a second pass so that the addresses handed to the
     * branches are those of the buffers as they sit in the "slots" vector, which
     * is not reallocated any further.
     */
    for(Slot & slot : slots)
    {
        if(slot.cls != nullptr)
        {
            slot.object = slot.cls->New();
            slot.blank = slot.cls->New();
            slot.element->SetObject(slot.object);
        }
        else
            slot.branch->SetAddress(slot.buffer.data());
    }
}

// Append a default-initialized entry to the output TTree.
void sys::GenieRecordWriter::fill_default()
{
    /**
     * @brief Swap in the pristine objects, append an entry, and swap back.
     * @details The pristine objects are never written into by the input tree, so
     * they remain in the state left by their default constructor for the whole
     * life of this class. Swapping them in is preferred over clearing the
     * objects that the input tree reads into, which would require destroying and
     * reconstructing a live buffer in place.
     */
    for(Slot & slot : slots)
    {
        if(slot.cls != nullptr)
            slot.element->SetObject(slot.blank);
        else
            std::memset(slot.buffer.data(), 0, slot.buffer.size());
    }

    output->Fill();

    for(Slot & slot : slots)
    {
        if(slot.cls != nullptr)
            slot.element->SetObject(slot.object);
    }
}

// Write the output TTree to the configured directory.
void sys::GenieRecordWriter::write()
{
    if(output == nullptr)
    {
        if(nfilled > 0)
            std::cerr << "Warning: No GENIE event record tree was found in any input file. "
                      << "No GENIE event records were written." << std::endl;
        return;
    }

    /**
     * @brief Check that the output tree stayed aligned with the selected
     * candidate tree.
     * @details The entire point of the default-initialized entries is that the
     * two trees can be read entry-by-entry in lockstep. A mismatch here would
     * silently associate records with the wrong candidates downstream, so it is
     * worth reporting loudly.
     */
    if((size_t) output->GetEntries() != nfilled)
        std::cerr << "Warning: The GENIE event record tree '" << name << "' has "
                  << output->GetEntries() << " entries but " << nfilled
                  << " selected candidates were processed. The trees are NOT aligned."
                  << std::endl;

    std::cout << "Wrote " << output->GetEntries() << " GENIE event records to '"
              << name << "' (" << ndefaulted << " default-initialized)." << std::endl;

    directory->WriteObject(output, name.c_str());
}
