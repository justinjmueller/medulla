/**
 * @file genie_record.h
 * @brief Header file for the GenieRecordWriter class.
 * @details This class copies the GENIE event records ("GenieEvtRecTree") that
 * belong to the selected signal candidates into a tree in the output file that
 * is entry-aligned with the selected candidate tree. The records are handled
 * generically through the branch structure of the input tree, so that this code
 * never needs to name (or link against) the GENIE classes themselves; ROOT
 * builds emulated classes from the StreamerInfo carried by the CAF files.
 * @author mueller@fnal.gov
 */
#ifndef GENIE_RECORD_H
#define GENIE_RECORD_H
#include <cstdint>
#include <string>
#include <vector>

class TBranch;
class TBranchElement;
class TClass;
class TDirectory;
class TTree;

namespace sys
{
    /**
     * @brief Check that the records held by a tree can actually be written out.
     * @details Copying a record entry-by-entry requires ROOT to stream the
     * object, and the GENIE classes stored in the CAF files all derive from
     * TObject. When no compiled dictionary is loaded for them, ROOT still builds
     * emulated classes from the StreamerInfo carried by the file -- which is
     * enough to read them -- but it also reports them as TObject-derived, and so
     * writing one dispatches a virtual Streamer() call against an object that
     * has no vtable. That is an outright segmentation fault rather than an
     * error, so the dictionaries have to be confirmed before any record is
     * copied. The GENIE dictionaries are loaded automatically by ROOT (through
     * the rootmap files that GENIE ships) whenever GENIE is set up in the
     * environment; note that setting up sbnana alone does not pull GENIE in.
     * @param source The "GenieEvtRecTree" to check.
     * @param reason Set to a description of the first problem found.
     * @return True if the records can be written, false otherwise.
     */
    bool records_are_writable(TTree * source, std::string & reason);

    /**
     * @class GenieRecordWriter
     * @brief A class that copies GENIE event records for selected candidates.
     * @details This class lazily clones the structure of the "GenieEvtRecTree"
     * found in the input CAF files and fills one entry per selected signal
     * candidate, in the same order as the selected candidate tree. Candidates
     * without an associated GENIE event record receive a default-initialized
     * entry so that the entry-by-entry correspondence between the two trees is
     * never broken.
     * @note The GENIE event record index stored in the CAF record
     * ("rec.mc.nu.genie_evtrec_idx") is relative to the file that the record
     * lives in, not to the position within a TChain spanning many files. The
     * caller is therefore expected to supply the "GenieEvtRecTree" belonging to
     * the file that is currently loaded, along with an identifier for that file
     * so that this class can re-point the input branches when the file changes.
     */
    class GenieRecordWriter
    {
        public:

        /**
         * @brief Constructor for the GenieRecordWriter class.
         * @details This constructor stores the name to give the output tree and
         * the directory that it will eventually be written to. No tree is
         * created until the first record is offered by @ref fill(), since the
         * structure of the output tree is cloned from the input.
         * @param name The name to give the output TTree.
         * @param directory The output directory that the TTree belongs to.
         */
        GenieRecordWriter(const std::string & name, TDirectory * directory);

        /**
         * @brief Destructor for the GenieRecordWriter class.
         * @details This destructor releases the record objects that the class
         * allocated to back the branches of the output TTree.
         */
        ~GenieRecordWriter();

        GenieRecordWriter(const GenieRecordWriter &) = delete;
        GenieRecordWriter & operator=(const GenieRecordWriter &) = delete;

        /**
         * @brief Append the GENIE event record at the requested index.
         * @details This method appends one entry to the output TTree. If the
         * requested index does not point at a valid entry of the supplied tree,
         * or if no tree is available for the current file, a default-initialized
         * entry is appended instead. Exactly one entry is appended per call, so
         * that the output tree stays aligned with the selected candidate tree.
         * @param source The "GenieEvtRecTree" of the file that is currently
         * loaded, or nullptr if the current file does not have one.
         * @param file_index An identifier for the currently loaded file, used to
         * detect when the input branches need to be re-pointed. In practice this
         * is the TChain tree number.
         * @param idx The index of the record within the source tree. A negative
         * value, or one beyond the end of the source tree, is treated as "no
         * record present".
         * @return void
         */
        void fill(TTree * source, int file_index, int64_t idx);

        /**
         * @brief Write the output TTree to the configured directory.
         * @details This method writes the output TTree to the directory that was
         * supplied to the constructor. If no entry was ever offered, no tree is
         * created and this method does nothing.
         * @return void
         */
        void write();

        /**
         * @brief Get the number of default-initialized entries written.
         * @details This method returns the number of entries that were written
         * with a default-initialized record because no GENIE event record could
         * be located for the corresponding selected candidate.
         * @return The number of default-initialized entries.
         */
        size_t get_ndefaulted() const { return ndefaulted; }

        private:

        /**
         * @brief A single top-level branch of the output TTree.
         * @details This struct holds the storage backing one top-level branch of
         * the output TTree. Object-valued branches (such as the GENIE event
         * record itself) are backed by an object allocated through the branch's
         * TClass, alongside a second, pristine object that is swapped in to
         * produce a default-initialized entry. Branches of fundamental type
         * (such as the "GENIEEntry" and "SourceFileHash" branches found in
         * structured CAF files) are backed by a local byte buffer instead.
         */
        struct Slot
        {
            TBranch * branch{nullptr}; ///< The output branch.
            TBranchElement * element{nullptr}; ///< The output branch, if object-valued.
            TClass * cls{nullptr}; ///< The class of the branch, if object-valued.
            void * object{nullptr}; ///< The object the branch reads and writes.
            void * blank{nullptr}; ///< A pristine, default-constructed object.
            std::vector<char> buffer; ///< Storage for a fundamental-type branch.
        };

        /**
         * @brief Create the output TTree and/or re-point the input branches.
         * @details This method clones the structure of the source tree on the
         * first call and takes ownership of the storage backing its branches. On
         * every call where the currently loaded file has changed, the branches of
         * the source tree are pointed at that storage so that reading an entry
         * from the source populates the output branches.
         * @param source The "GenieEvtRecTree" of the currently loaded file.
         * @param file_index An identifier for the currently loaded file.
         * @return void
         */
        void connect(TTree * source, int file_index);

        /**
         * @brief Take ownership of the storage backing the output branches.
         * @details TTree::CloneTree() leaves the cloned branches sharing the
         * buffers of the input tree, which are destroyed as soon as a TChain
         * moves on to the next file. This method replaces those buffers with
         * storage owned by this class so that the output tree survives the input
         * file rotating, and allocates the pristine objects used to produce
         * default-initialized entries.
         * @return void
         */
        void take_ownership();

        /**
         * @brief Append a default-initialized entry to the output TTree.
         * @details This method temporarily swaps the pristine objects into the
         * output branches, appends an entry, and then restores the objects that
         * the input tree reads into.
         * @return void
         */
        void fill_default();

        std::string name; ///< The name to give the output TTree.
        TDirectory * directory{nullptr}; ///< The output directory.
        TTree * output{nullptr}; ///< The output TTree.
        int connected{-1}; ///< The file the input branches are currently pointed at.
        std::vector<Slot> slots; ///< The storage backing the output branches.
        size_t pending{0}; ///< Default entries owed from before the tree existed.
        size_t nfilled{0}; ///< The number of entries that have been requested.
        size_t ndefaulted{0}; ///< The number of default-initialized entries.
    };
} // namespace sys
#endif // GENIE_RECORD_H
