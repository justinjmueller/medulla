/**
 * @file test.h
 * @brief Header file for the framework test library.
 * @details This file contains the declarations of functions and types used for
 * testing the framework functionality, including interaction- and particle-
 * level variables, matching between truth and reconstructed objects, and all
 * manner of cases surrounding these functionalities. This is, essentially,
 * intended to be a "brute force" test of the framework, ensuring that all
 * functionality works as expected. The alternative is a full "logic" test
 * of the framework, which does tend to be more error-prone due to the 
 * variety of possible cases and the complexity of the logic involved.
 * @author mueller@fnal.gov
 */
#include <vector>
#include <array>

#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"
#include "sbnanaobj/StandardRecord/SRParticleDLP.h"
#include "sbnanaobj/StandardRecord/SRParticleTruthDLP.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

// Define a scale for the energy values used in the test. This is basically
// a "nuance offset" to ensure that the values are not zero, which can result
// in particles falling below various thresholds in the framework.
#define ENERGY_SCALE 200.0

// Define a type for multiplicity, which is used to define the number of
// particles of each type in an interaction. This is a fixed-size array
// of integers, where each element corresponds to a different particle type.
typedef std::array<int64_t, 5> multiplicity_t;

// Define a typedef for a row of information.
typedef std::map<std::string, double> row_t;

// Defined a typedef for a "condition," or a named row of information that
// defines a test condition.
typedef std::pair<std::string, row_t> condition_t;

// Simple NaN value
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

/**
 * @brief Collect the offset for a given object.
 * @tparam T The type of the object for which to collect the offset.
 * @param obj The object for which to collect the offset.
 * @return The offset for the object, which is the number of particles
 * in the interaction or the number of particles in the vector of interactions.
 * This is intended to be used to ensure continuity in the particle IDs
 * across interactions, so that the IDs are unique and do not overlap.
 * @return size_t The offset for the object.
 */
template<typename T>
size_t collect_offset(const T & obj);

/**
 * @brief Generate a particle of the specified type.
 * @tparam T The type of the particle to generate.
 * @param id The ID of the particle.
 * @param pid The PID of the particle.
 * @return A particle of the specified type with the given ID and PID.
 * This function is used to generate particles for testing purposes,
 * ensuring that the particles have unique IDs and PIDs.
 * @return T A particle of the specified type.
 */
template<typename T>
T generate_particle(int64_t id, int64_t pid);

/**
 * @brief Generate an interaction of the specified type with the requested
 * particle multiplicity.
 * @tparam T The type of the interaction to generate.
 * @param id The ID of the interaction.
 * @param poffset The offset for the particles in the interaction.
 * @param mult The multiplicity of particles in the interaction, which is a
 * fixed-size array of integers representing the number of particles of each
 * type in the interaction.
 * @param assign_fm Whether to assign a flash match to the interaction.
 * @return An interaction of the specified type with the given ID and
 * particle multiplicity.
 */
template<typename T>
T generate_interaction(int64_t id,
                       int64_t poffset,
                       multiplicity_t mult = {2, 2, 2, 2, 2},
                       bool assign_fm = true);

/**
 * @brief Pair two interactions together, setting the match IDs.
 * @tparam T The type of the interaction (left) to pair.
 * @tparam U The type of the interaction (right) to pair.
 * @param left The interaction to pair.
 * @param right The true interaction to pair.
 * This function is used to create a match between the left and right
 * interactions, setting the match IDs accordingly.
 * @return void
 */
template<typename T, typename U>
void pair(T & left, U & right);

/**
 * @brief Mark the particles as contained.
 * @details This function marks all particles in the interaction as
 * contained, which is used to create a set of "passing" particles (w.r.t.
 * containment) in the interaction.
 * @param reco_interaction The interaction to mark the particles as contained.
 * @param true_interaction The true interaction to mark the particles as contained.
 * @return void
 * @note The true interaction may be null, in which case only the
 * reco_interaction is marked as contained.
 */
void mark_contained(caf::SRInteractionDLP * reco_interaction,
                    caf::SRInteractionTruthDLP * true_interaction = nullptr);

/**
 * @brief Set the event metadata in the StandardRecord header.
 * @details This function sets the run, subrun, and event number for the
 * StandardRecord header in the test, along with the exposure information
 * and the size of the interaction vectors. Finally, this fills the exposure
 * histograms and fills the tree. This is a one-liner shortcut to improve
 * readability and reduce boilerplate code in the test code.
 * @param rec The StandardRecord pointer to set the event numbers in.
 * @param run The run number to set.
 * @param subrun The subrun number to set.
 * @param event_num The event number to set.
 * @param pot The histogram for the total POT (optional).
 * @param nevt The histogram for the total number of events (optional).
 * @param t The TTree to fill with the event data (optional).
 * @param trigger_time The global trigger time to set in the header (default is
 * 500).
 * @return void
 */
void write_event(caf::StandardRecord * rec, int64_t run, int64_t subrun,
                 int64_t event_num, TH1F * pot, TH1F * nevt, TTree * t,
                 int32_t trigger_time = 500);

/**
 * @brief Read the event data from the TTree at the specified path.
 * @details This function reads the event data from the TTree at the specified
 * path, extracting the run, subrun, event number, and all associated
 * additional variables.
 * @param name The name of the TTree to read from.
 * @return std::vector<row_t> A vector of rows, where each row is a map
 * containing the event data.
 */
std::vector<row_t> read_event_data(const std::string & name);

/**
 * @brief Match the (Run, Subrun, Evt) metadata between a row_t object and a
 * condition_t object.
 * @details This function checks if the (Run, Subrun, Evt) metadata in two
 * row_t objects match. This is used to check if the condition_t being
 * validated against the row is at least from the correct event for the
 * condition.
 * @param row The row_t object to check.
 * @param condition The condition_t object to check against.
 * @return bool True if the (Run, Subrun, Evt) metadata match, false otherwise.
 */
bool match_metadata(const row_t & row, const condition_t & condition);

/**
 * @brief Match a set of row_t objects against a set of condition_t objects.
 * @details Check if each condition_t entry is met by the corresponding (via
 * Run, Subrun, Evt metadata) row_t object. This is used to validate the
 * results of the test against the expected results.
 * @param rows The vector of row_t objects to check.
 * @param conditions The vector of condition_t objects to check against.
 */
void match_conditions(const std::vector<row_t> & rows,
                      const std::vector<condition_t> & conditions);