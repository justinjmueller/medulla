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
 * @return An interaction of the specified type with the given ID and
 * particle multiplicity.
 */
template<typename T>
T generate_interaction(int64_t id, int64_t poffset, multiplicity_t mult = {2, 2, 2, 2, 2});

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
 * @return void
 */
void write_event(caf::StandardRecord * rec, int64_t run, int64_t subrun,
                 int64_t event_num, TH1F * pot, TH1F * nevt, TTree * t);