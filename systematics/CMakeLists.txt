cmake_minimum_required(VERSION 3.12)
include(FetchContent)
project(systematics LANGUAGES CXX)

# -O3 rather than -Ofast: -Ofast implies -ffast-math, and therefore
# -ffinite-math-only, which lets the compiler assume no value is ever NaN. That
# silently deletes NaN handling. Concretely, it compiled away the
# std::isnan(nu_id) check in trees.cc that routes cosmic-overlay interactions
# (true_neutrino_id = NaN, the kNoMatchValue convention) to the non-matched
# tree, and simultaneously made "nu_id >= 0" true for those entries so they
# were offered to the candidate map, where a NaN key compares equivalent to a
# real one and the insert was dropped. The candidate then appeared in neither
# output tree: silent data loss, at ~32% of jobs in samples with many
# cosmic-overlay candidates. Any float-sensitive check here needs NaNs to keep
# behaving like NaNs, so the unsafe math is not worth the speed.
add_compile_options(-Wall -Werror -Wno-error=maybe-uninitialized -O3 -g)

# Find packages
find_package(ROOT REQUIRED)
find_package(sbnanaobj)


add_library(weight_reader SHARED src/weight_reader.cc)
target_link_directories(weight_reader PRIVATE ${SBNANAOBJ_LIB})
target_link_libraries(weight_reader PRIVATE ${ROOT_LIBRARIES} shared sbnanaobj_StandardRecordProxy)
target_include_directories(weight_reader PRIVATE include/ ${ROOT_INCLUDE_DIRS} ${SBNANAOBJ_INC} ${SRPROXY_INC})

# Library for generic systematics
add_library(systematics SHARED src/systematic.cc)
target_link_libraries(systematics PRIVATE ${ROOT_LIBRARIES} shared)
target_include_directories(systematics PRIVATE include/ ${ROOT_INCLUDE_DIRS} ${SBNANAOBJ_INC})

# Library for detector systematics
add_library(detsys SHARED src/detsys.cc)
target_link_libraries(detsys PRIVATE ${ROOT_LIBRARIES} shared)
target_include_directories(detsys PRIVATE include/ ${ROOT_INCLUDE_DIRS} ${SBNANAOBJ_INC})

# Library for copying GENIE event records
add_library(genie_record SHARED src/genie_record.cc)
target_link_libraries(genie_record PRIVATE ${ROOT_LIBRARIES} shared)
target_include_directories(genie_record PRIVATE include/ ${ROOT_INCLUDE_DIRS})

# Library for tree-handling code
add_library(trees SHARED src/trees.cc)
target_link_libraries(trees PRIVATE ${ROOT_LIBRARIES} shared weight_reader systematics detsys genie_record)
target_include_directories(trees PRIVATE include/ ${ROOT_INCLUDE_DIRS} ${SBNANAOBJ_INC})

# Library for CAFAna
add_library(sbnanaobj_standardrecord SHARED IMPORTED)
set_target_properties(sbnanaobj_standardrecord PROPERTIES
    IMPORTED_LOCATION "${SBNANAOBJ_LIB}/libsbnanaobj_StandardRecord.so"
    INTERFACE_INCLUDE_DIRECTORIES "${SBNANAOBJ_INC}"
)

# Add the executable target
add_executable(run_systematics src/main.cc)

# Link the ROOT libraries to the target
target_link_libraries(run_systematics PRIVATE ${ROOT_LIBRARIES} sbnanaobj_standardrecord shared detsys trees)
target_include_directories(run_systematics PRIVATE include/ ${SBNANAOBJ_INC} ${ROOT_INCLUDE_DIRS})

# Add ROOT definitions
add_definitions(${ROOT_CXX_FLAGS})