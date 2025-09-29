cmake_minimum_required(VERSION 3.12)
include(FetchContent)
project(systematics LANGUAGES CXX)

add_compile_options(-Wall -Werror -Wno-error=maybe-uninitialized -Ofast -g)

# Find packages
find_package(ROOT REQUIRED)
find_package(sbnanaobj)


add_library(weight_reader SHARED src/weight_reader.cc)
target_link_libraries(weight_reader PRIVATE ${ROOT_LIBRARIES} shared)
target_include_directories(weight_reader PRIVATE include/ ${ROOT_INCLUDE_DIRS} ${SBNANAOBJ_INC})

# Library for generic systematics
add_library(systematics SHARED src/systematic.cc)
target_link_libraries(systematics PRIVATE ${ROOT_LIBRARIES} shared)
target_include_directories(systematics PRIVATE include/ ${ROOT_INCLUDE_DIRS} ${SBNANAOBJ_INC})

# Library for detector systematics
add_library(detsys SHARED src/detsys.cc)
target_link_libraries(detsys PRIVATE ${ROOT_LIBRARIES} shared)
target_include_directories(detsys PRIVATE include/ ${ROOT_INCLUDE_DIRS} ${SBNANAOBJ_INC})

# Library for tree-handling code
add_library(trees SHARED src/trees.cc)
target_link_libraries(trees PRIVATE ${ROOT_LIBRARIES} shared weight_reader systematics detsys)
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