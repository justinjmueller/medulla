cmake_minimum_required(VERSION 3.12)
include(FetchContent)
project(run_systematics)

add_compile_options(-Wall -Werror -Wno-error=maybe-uninitialized -Ofast -g)

# Find packages
find_package(ROOT REQUIRED)
find_package(sbnanaobj)


set(SBNANAOBJ_INCLUDE_DIRS "$ENV{MRB_INSTALL}/sbnanaobj/$ENV{SBNANAOBJ_VERSION}/include")
set(SBNANAOBJ_LIBRARY_DIRS "$ENV{MRB_INSTALL}/sbnanaobj/$ENV{SBNANAOBJ_VERSION}/slf7.x86_64.e26.prof/lib")

file(GLOB SYSINC "include/*.h")


add_library(weight_reader SHARED src/weight_reader.cc)
target_link_libraries(weight_reader PRIVATE ${ROOT_LIBRARIES} shared)
target_include_directories(weight_reader PRIVATE include/ ${ROOT_INCLUDE_DIRS})

# Library for generic systematics
add_library(systematics SHARED src/systematic.cc)
target_link_libraries(systematics PRIVATE ${ROOT_LIBRARIES} shared)
target_include_directories(systematics PRIVATE include/ ${ROOT_INCLUDE_DIRS})

# Library for detector systematics
add_library(detsys SHARED src/detsys.cc)
target_link_libraries(detsys PRIVATE ${ROOT_LIBRARIES} shared)
target_include_directories(detsys PRIVATE include/ ${ROOT_INCLUDE_DIRS})

# Library for tree-handling code
add_library(trees SHARED src/trees.cc)
target_link_libraries(trees PRIVATE ${ROOT_LIBRARIES} shared weight_reader systematics detsys)
target_include_directories(trees PRIVATE include/ ${ROOT_INCLUDE_DIRS})

# Add the executable target
add_executable(run_systematics src/main.cc ${SYSINC})

# Link the ROOT libraries to the target
target_link_libraries(run_systematics ${ROOT_LIBRARIES} ${SBNANAOBJ_LIBRARY_DIRS}/libsbnanaobj_StandardRecord.so shared detsys trees)

# Include the ROOT headers
include_directories(${ROOT_INCLUDE_DIRS} ${SBNANAOBJ_INCLUDE_DIRS} include/)

# Add ROOT definitions
add_definitions(${ROOT_CXX_FLAGS})