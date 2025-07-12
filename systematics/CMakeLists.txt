cmake_minimum_required(VERSION 3.12)
include(FetchContent)
project(run_systematics)

add_compile_options(-Wall -Werror -Wno-error=maybe-uninitialized -Ofast -g)

# Find packages
find_package(ROOT REQUIRED)
find_package(sbnanaobj)

FetchContent_Declare(
    tomlplusplus
    GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
    GIT_TAG v3.4.0
)
FetchContent_MakeAvailable(tomlplusplus)

set(SBNANAOBJ_INCLUDE_DIRS "$ENV{MRB_INSTALL}/sbnanaobj/$ENV{SBNANAOBJ_VERSION}/include")
set(SBNANAOBJ_LIBRARY_DIRS "$ENV{MRB_INSTALL}/sbnanaobj/$ENV{SBNANAOBJ_VERSION}/slf7.x86_64.e26.prof/lib")

file(GLOB SYSINC "include/*.h")

# Library for configuration details
add_library(configuration SHARED src/configuration.cc)
target_link_libraries(configuration PRIVATE tomlplusplus::tomlplusplus)
target_include_directories(configuration PRIVATE include/ tomlplusplus/include)

# Library for WeightReader class
add_library(weight_reader SHARED src/weight_reader.cc)
target_link_libraries(weight_reader PRIVATE ${ROOT_LIBRARIES} configuration)
target_include_directories(weight_reader PRIVATE include/ ${ROOT_INCLUDE_DIRS})

# Library for generic systematics
add_library(systematics SHARED src/systematic.cc)
target_link_libraries(systematics PRIVATE ${ROOT_LIBRARIES} configuration tomlplusplus::tomlplusplus)
target_include_directories(systematics PRIVATE include/ ${ROOT_INCLUDE_DIRS} tomlplusplus/include)

# Library for detector systematics
add_library(detsys SHARED src/detsys.cc)
target_link_libraries(detsys PRIVATE ${ROOT_LIBRARIES} configuration tomlplusplus::tomlplusplus)
target_include_directories(detsys PRIVATE include/ ${ROOT_INCLUDE_DIRS} tomlplusplus/include)

# Library for tree-handling code
add_library(trees SHARED src/trees.cc)
#add_library(trees SHARED src/trees_gundam.cc)
target_link_libraries(trees PRIVATE ${ROOT_LIBRARIES} configuration weight_reader systematics detsys tomlplusplus::tomlplusplus)
target_include_directories(trees PRIVATE include/ ${ROOT_INCLUDE_DIRS} tomlplusplus/include)

# Add the executable target
add_executable(run_systematics src/main.cc ${SYSINC})
#add_executable(run_systematics src/main_gundam.cc ${SYSINC})

# Link the ROOT libraries to the target
target_link_libraries(run_systematics ${ROOT_LIBRARIES} ${SBNANAOBJ_LIBRARY_DIRS}/libsbnanaobj_StandardRecord.so tomlplusplus::tomlplusplus configuration detsys trees)

# Include the ROOT headers
include_directories(${ROOT_INCLUDE_DIRS} ${SBNANAOBJ_INCLUDE_DIRS} include/ tomlplusplus/include)

# Add ROOT definitions
add_definitions(${ROOT_CXX_FLAGS})