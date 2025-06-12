cmake_minimum_required(VERSION 3.14)
project(cafana LANGUAGES CXX)

# Central environment path resolution
set(SBNANA_BASE "$ENV{MRB_INSTALL}/sbnana/$ENV{SBNANA_VERSION}")
set(SBNANAOBJ_BASE "$ENV{MRB_INSTALL}/sbnanaobj/$ENV{SBNANAOBJ_VERSION}")

# Set C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set optimization flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
#set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Tree)
include(${ROOT_USE_FILE})

# TODO: Move environment path resolution to top-level CMake or shared config module
set(SRPROXY_INC "$ENV{SRPROXY_INC}")
set(OSCLIB_INC "$ENV{OSCLIB_INC}")
set(OSCLIB_LIB "$ENV{OSCLIB_LIB}")
set(EIGEN_INC "$ENV{EIGEN_INC}")
set(SBNANA_INC "$ENV{MRB_INSTALL}/sbnana/$ENV{SBNANA_VERSION}/include")
set(SBNANAOBJ_INC "$ENV{MRB_INSTALL}/sbnanaobj/$ENV{SBNANAOBJ_VERSION}/include")
set(SBNANA_LIB "$ENV{MRB_INSTALL}/sbnana/$ENV{SBNANA_VERSION}/slf7.x86_64.e26.prof/lib")
set(SBNANAOBJ_LIB "$ENV{MRB_INSTALL}/sbnanaobj/$ENV{SBNANAOBJ_VERSION}/slf7.x86_64.e26.prof/lib")

# Library for the framework
add_library(framework SHARED src/framework.cc)
target_link_directories(framework PRIVATE ${SBNANA_LIB} ${SBNANAOBJ_LIB})
target_link_libraries(framework PRIVATE shared CAFAnaCore sbnanaobj_StandardRecordProxy)
target_include_directories(framework PRIVATE include/ ${SBNANA_INC} ${SBNANAOBJ_INC} ${SRPROXY_INC})
target_compile_features(framework PRIVATE cxx_std_17)

# Library for the framework testing code
add_library(test SHARED src/test.cc)
target_link_directories(test PRIVATE ${SBNANA_LIB} ${SBNANAOBJ_LIB})
target_link_libraries(test PRIVATE shared CAFAnaCore sbnanaobj_StandardRecord)
target_include_directories(test PRIVATE include/ ${SBNANA_INC} ${SBNANAOBJ_INC} ${SRPROXY_INC})
target_compile_features(test PRIVATE cxx_std_17)

# Add a common library
add_library(common INTERFACE)

# Include directories
target_include_directories(common INTERFACE
    ${ROOT_INCLUDE_DIRS}
    ${SRPROXY_INC}
    ${OSCLIB_INC}
    ${EIGEN_INC}
    ${SBNANA_INC}
    ${SBNANAOBJ_INC}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link directories
target_link_directories(common INTERFACE
    ${OSCLIB_LIB}
    ${SBNANA_LIB}
    ${SBNANAOBJ_LIB}
)

# Link required libraries
target_link_libraries(common INTERFACE
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Tree
    Minuit2
    sbnanaobj_StandardRecordProxy
    OscLib
    CAFAnaCore
    CAFAnaUnfold
    CAFAnaVars
    CAFAnaCuts
    CAFAnaExtrap
    CAFAnaExperiment
    CAFAnaAnalysis
    SBNAnaCuts
    shared
    framework
)

target_compile_features(common INTERFACE cxx_std_17)

add_executable(main src/main.cc)
target_link_libraries(main PRIVATE shared framework common)
target_include_directories(main PRIVATE . include/ ${SBNANA_INC} ${SBNANAOBJ_INC} ${SRPROXY_INC})

add_executable(validate src/validate.cc)
target_link_libraries(validate PRIVATE test common sbnanaobj_StandardRecord)
target_include_directories(validate PRIVATE include/ ${SBNANA_INC} ${SBNANAOBJ_INC} ${SRPROXY_INC})

option(BUILD_DOCS "Build documentation" OFF)

if(BUILD_DOCS)
    find_package(Doxygen QUIET)
    find_package(Sphinx QUIET)
    if(DOXYGEN_FOUND)
        message(STATUS "Doxygen found: enabling C++ API documentation.")
    else()
        message(WARNING "Doxygen not found: skipping C++ API documentation.")
    endif()

    if(SPHINX_FOUND)
        message(STATUS "Sphinx found: enabling Python + unified documentation.")
    else()
        message(WARNING "Sphinx not found: skipping Sphinx-based documentation.")
    endif()

    if(DOXYGEN_FOUND AND SPHINX_FOUND)
        add_subdirectory(docs)
    endif()
endif()