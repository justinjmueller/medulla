cmake_minimum_required(VERSION 3.14)
project(spine_anaplot_tools LANGUAGES C CXX)

# Get the absolute paths
set(MEDULLA_SOURCE_DIR ${CMAKE_SOURCE_DIR})
set(MEDULLA_BUILD_DIR ${CMAKE_BINARY_DIR})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(CheckEnvironment) # Validate first
include(SetupEnvironment) # Set paths to external dependencies

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optional: CMake helper files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Shared library first
add_subdirectory(shared)

# Subprojects
add_subdirectory(selection)
add_subdirectory(systematics)

# Python unit tests
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_FOUND)
    add_custom_target(pytest
        COMMAND ${CMAKE_COMMAND} -E env
            MEDULLA_BUILD_DIR=${CMAKE_BINARY_DIR}/selection
            ${Python3_EXECUTABLE} -m pytest
        DEPENDS medulla validate
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Running Python unit tests"
        VERBATIM
    )
    message(STATUS "Python3 found (${Python3_EXECUTABLE}): 'make pytest' target enabled.")
else()
    message(WARNING "Python3 interpreter not found: 'make pytest' target disabled.")
endif()

# Docs (optional)
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

# Generate setup script
configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/setup_medulla.sh.in
    ${CMAKE_BINARY_DIR}/setup_medulla.sh
    @ONLY
)

# Make it executable
execute_process(
    COMMAND chmod 755 ${CMAKE_BINARY_DIR}/setup_medulla.sh
)

message(STATUS "Setup script generated: ${CMAKE_BINARY_DIR}/setup_medulla.sh")
message(STATUS "Source it with: source ${CMAKE_BINARY_DIR}/setup_medulla.sh")
