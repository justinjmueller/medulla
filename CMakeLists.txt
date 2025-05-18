cmake_minimum_required(VERSION 3.14)
project(spine_anaplot_tools LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optional: CMake helper files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Shared library first
add_subdirectory(shared)

# Subprojects
add_subdirectory(cafana)
add_subdirectory(systematics)

# Docs (optional)
option(BUILD_DOCS "Build documentation" OFF)

if(BUILD_DOCS)
    find_package(Doxygen QUIET)
    find_package(Sphinx QUIET)
    if(DOXYGEN_FOUND AND SPHINX_FOUND)
        add_subdirectory(docs)
    else()
        message(WARNING "Documentation requested, but Doxygen or Sphinx not found. Skipping docs build.")
    endif()
endif()