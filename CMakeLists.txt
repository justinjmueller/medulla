cmake_minimum_required(VERSION 3.14)
project(spine_anaplot_tools LANGUAGES C CXX)

# Central environment path resolution
set(SBNANA_BASE "$ENV{MRB_INSTALL}/sbnana/$ENV{SBNANA_VERSION}")
set(SBNANAOBJ_BASE "$ENV{MRB_INSTALL}/sbnanaobj/$ENV{SBNANAOBJ_VERSION}")
set(SRPROXY_INC "$ENV{SRPROXY_INC}")
set(OSCLIB_INC "$ENV{OSCLIB_INC}")
set(OSCLIB_LIB "$ENV{OSCLIB_LIB}")
set(EIGEN_INC "$ENV{EIGEN_INC}")

set(SBNANA_INC "${SBNANA_BASE}/include" CACHE INTERNAL "")
set(SBNANA_LIB "${SBNANA_BASE}/slf7.x86_64.e26.prof/lib" CACHE INTERNAL "")
set(SBNANAOBJ_INC "${SBNANAOBJ_BASE}/include" CACHE INTERNAL "")
set(SBNANAOBJ_LIB "${SBNANAOBJ_BASE}/slf7.x86_64.e26.prof/lib" CACHE INTERNAL "")

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