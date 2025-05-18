cmake_minimum_required(VERSION 3.14)

project(shared LANGUAGES CXX)

# Fetch tomlplusplus (shared for all projects)
include(FetchContent)
FetchContent_Declare(
    tomlplusplus
    GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
    GIT_TAG v3.4.0
)
FetchContent_MakeAvailable(tomlplusplus)

# Shared configuration library
add_library(shared SHARED
    src/configuration.cc
)

target_include_directories(shared
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(shared
    PUBLIC tomlplusplus::tomlplusplus
)

target_compile_features(shared PUBLIC cxx_std_17)