cmake_minimum_required(VERSION 3.22)

# ------------------------------------------------------------------------------
#! Disallow building cmake project inside the source dir.
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
message(FATAL_ERROR
  "In-source builds are not permitted in this repository!\n"
  "Make a separate folder for building, and make sure it's outside the "
  "source dir. "
  "But first you'll need to remove the CMake cache files already created:\n"
  "cd ${CMAKE_SOURCE_DIR} && rm -r CMakeCache.txt CMakeFiles\n")
endif()

# These options to CMake are undocumented and cannot be relied upon.
# In versions that support them, they cause CMake to disallow in-source builds
# and modifications to the source directory respectively.
# They must be set *before* the 'project(...)' command is called.
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)


# ------------------------------------------------------------------------------
# Default to using ninja for building
set(CMAKE_GENERATOR "Ninja")

# ------------------------------------------------------------------------------
# if user specifies the build type, use theirs, otherwise RelWithDebInfo.
# Note: this needs to be defined before the project() call, otherwise
# it remains unset.
if(NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "")
endif()

# Constrain build type drop-down to contain only valid build types
set_property(CACHE CMAKE_BUILD_TYPE
  PROPERTY STRINGS
  Release Debug RelWithDebInfo MinSizeRel Coverage
)

# ------------------------------------------------------------------------------
project(nuSTORM LANGUAGES CXX)

# C++ standard for the entire project
set(CMAKE_CXX_STANDARD 17 CACHE STRING "")
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "")


add_subdirectory(analysis)