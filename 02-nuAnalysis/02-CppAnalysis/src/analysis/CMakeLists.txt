cmake_minimum_required(VERSION 3.22)

# --------------------------------------------------------------
project( analysis )

# Locate the ROOT package and define a number of useful targets and variables.
find_package(ROOT REQUIRED COMPONENTS RIO Net)



# Create library for the nuSTORM analysis code
add_library(nuStorm_analysis)


target_sources(nuStorm_analysis
  PRIVATE
  src/nu_analysis.cc
  src/run_control.cc
)

target_include_directories(nuStorm_analysis
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)


target_link_libraries(nuStorm_analysis
  PUBLIC
    ROOT::RIO 
    ROOT::Net)

# add_subdirectory(tests)