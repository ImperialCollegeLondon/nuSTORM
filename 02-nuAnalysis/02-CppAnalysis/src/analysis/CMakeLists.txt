cmake_minimum_required(VERSION 3.22)

# --------------------------------------------------------------
project( nuStorm-analysis )

# Locate the ROOT package and define a number of useful targets and variables.
find_package(ROOT REQUIRED COMPONENTS
              Core
              RIO
              Net
              Hist
              Graf
              Graf3d
              Gpad
              Tree
              Rint
              Postscript
              Matrix
              Physics
              MathCore
              Thread
              MultiProc
)


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

# Add all of roots targets - this definitely isn't need and should be cut down
target_link_libraries(nuStorm_analysis
  PUBLIC
    ROOT::Core
    ROOT::RIO
    ROOT::Net
    ROOT::Hist
    ROOT::Graf
    ROOT::Graf3d
    ROOT::Gpad
    ROOT::Tree
    ROOT::Rint
    ROOT::Postscript
    ROOT::Matrix
    ROOT::Physics
    ROOT::MathCore
    ROOT::Thread
    ROOT::MultiProc
)

install(
  TARGETS
    nuStorm_analysis
  EXPORT nuStorm-analysis
)

add_subdirectory(apps)
# add_subdirectory(tests)