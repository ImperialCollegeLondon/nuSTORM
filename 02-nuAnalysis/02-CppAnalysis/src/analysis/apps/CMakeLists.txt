# ----------------------------------------------------------------------

add_executable(nuAnalysis
               nuAnalysis-skeleton.cc
)

target_link_libraries(nuAnalysis
  PUBLIC
    ROOT::RIO 
    ROOT::Net
  PRIVATE
    nuStorm_analysis
)

install(
  TARGETS
    nuAnalysis
  EXPORT analysis
)
