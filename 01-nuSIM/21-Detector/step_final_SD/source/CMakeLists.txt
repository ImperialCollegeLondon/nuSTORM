# Setup the project

cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
project(triangles)

# Find Geant4 package and activate everything we want

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
	find_package(Geant4 REQUIRED ui_all vis_all)
else()
	find_package(Geant4 REQUIRED)
endif()

#Setup Geant4 include directories and compile definitions
# Setup include directory for this project
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources 
	${PROJECT_SOURCE_DIR}/src/stepHit.cc
	${PROJECT_SOURCE_DIR}/src/stepSD.cc
	${PROJECT_SOURCE_DIR}/src/trianglesActionInitialization.cc 
	${PROJECT_SOURCE_DIR}/src/trianglesDetectorConstruction.cc 
	${PROJECT_SOURCE_DIR}/src/trianglesEventAction.cc 
	${PROJECT_SOURCE_DIR}/src/trianglesMagneticField.cc 
	${PROJECT_SOURCE_DIR}/src/trianglesScintHit.cc 
	${PROJECT_SOURCE_DIR}/src/trianglesScintSD.cc 
	${PROJECT_SOURCE_DIR}/src/trianglesPrimaryGeneratorAction.cc 
	${PROJECT_SOURCE_DIR}/src/trianglesRunAction.cc)

file(GLOB headers 
	${PROJECT_SOURCE_DIR}/include/stepHit.hh
	${PROJECT_SOURCE_DIR}/include/stepSD.hh
	${PROJECT_SOURCE_DIR}/include/trianglesActionInitialization.hh 
	${PROJECT_SOURCE_DIR}/include/trianglesConstants.hh
	${PROJECT_SOURCE_DIR}/include/trianglesDetectorConstruction.hh 
	${PROJECT_SOURCE_DIR}/include/trianglesEventAction.hh 
	${PROJECT_SOURCE_DIR}/include/trianglesMagneticField.hh 
	${PROJECT_SOURCE_DIR}/include/trianglesScintHit.hh 
	${PROJECT_SOURCE_DIR}/include/trianglesScintSD.hh 
	${PROJECT_SOURCE_DIR}/include/trianglesPrimaryGeneratorAction.hh 
	${PROJECT_SOURCE_DIR}/include/trianglesRunAction.hh)

add_executable(triangles triangles.cc ${sources} ${headers})
target_link_libraries(triangles ${Geant4_LIBRARIES})
#Copy all the necessary scripts to the build directory

set(TRIANGLES_SCRIPTS init_vis.mac vis.mac run_100eve.mac)

foreach(_script ${TRIANGLES_SCRIPTS})
	configure_file(${PROJECT_SOURCE_DIR}/${_script} ${PROJECT_BINARY_DIR}/${_script} COPYONLY)
endforeach()

add_custom_target(Triangles DEPENDS triangles)

install(TARGETS triangles DESTINATION bin)

