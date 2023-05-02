cmake_minimum_required(VERSION 3.0.0)
project(waterKernel VERSION 0.1.0)

set(CMAKE_BIULD_TYPE Debug)

include(CTest)
enable_testing()

# Add phase information
# Phase 0 is to obtain vanilla hits collection
# Phase 1 is to obtain vanilla dose distribution
# Phase 2 is to obtain the point dose kernel
set(PHASE "0")
configure_file("${PROJECT_SOURCE_DIR}/include/config.h.in" "${PROJECT_SOURCE_DIR}/include/config.h")

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

add_executable(waterKernel main.cpp ${sources} ${headers})

find_package(Boost REQUIRED COMPONENTS program_options)
target_include_directories(waterKernel PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(waterKernel PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(waterKernel ${Boost_LIBRARIES})

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
target_include_directories(waterKernel PRIVATE ${Geant4_INCLUDE_DIR})
target_link_libraries(waterKernel ${Geant4_LIBRARIES})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
