cmake_minimum_required(VERSION 3.0.0)
project(waterPoint VERSION 0.1.0)

# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

# which phantom setting to use
# set(PHASE 0)  # to examine the tracking functionality
set(PHASE 1)  # to use the tracking functionality to score dose
configure_file("${PROJECT_SOURCE_DIR}/include/config.h.in" "${PROJECT_SOURCE_DIR}/include/config.h")

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

add_executable(waterPoint main.cpp ${sources} ${headers})

find_package(Boost REQUIRED COMPONENTS program_options filesystem)
target_include_directories(waterPoint PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(waterPoint PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(waterPoint ${Boost_LIBRARIES})

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" OFF)
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
target_include_directories(waterPoint PRIVATE ${Geant4_INCLUDE_DIR})
target_link_libraries(waterPoint ${Geant4_LIBRARIES})


# For point dose convolution
add_executable(slabPointTab slabPointTab.cpp ${sources} ${headers})
target_include_directories(slabPointTab PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(slabPointTab PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(slabPointTab ${Boost_LIBRARIES})
target_include_directories(slabPointTab PRIVATE ${Geant4_INCLUDE_DIR})
target_link_libraries(slabPointTab ${Geant4_LIBRARIES})