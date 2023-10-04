cmake_minimum_required(VERSION 3.0.0)
project(boxScore VERSION 0.1.0)

# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

add_executable(boxScore main.cpp ${sources} ${headers})

find_package(Boost REQUIRED COMPONENTS program_options filesystem)
target_include_directories(boxScore PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(boxScore PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(boxScore ${Boost_LIBRARIES})

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})

target_include_directories(boxScore PRIVATE ${Geant4_INCLUDE_DIR})
target_link_libraries(boxScore ${Geant4_LIBRARIES})