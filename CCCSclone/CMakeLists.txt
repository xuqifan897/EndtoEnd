cmake_minimum_required(VERSION 3.0.0)
project(CCCSclone VERSION 0.1.0 LANGUAGES CXX C CUDA)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_SOURCE_DIR ${CMAKE_SOURCE_DIR})
configure_file(${CMAKE_SOURCE_DIR}/include/old/configure.h.in ${CMAKE_SOURCE_DIR}/include/old/configure.h)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/**/*.cpp ${PROJECT_SOURCE_DIR}/src/**/*.cu)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/**/*.h ${PROJECT_SOURCE_DIR}/include/**/*.cuh)

add_executable(CCCSclone main.cu ${sources} ${headers})
target_include_directories(CCCSclone PRIVATE 
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/include/old
    ${PROJECT_SOURCE_DIR}/include/new)

find_package(CUDA REQUIRED)
set(CUDA_ARCHITECTURES 75)
target_include_directories(CCCSclone PRIVATE ${CUDA_INCLUDE_DIRS})

find_package(Boost REQUIRED COMPONENTS program_options filesystem)
target_include_directories(CCCSclone PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(CCCSclone PRIVATE Boost::program_options Boost::filesystem)

find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
target_include_directories(CCCSclone PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(CCCSclone PRIVATE ${HDF5_CXX_LIBRARIES})