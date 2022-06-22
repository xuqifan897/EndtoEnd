cmake_minimum_required(VERSION 3.23.1)
project(EndtoEnd3 VERSION 0.1.0 LANGUAGES CXX CUDA)

set(CMAKE_BUILD_TYPE Debug)
include(CTest)
enable_testing()

set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS}
    -gencode=arch=compute_61,code=sm_61)     # GeForce GTX 1080 TI

set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS}
    -g; -G;
    )
string(REPLACE ";" " " CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS}")

add_subdirectory(args)
add_subdirectory(geometry)

FIND_PACKAGE(Boost 1.65 COMPONENTS program_options REQUIRED)

add_executable(EndtoEnd3 main.cpp)
target_link_libraries(EndtoEnd3 PRIVATE args geom ${Boost_LIBRARIES})
target_include_directories(EndtoEnd3 PRIVATE ${Boost_INCLUDE_DIR})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
