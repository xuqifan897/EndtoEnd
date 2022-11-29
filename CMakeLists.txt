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
add_subdirectory(optimize)

FIND_PACKAGE(Boost 1.65 COMPONENTS system filesystem program_options REQUIRED)

# add_executable(EndtoEnd3 main.cpp)
# target_link_libraries(EndtoEnd3 PRIVATE args geom optim ${Boost_LIBRARIES})
# target_include_directories(EndtoEnd3 PRIVATE ${Boost_INCLUDE_DIR})

# add_executable(dosecalc dosecalc.cpp)
# target_link_libraries(dosecalc PRIVATE args geom optim ${Boost_LIBRARIES})
# target_include_directories(dosecalc PRIVATE ${Boost_INCLUDE_DIR})

add_executable(dose_calculation dose_calculation.cpp)
target_link_libraries(dose_calculation PRIVATE args geom optim ${Boost_LIBRARIES})
target_include_directories(dose_calculation PRIVATE ${Boost_INCLUDE_DIR})

add_executable(optimize_stationary optimize_stationary.cpp)
target_link_libraries(optimize_stationary PRIVATE args geom optim ${Boost_LIBRARIES})
target_include_directories(dose_calculation PRIVATE ${Boost_INCLUDE_DIR})

add_executable(optimize_stationary_smoothness optimize_stationary_smoothness.cpp)
target_link_libraries(optimize_stationary_smoothness PRIVATE args geom optim ${Boost_LIBRARIES})
target_include_directories(optimize_stationary_smoothness PRIVATE ${Boost_INCLUDE_DIR})

add_executable(optimize_dynamic optimize_dynamic.cpp)
target_link_libraries(optimize_dynamic PRIVATE args geom optim ${Boost_LIBRARIES})
target_include_directories(optimize_dynamic PRIVATE ${Boost_INCLUDE_DIR})

add_executable(optimize_dynamic_random optimize_dynamic_random.cpp)
target_link_libraries(optimize_dynamic_random PRIVATE args geom optim ${Boost_LIBRARIES})
target_include_directories(optimize_dynamic_random PRIVATE ${Boost_INCLUDE_DIR})

add_executable(optimize_annealing optimize_annealing.cpp)
target_link_libraries(optimize_annealing PRIVATE args geom optim ${Boost_LIBRARIES})
target_include_directories(optimize_annealing PRIVATE ${Boost_INCLUDE_DIR})

add_executable(optimize_annealing_schedule optimize_annealing_schedule.cpp)
target_link_libraries(optimize_annealing_schedule PRIVATE args geom optim ${Boost_LIBRARIES})
target_include_directories(optimize_annealing_schedule PRIVATE ${Boost_INCLUDE_DIR})

add_executable(optimize_annealing_constrain optimize_annealing_constrain.cpp)
target_link_libraries(optimize_annealing_constrain PRIVATE args geom optim ${Boost_LIBRARIES})
target_include_directories(optimize_annealing_constrain PRIVATE ${Boost_INCLUDE_DIR})

add_executable(benchmarking benchmarking.cpp)
target_link_libraries(benchmarking PRIVATE args geom optim ${Boost_LIBRARIES})
target_include_directories(benchmarking PRIVATE ${Boost_INCLUDE_DIR})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
