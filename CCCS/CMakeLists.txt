cmake_minimum_required(VERSION 3.14)
set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
set(CMAKE_CUDA_ARCHITECTURES 89)
project(DOSECALC_GPU LANGUAGES C CXX CUDA)

include(CMakeDefs.cmake)

# recurse into subdirs
add_subdirectory(Utilities)
add_subdirectory(RTClasses)
add_subdirectory(DoseCalcIO)
add_subdirectory(DoseCalcAlgo)
add_subdirectory(CudaUtilities)
add_subdirectory(dosecalc-preprocess)
add_subdirectory(dosecalc-beamlet)
add_subdirectory(dosecalc-beam)

enable_testing()
add_subdirectory(unit_tests)

install(DIRECTORY data/ctlut
                  data/dsa
                  data/spectra
        TYPE DATA
        )
install(DIRECTORY doc/
        TYPE DOC
        )
install(DIRECTORY extern
        DESTINATION .
        )
