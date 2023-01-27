# Project Build Settings common to all subprojects

# include (header) guard
if(RS4PI_DOSECALC_BUILDDEFS_INCLUDED)
    return()
endif(RS4PI_DOSECALC_BUILDDEFS_INCLUDED)
set(RS4PI_DOSECALC_BUILDDEFS_INCLUDED true)

# Find Dependent Packages
## REPLACED BY CMAKE LANGUAGE SUPPORT FOR CUDA (CMAKE >= v3.8)
#SET(CUDA_TOOLKIT_ROOT_DIR /usr/local/cuda)
# find_package(CUDA)
# IF(NOT CUDA_FOUND)
#     message(FATAL_ERROR
#             "CUDA Toolkit not found: explicitly set install dir using:\n"
#             "  cmake <src-dir> -DCUDA_TOOLKIT_ROOT_DIR=<cuda-install-dir>"
#        )
# ENDIF()

find_package( Boost REQUIRED COMPONENTS filesystem unit_test_framework)


#SET(DCMTK_DIR /opt/dcmtk-3.6.1)
find_package(DCMTK)
IF(NOT DCMTK_FOUND)
    message(FATAL_ERROR
            "DCMTK not found: explicitly set install dir using:\n"
            "  cmake <src-dir> -DDCMTK_DIR=<dcmtk-install-dir>"
	   )
ENDIF()
find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
IF(NOT HDF5_FOUND)
    message(FATAL_ERROR
            "HDF5 not found: explicitly set install dir using:\n"
            "  cmake <src-dir> -DHDF5_DIR=<hdf5-install-dir"
        )
ENDIF()

# generate shortname
get_filename_component(CUDA_ROOT_DIR "${CMAKE_CUDA_COMPILER}/../../" ABSOLUTE)
if(NOT IS_DIRECTORY ${CUDA_ROOT_DIR})
    message(FATAL_ERROR "CUDA_ROOT_DIR=\"${CUDA_ROOT_DIR}\": directory does not exist" )
endif()

# project-wide include dirs
include_directories(
    ./include
    ./include/cuda-common
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    # provides "cuda_runtime.h"
    ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    )
# project-wide lib directories
link_directories(
    ${HDF5_LIBRARY_DIRS}
    )

# CUDA NVCC Args
set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS}
    -cudart=static
    )
# use c++11 in cuda code (passed to nvcc)
set(CMAKE_CUDA_STANDARD 11)

set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS}
    # at least one PTX backend (code=compute_XX) should be included to maintain compatability with newer \
    # architectures (CC) using JIT compilation
    #-arch=sm_30;  # shorthand for following two lines
    -gencode=arch=compute_35,code=compute_35 # GeForce GTX >=600 series PTX

    # adds additional cubin (backend) compilation targets, bypassing JIT compilation if executing hw matches CC
    #-gencode=arch=compute_20,code=sm_20
    # -gencode=arch=compute_30,code=sm_30     # GeForce GTX 600 series cubin
    #-gencode=arch=compute_35,code=sm_35
    #-gencode=arch=compute_50,code=sm_50
    -gencode=arch=compute_52,code=sm_52     # GeForce GTX Titan X cubin
    #-gencode=arch=compute_53,code=sm_53    # Titan X cubin (w/ support for 16bit floating point "half" types)
    # -gencode=arch=compute_70,code=sm_70     # volta v100 cubin
    )

# set build type
IF(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
ENDIF()
message("\nBuild-type: \"${CMAKE_BUILD_TYPE}\"\n"
        "  (To specify a build type please use the cmake argument \"-DCMAKE_BUILD_TYPE=<type>\" and choose from the list:)\n"
        "  ['Release', 'Debug', 'Debug-GPU', 'Memcheck', 'Profile']\n")
IF(CMAKE_BUILD_TYPE MATCHES Debug-GPU)
    # generate debug information
    set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS}
        -g; -G;
        -Xptxas; -O0;
        -keep
        )
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}
        -g; -O0;
        )
    message(">> Warning: Device code optimizations disabled")
    message(">> Warning: Execution will be much slower in Debug mode")
ELSEIF(CMAKE_BUILD_TYPE MATCHES Debug)
    # generate debug information
    set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS}
        -g
        -keep
        -Xcompiler -O0; -Xlinker -O0
        )
    message(">> Warning: Device code optimizations disabled")
    message(">> Warning: Execution will be much slower in Debug mode")
ELSEIF(CMAKE_BUILD_TYPE MATCHES Memcheck)
    set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS}
        # valgrind host-side memcheck
        # cuda-memcheck GPU memcheck
        -g; -lineinfo;
        )
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}
        -g; -O0;
        )
ELSEIF(CMAKE_BUILD_TYPE MATCHES Profile)
    set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS}
        # valgrind host-side memcheck; callgrind/gprof CPU profiling
        # nvprof/nvvp GPU profiling
        -g; -pg; -lineinfo
        )
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}
        -g; -pg; -O3;
        )
    add_definitions(
        -DUSE_NVTX
        )
ELSEIF(CMAKE_BUILD_TYPE MATCHES Release)
    set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS}
        # -use_fast_math;
        )
ELSE()
    message(FATAL_ERROR "Invalid build-type: \"${CMAKE_BUILD_TYPE}\" specified")
ENDIF()

# convert cmake lists (;-separated) to strings (space-separated)
string(REPLACE ";" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
string(REPLACE ";" " " CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS}")
message(">> CMAKE_CXX_FLAGS =  ${CMAKE_CXX_FLAGS}")
message(">> CMAKE_CUDA_FLAGS = ${CMAKE_CUDA_FLAGS}")


# include git commit hash in version number
IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  FIND_PACKAGE(Git)
  IF(GIT_FOUND)
    EXECUTE_PROCESS(
      # COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
      COMMAND bash "-c" "${GIT_EXECUTABLE} describe --tags --always | sed 's/.*_v//; s/-[0-9]\\+-/-/'"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE "dosecalc_BUILD_VERSION"
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE( STATUS "Git version: v${dosecalc_BUILD_VERSION}" )
  ELSE(GIT_FOUND)
    SET(dosecalc_BUILD_VERSION 0)
  ENDIF(GIT_FOUND)
ENDIF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/include/version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/include/version.h @ONLY)
