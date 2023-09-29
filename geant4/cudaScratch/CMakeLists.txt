# Specify the minimum CMake version required
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

# Create a project and specify its name
project(MyCudaProject)

# Find the CUDA package
find_package(CUDA REQUIRED)

# Set the CUDA compute capability (optional but recommended)
# For example, use compute capability 7.0 (Volta architecture)
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_89")

# Add the CUDA source file to the project
# Replace 'main.cu' with the actual name of your CUDA source file(s)
cuda_add_executable(my_cuda_program main.cpp)