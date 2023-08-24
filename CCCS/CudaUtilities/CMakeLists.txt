file(GLOB SOURCES *.cpp *.cc *.cxx *.cu )
file(GLOB HEADERS *.h *.hpp *.cuh )

add_library(CudaUtilities STATIC
    ${SOURCES} ${HEADERS}
    )
target_compile_features(CudaUtilities PUBLIC cxx_std_11)
# set_target_properties(CudaUtilities PROPERTIES LINKER_LANGUAGE CUDA)
set_target_properties(CudaUtilities PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
