file(GLOB SOURCES *.cpp *.cc *.cxx *.cu )
file(GLOB HEADERS *.h *.hpp *.cuh )

add_library(DoseCalcAlgo STATIC
    ${SOURCES} ${HEADERS}
    )
target_compile_features(DoseCalcAlgo PUBLIC cxx_std_11)
# set_target_properties(DoseCalcAlgo PROPERTIES LINKER_LANGUAGE CUDA)
set_target_properties(DoseCalcAlgo PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(DoseCalcAlgo
    CudaUtilities
    )
