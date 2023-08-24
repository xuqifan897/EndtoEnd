find_package( OpenCV REQUIRED )
add_compile_definitions(OPENCV_VERSION=${OpenCV_VERSION} )
if(${OpenCV_VERSION} VERSION_LESS 3.0)
  add_compile_definitions(WITH_OPENCV2)
endif()

add_executable(dosecalc-preprocess
    main.cu
    fmapProcessing.cpp fmapProcessing.hpp
    cudaMemFunctions.cu cudaMemFunctions.cuh
    kernels.cuh
    debug.h
    argparser.h
    )

target_compile_features(dosecalc-preprocess PUBLIC cxx_std_11)
set_target_properties(dosecalc-preprocess PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# Project libs
target_link_libraries(dosecalc-preprocess PRIVATE
    ${OpenCV_LIBS}
    Utilities
    RTClasses
    DoseCalcIO
    DoseCalcAlgo
    CudaUtilities
    -lpthread
    )

install(TARGETS dosecalc-preprocess
        RUNTIME DESTINATION bin
        )
