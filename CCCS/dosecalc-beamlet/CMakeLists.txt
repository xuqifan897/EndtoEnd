add_executable(dosecalc-beamlet
    main.cu
    sparsify_manager.cpp sparsify_manager.h
    )
target_compile_features(dosecalc-beamlet PUBLIC cxx_std_11)
set_target_properties(dosecalc-beamlet PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(dosecalc-beamlet
    Utilities
    DoseCalcIO
    DoseCalcAlgo
    CudaUtilities
    -lpthread
    -lnvToolsExt
    )

install(TARGETS dosecalc-beamlet
        RUNTIME DESTINATION bin
        )
