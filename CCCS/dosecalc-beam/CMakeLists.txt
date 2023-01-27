add_executable(dosecalc-beam
    main.cu
    )
target_compile_features(dosecalc-beam PUBLIC cxx_std_11)
set_target_properties(dosecalc-beam PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(dosecalc-beam
    Utilities
    DoseCalcIO
    DoseCalcAlgo
    CudaUtilities
    ${Boost_LIBRARIES}
    -lpthread
    )

install(TARGETS dosecalc-beam
        RUNTIME DESTINATION bin
        )
