file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach(testSrc ${TEST_SRCS})
    get_filename_component(testName ${testSrc} NAME_WE)
    add_executable(${testName} ${testSrc} "../dosecalc-preprocess/fmapProcessing.cpp")
    target_compile_features(${testName} PUBLIC cxx_std_11)
    target_link_libraries(${testName} PRIVATE
        ${Boost_LIBRARIES}
        Utilities
        DoseCalcIO
        DoseCalcAlgo
        CudaUtilities
        )

    # set_target_properties(${testName} PROPERTIES
    #     RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin)

    add_test(NAME ${testName}
        COMMAND ./${testName})
    # add_test(NAME ${testName}
    #     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin
    #     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testName})
endforeach(testSrc)
