file(GLOB SOURCES *.cpp *.cc *.cxx )
file(GLOB HEADERS *.h *.hpp )

add_library(DoseCalcIO STATIC
    ${SOURCES} ${HEADERS}
    )
target_compile_features(DoseCalcIO PUBLIC cxx_std_11)

target_include_directories(DoseCalcIO PUBLIC
    )
target_link_libraries(DoseCalcIO
    Utilities
    RTClasses

    ${Boost_LIBRARIES}

    # Link hdf5 library for sparse data I/O
    ${HDF5_CXX_HL_LIBRARIES} # -lhdf5_hl_cpp
    ${HDF5_CXX_LIBRARIES} # -lhdf5_cpp
    ${HDF5_C_HL_LIBRARIES} # -lhdf5_hl or -lhdf5_serial_hl
    ${HDF5_C_LIBRARIES} # -lhdf5 or -lhdf5_serial
    -lz
    )
