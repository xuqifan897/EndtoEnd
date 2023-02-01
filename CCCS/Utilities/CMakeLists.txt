file(GLOB SOURCES *.cpp *.cc *.cxx )
file(GLOB HEADERS *.h *.hpp )

add_library(Utilities ${SOURCES} ${HEADERS})
target_compile_features(Utilities PUBLIC cxx_std_11)
