file(GLOB SOURCES ./*.cpp)
file(GLOB HEADERS ./*.h)
add_library(args ${SOURCES} ${HEADERS})

target_include_directories(args 
    PUBLIC .)