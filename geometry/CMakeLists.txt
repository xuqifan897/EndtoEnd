file(GLOB SOURCES ./*.cpp)
file(GLOB HEADERS ./*.h)
add_library(geom ${SOURCES} ${HEADERS})

target_link_libraries(geom PRIVATE args)
target_include_directories(geom 
    PUBLIC .)