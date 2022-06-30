file(GLOB SOURCES ./*.cpp ./*.cu)
file(GLOB HEADERS ./*.h)
add_library(geom ${SOURCES} ${HEADERS})

find_package(CUDA REQUIRED)

target_link_libraries(geom PRIVATE args optim ${CUDA_LIBRARIES})
target_include_directories(geom 
    PUBLIC .)