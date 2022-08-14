file(GLOB SOURCES ./*.cpp ./*.cu)
file(GLOB HEADERS ./*.h)
add_library(optim ${SOURCES} ${HEADERS})

find_package(CUDA REQUIRED)

target_link_libraries(optim PRIVATE args geom ${CUDA_LIBRARIES})
target_include_directories(optim PUBLIC .)