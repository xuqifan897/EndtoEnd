file(GLOB SOURCES ./*.cpp ./*.cu)
file(GLOB HEADERS ./*.h)
add_library(args ${SOURCES} ${HEADERS})

find_package(CUDA REQUIRED)

target_link_libraries(args PUBLIC ${CUDA_LIBRARIES})
target_include_directories(args 
    PUBLIC .)