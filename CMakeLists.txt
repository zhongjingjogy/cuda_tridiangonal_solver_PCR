cmake_minimum_required(VERSION 3.0)
project(cuda_test LANGUAGES CXX)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set (CMAKE_CXX_STANDARD 11)

find_package(CUDA REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})

list(APPEND CUDA_NVCC_FLAGS "-std=c++11")

CUDA_ADD_EXECUTABLE(pcr-test PCR_Class.cu PCR_Device_Functions.cu PCR_TEST.cu)

target_link_libraries(pcr-test
    ${CUDA_cusparse_LIBRARY} 
)

add_executable(tdma-test TDMA_TEST.cpp)
