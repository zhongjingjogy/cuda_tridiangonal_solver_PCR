cmake_minimum_required(VERSION 3.0)
project(cuda_test)

find_package(CUDA REQUIRED)
list(APPEND CUDA_NVCC_FLAGS "-std=c++11") 
include_directories(${CUDA_INCLUDE_DIRS})

CUDA_ADD_EXECUTABLE(test_cuda Test_Main.cu PCR_Class.cu PCR_Device_Functions.cu)
set_property(TARGET test_cuda PROPERTY CUDA_STANDARD 11)
target_link_libraries(test_cuda ${CUDA_cusparse_LIBRARY} )

# set_property(TARGET test_cuda PROPERTY CXX_STANDARD 11)