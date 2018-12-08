cmake_minimum_required(VERSION 3.0)
project(pcr-tdma)
list(APPEND CUDA_NVCC_FLAGS "-std=c++11")

CUDA_ADD_EXECUTABLE(pcr-tdma Test_Main.c)
set_property(TARGET pcr-tdma PROPERTY CUDA_STANDARD 11)
target_link_libraries(pcr-tdma
    PCR::PCR
    TDMA::TDMA
    utils::utils
)