cmake_minimum_required(VERSION 3.0)
project(tdma-only)

add_executable(tdma-only Test_Main.cpp)
target_link_libraries(tdma-only
    TDMA::TDMA
    utils::utils
)

target_compile_features(tdma-only PUBLIC cxx_std_11)