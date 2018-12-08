cmake_minimum_required(VERSION 3.0)
project(cuda_test)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set (CMAKE_CXX_STANDARD 11)

# add_subdirectory(./pcr)
add_subdirectory(./tdma)
add_subdirectory(./utils)
add_subdirectory(./test/tdma-only)