cmake_minimum_required(VERSION 3.5)
project(tdma VERSION 1.0.0 LANGUAGES CXX)

##############################################
# Declare dependencies

##############################################
# Create target and set properties

add_library(tdma
    INTERFACE
)

#Add an alias so that library can be used inside the build tree, e.g. when testing
add_library(TDMA::TDMA ALIAS tdma)

#Set target properties
target_include_directories(tdma
    INTERFACE
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)