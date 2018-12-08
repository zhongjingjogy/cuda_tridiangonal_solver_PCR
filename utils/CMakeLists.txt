cmake_minimum_required(VERSION 3.5)
project(utils VERSION 1.0.0 LANGUAGES CXX)

##############################################
# Declare dependencies

##############################################
# Create target and set properties

add_library(utils
    INTERFACE
)

#Add an alias so that library can be used inside the build tree, e.g. when testing
add_library(utils::utils ALIAS utils)

#Set target properties
target_include_directories(utils
    INTERFACE
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)