cmake_minimum_required(VERSION 3.5)
project(pcr VERSION 1.0.0 LANGUAGES CXX)

#Make sure that custom modules like FindRapidJSON are found
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)

##############################################
# Declare dependencies

##############################################
# Create target and set properties

add_library(pcr
    src/PCR_Class.cu
    src/PCR_Device_Functions.cu
)

#Add an alias so that library can be used inside the build tree, e.g. when testing
add_library(PCR::PCR ALIAS pcr)

#Set target properties
target_include_directories(pcr
    INTERFACE
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_compile_features(pcr INTERFACE cxx_auto_type)

##############################################
# Installation instructions

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/PCR)

install(TARGETS pcr
    EXPORT pcr-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

#This is required so that the exported target has the name PCR and not pcr
set_target_properties(pcr PROPERTIES EXPORT_NAME PCR)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

#Export the targets to a script
install(EXPORT pcr-targets
  FILE
    PCRTargets.cmake
  NAMESPACE
    PCR::
  DESTINATION
    ${INSTALL_CONFIGDIR}
)

#Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/PCRConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/PCRConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/PCRConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

#Install the config, configversion and custom find modules
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/PCRConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/PCRConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

##############################################
## Exporting from the build tree

export(EXPORT pcr-targets FILE ${CMAKE_CURRENT_BINARY_DIR}/PCRTargets.cmake NAMESPACE PCR::)

#Register package in user's package registry
export(PACKAGE PCR)

