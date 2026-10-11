cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(polyhedral-template-matching VERSION 1.0.0 LANGUAGES CXX)

if(NOT TARGET coretoolkit::coretoolkit)
    find_package(coretoolkit REQUIRED)
endif()
if(NOT TARGET yaml-cpp::yaml-cpp)
    find_package(yaml-cpp REQUIRED)
endif()
if(NOT TARGET structure-identification::structure-identification)
    find_package(structure-identification REQUIRED)
endif()
if(NOT TARGET common-neighbor-analysis::common-neighbor-analysis)
    find_package(common-neighbor-analysis REQUIRED)
endif()

file(GLOB PTM_VENDOR_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/ptm/*.cpp)

volt_add_plugin(polyhedral-template-matching
    MAIN_SOURCE src/plugin.cpp
    VENDOR_SOURCES ${PTM_VENDOR_SOURCES}
    INCLUDE_DIRS dependencies/ptm
    DEPENDENCIES
        common-neighbor-analysis::common-neighbor-analysis
        structure-identification::structure-identification
        yaml-cpp::yaml-cpp
    INSTALL_VENDOR_HEADERS dependencies/ptm
)
