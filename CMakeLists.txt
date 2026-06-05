cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(polyhedral-template-matching VERSION 1.0.0 LANGUAGES CXX)

find_package(coretoolkit REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(structure-identification REQUIRED)
find_package(common-neighbor-analysis REQUIRED)

file(GLOB PTM_VENDOR_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/dependencies/ptm/*.cpp)

volt_add_plugin(polyhedral-template-matching
    SOURCES
        src/analysis/ptm_service.cpp
        src/analysis/ptm.cpp
        src/analysis/ptm_crystal_info_provider.cpp
        src/analysis/ptm_structure_analysis.cpp
    MAIN_SOURCE src/plugin.cpp
    VENDOR_SOURCES ${PTM_VENDOR_SOURCES}
    INCLUDE_DIRS dependencies/ptm
    DEPENDENCIES
        common-neighbor-analysis::common-neighbor-analysis
        structure-identification::structure-identification
        yaml-cpp::yaml-cpp
    INSTALL_VENDOR_HEADERS dependencies/ptm
)
