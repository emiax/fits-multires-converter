cmake_minimum_required(VERSION 2.8)
include(ExternalProject)

project(FitsConverter)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${FitsConverter_SOURCE_DIR})

add_subdirectory(${FitsConverter_SOURCE_DIR}/ext/ghoul)
include(${FitsConverter_SOURCE_DIR}/ext/ghoul/ext/CopySharedLibraries.cmake)

add_subdirectory(${FitsConverter_SOURCE_DIR}/ext/basicj2k)

add_subdirectory(src)

set(CMAKE_BUILD_TYPE Debug)