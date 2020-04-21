cmake_minimum_required(VERSION 3.9.0)
project(geoEff)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${geoEff_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${geoEff_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${geoEff_SOURCE_DIR}/bin)

add_subdirectory(pybind11)
add_subdirectory(src)
add_subdirectory(app)

#
#pybind11_add_module(src/example.cpp)
