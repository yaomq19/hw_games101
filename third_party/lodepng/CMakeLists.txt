cmake_minimum_required(VERSION 3.10)

project(lodepng LANGUAGES CXX)

# Set the source files
set(SOURCES lodepng.cpp lodepng.h)

# Add the static library target
add_library(lodepng STATIC ${SOURCES})

# Include the header files
target_include_directories(lodepng PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Set the library to be built with debug information
set_target_properties(lodepng PROPERTIES DEBUG_POSTFIX "_d")