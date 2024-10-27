cmake_minimum_required(VERSION 3.10)
project(hw1)

set(CMAKE_CXX_STANDARD 17)

set(src_directory ".")
set(
    SOURCES
    ${src_directory}/main.cpp
    ${src_directory}/rasterizer.hpp
    ${src_directory}/rasterizer.cpp
    ${src_directory}/Triangle.hpp
    ${src_directory}/Triangle.cpp
)
include_directories(../../third_party)

add_executable(
    ${PROJECT_NAME}
    ${SOURCES}
)

target_link_libraries(${PROJECT_NAME} PRIVATE lodepng)
