cmake_minimum_required(VERSION 3.10)
project(Rasterizer)

set(CMAKE_CXX_STANDARD 17)

set(src_directory "source")
set(
    SOURCES
    ${src_directory}/main.cpp
    ${src_directory}/rasterizer.hpp
    ${src_directory}/rasterizer.cpp
    ${src_directory}/Triangle.hpp
    ${src_directory}/Triangle.cpp
)
include_directories(third_party)

add_executable(
    Rasterizer 
    ${SOURCES}
)
add_subdirectory(third_party/lodepng)
target_link_libraries(Rasterizer PRIVATE lodepng)
