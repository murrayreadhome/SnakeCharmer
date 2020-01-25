cmake_minimum_required(VERSION 3.0)
project(SnakeCharmer)
set(SOURCES 
    SnakeCharmer.cpp
)
# https://cmake.org/cmake/help/latest/prop_tgt/CXX_STANDARD.html
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_definitions(TESTING_AT_HOME)

add_executable(SnakeCharmer ${SOURCES})
