cmake_minimum_required(VERSION 3.28)
project(snowflake_automata_cpp)

set(CMAKE_CXX_STANDARD 17)

add_library(snowflake_automata_cpp SHARED library.cpp)
