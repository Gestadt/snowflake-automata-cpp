cmake_minimum_required(VERSION 3.28)
project(snowflake_automata_cpp)

set(CMAKE_CXX_STANDARD 17)


if (DEFINED EMSCRIPTEN)
	add_executable(snowflake_automata library.cpp library.h)
	set(EMSCRIPTEN_LINK_FLAGS "-s EXTRA_EXPORTED_RUNTIME_METHODS=['cwrap','getValue']")


	set_target_properties(snowflake_automata PROPERTIES COMPILE_FLAGS "-Os -s SIDE_MODULE=1 ")
	set_target_properties(snowflake_automata PROPERTIES LINK_FLAGS    "-Os -s WASM=1 --no-entry ${EMSCRIPTEN_LINK_FLAGS}")
else()
    add_library(snowflake_automata_cpp SHARED library.cpp)
endif()
