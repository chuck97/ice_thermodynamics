cmake_minimum_required(VERSION 3.17.0)

project(ice_thermodynamics)

enable_language (CXX)
set(CMAKE_CXX_COMPILER g++)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory("ThermoView")
add_subdirectory("Headers")

include_directories("Headers")

set(INMOST_SRC "...")
set(INMOST_BUILD "...")


include_directories(${INMOST_SRC}/Source/Headers)
include_directories(${INMOST_BUILD})

set(
	SOURCE
	${SOURCE}
	${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
	${HEADER}
)

set(EXEC_NAME ${CMAKE_PROJECT_NAME})

add_executable(${EXEC_NAME} ${SOURCE})

target_link_libraries(${EXEC_NAME} ${INMOST_BUILD}/libinmost.a)
