cmake_minimum_required (VERSION 3.17.0)

project(itinmcm1d_example)

# C/C++, Fortran compilers
enable_language (Fortran)

add_executable(itinmcm1d_example main.f90)

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
	message(STATUS "Compiling itinmcm0d example with Intel")
	set(CMAKE_CXX_LINK_EXECUTABLE "<CMAKE_Fortran_COMPILER> <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>")
	target_link_libraries(itinmcm1d_example PRIVATE icethermo -lstdc++)
else()
	message(STATUS"Compiling itinmcm example with GNU")
	target_link_libraries(itinmcm1d_example PRIVATE icethermo)
endif()
