cmake_minimum_required (VERSION 3.17.0)

project(itinmcm0d_example)

# C/C++, Fortran compilers
enable_language (Fortran)

add_executable(itinmcm0d_example main.f90)

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
	message(STATUS "Compiling itinmcm0d example with Intel")
	target_link_options(itinmcm0d_example PRIVATE -fPIC -nostartfiles)
	target_link_libraries(itinmcm0d_example PRIVATE icethermo -lstdc++)
else()
	message(STATUS"Compiling itinmcm example with GNU")
	target_link_libraries(itinmcm0d_example PRIVATE icethermo)
endif()
