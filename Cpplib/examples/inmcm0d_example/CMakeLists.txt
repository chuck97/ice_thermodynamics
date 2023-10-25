cmake_minimum_required (VERSION 3.17.0)

project(itinmcm0d_example)

# C/C++, Fortran compilers
enable_language (Fortran)

add_executable(itinmcm0d_example main.f90)

target_link_libraries(itinmcm0d_example icethermo)