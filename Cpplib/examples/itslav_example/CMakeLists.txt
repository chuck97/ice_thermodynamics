cmake_minimum_required (VERSION 3.17.0)

project(itslav_example)

# C/C++, Fortran compilers
enable_language (Fortran)

add_executable(itslav_example main.f90)

target_link_libraries(itslav_example icethermo)