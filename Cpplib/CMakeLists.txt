cmake_minimum_required (VERSION 3.17.0)

project(IceThermo)

# C/C++ compilers
enable_language (CXX C)

# set C++ standart
if(NOT "${CMAKE_CXX_STANDARD}")
  set(CMAKE_CXX_STANDARD 20)
endif()

include_directories("source/headers")

# options
option(USE_FP64 "Use 64-bit floating-point types (usually on)" ON)
option(COMPILE_TESTS "Compile tests" ON)
option(COMPILE_EXAMPLES "Compile examples" OFF)

if (COMPILE_TESTS)
  MESSAGE(STATUS "Test compiling ENABLED")
endif()

if (COMPILE_EXAMPLES)
  MESSAGE(STATUS "Examples compiling ENABLED")
endif()

# add subdirectories
add_subdirectory(source)

# make library
add_library(icethermo STATIC ${SOURCE} ${HEADER})

# compile tests
if (COMPILE_TESTS)
  add_subdirectory(tests)
endif()

# compile examples
if (COMPILE_EXAMPLES)
  add_subdirectory(examples)
endif()