cmake_minimum_required (VERSION 3.17.0)

project(IceThermo VERSION 1.0.0)

# C/C++ compilers
enable_language (CXX C)

# set C++ standart
if(NOT "${CMAKE_CXX_STANDARD}")
  set(CMAKE_CXX_STANDARD 17)
endif()

# options
option(USE_JSON_OUTPUT "Link library with nlohmann::json" OFF)
option(COMPILE_TESTS "Compile tests" OFF)
option(COMPILE_EXAMPLES "Compile examples" OFF)
option(COMPILE_ITSLAV "Compile IceThermo intefaces for SLAV model" OFF)
option(COMPILE_ITINMCM "Compile IceThermo intefaces for INMCM model" OFF)

# add Fortran support if SLAV interfaces are compiled
if (COMPILE_ITSLAV)
  MESSAGE(STATUS "SLAV interfaces compilation ENABLED")
  enable_language (Fortran)
endif()

# add Fortran support if INMCM interfaces are compiled
if (COMPILE_ITINMCM)
  MESSAGE(STATUS "INMCM interfaces compilation ENABLED")
  enable_language (Fortran)
endif()

# add subdirectory with source code
add_subdirectory(source)

if (COMPILE_ITSLAV)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()

# make library
add_library(icethermo STATIC ${SOURCE} ${HEADER})

# include directories
target_include_directories(icethermo PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(icethermo PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/source/headers")

# to install library
set_property(TARGET icethermo PROPERTY PUBLIC_HEADER
  "${PROJECT_SOURCE_DIR}/source/headers/constants.hpp"
  "${PROJECT_SOURCE_DIR}/source/headers/defines.hpp"
  "${PROJECT_SOURCE_DIR}/source/headers/icethermo.hpp"
  "${PROJECT_SOURCE_DIR}/source/headers/matvec.hpp"
  "${PROJECT_SOURCE_DIR}/source/headers/mesh.hpp"
  "${PROJECT_SOURCE_DIR}/source/headers/solver.hpp"
  "${PROJECT_SOURCE_DIR}/source/headers/tools.hpp"
  "${PROJECT_BINARY_DIR}/itslav.mod")

install(TARGETS icethermo 
	      ARCHIVE DESTINATION lib
	      LIBRARY DESTINATION lib
	      RUNTIME DESTINATION bin
	      PUBLIC_HEADER DESTINATION include)

if (COMPILE_ITSLAV)
  install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION fortran)
endif()


# link library with nlohmann::json
if (USE_JSON_OUTPUT)
  target_compile_definitions(icethermo PUBLIC USE_JSON_OUTPUT)
  find_package(nlohmann_json 3.2.0 REQUIRED)
  target_link_libraries(icethermo PUBLIC nlohmann_json::nlohmann_json)
endif()

# compile tests
if(COMPILE_TESTS)
  MESSAGE(STATUS "Test compiling ENABLED")
  add_subdirectory(tests)
endif(COMPILE_TESTS)

# compile examples
if (COMPILE_EXAMPLES)
  MESSAGE(STATUS "Examples compiling ENABLED")
  add_subdirectory(examples)
endif()
