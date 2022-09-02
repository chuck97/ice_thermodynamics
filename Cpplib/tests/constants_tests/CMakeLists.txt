project(constants_tests)

set(SOURCE main.cpp)
add_executable(constants_tests ${SOURCE})
target_link_libraries(constants_tests icethermo)