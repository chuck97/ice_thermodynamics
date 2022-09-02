project(secant_test)

set(SOURCE main.cpp)
add_executable(secant_test ${SOURCE})
target_link_libraries(secant_test icethermo)