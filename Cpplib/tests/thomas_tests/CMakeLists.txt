project(thomas_test)

set(SOURCE main.cpp)
add_executable(thomas_test ${SOURCE})
target_link_libraries(thomas_test icethermo)