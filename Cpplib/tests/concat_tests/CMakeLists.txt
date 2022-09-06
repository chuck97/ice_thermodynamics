project(concat_test)

set(SOURCE main.cpp)
add_executable(concat_test ${SOURCE})
target_link_libraries(concat_test icethermo)