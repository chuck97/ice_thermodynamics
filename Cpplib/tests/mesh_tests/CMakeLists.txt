project(mesh_test)

set(SOURCE main.cpp)
add_executable(mesh_test ${SOURCE})
target_link_libraries(mesh_test icethermo)