project(solver_test)

# 1d ice solver 
add_executable(ice1d ice1d.cpp)
target_link_libraries(ice1d icethermo)