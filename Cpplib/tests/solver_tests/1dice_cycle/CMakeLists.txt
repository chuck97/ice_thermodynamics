project(1dice_cycle)

# 1d ice solver 
add_executable(1dice_cycle 1dice_cycle.cpp)
target_link_libraries(1dice_cycle icethermo)