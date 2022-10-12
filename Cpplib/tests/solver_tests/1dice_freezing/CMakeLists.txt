project(1dice_freezing)

# 1d ice solver 
add_executable(1dice_freezing 1dice_freezing.cpp)
target_link_libraries(1dice_freezing icethermo)