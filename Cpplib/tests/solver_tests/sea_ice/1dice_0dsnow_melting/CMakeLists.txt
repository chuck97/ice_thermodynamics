project(1dice_0dsnow_melting)

# 1d ice solver 
add_executable(1dice_0dsnow_melting 1dice_0dsnow_melting.cpp)
target_link_libraries(1dice_0dsnow_melting icethermo)