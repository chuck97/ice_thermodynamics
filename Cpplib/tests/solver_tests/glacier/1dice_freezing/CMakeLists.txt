project(1dglacier_freezing)

# 1d ice solver 
add_executable(1dglacier_freezing 1dglacier_freezing.cpp)
target_link_libraries(1dglacier_freezing icethermo)