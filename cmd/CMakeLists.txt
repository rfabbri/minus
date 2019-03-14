# minus/cmd/CMakeLists.txt
find_package(Threads REQUIRED)
#link_libraries(minus)
add_executable(minus-chicago minus-chicago.cxx)
target_link_libraries(minus-chicago Threads::Threads)
