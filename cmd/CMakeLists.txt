# minus/cmd/CMakeLists.txt
find_package(Threads REQUIRED)
link_libraries(minus)
add_executable(minus-chicago minus-chicago.cxx)
target_link_libraries(minus-chicago Threads::Threads)

add_executable(minus-cleveland minus-cleveland.cxx)
target_link_libraries(minus-cleveland Threads::Threads)

add_executable(minus-2v5p minus-2v5p.cxx)
target_link_libraries(minus-2v5p Threads::Threads)
