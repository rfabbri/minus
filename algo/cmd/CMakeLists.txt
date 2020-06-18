# minus/cmd/CMakeLists.txt
find_package(Threads REQUIRED)
link_libraries(minus)

add_executable(expminus-chicago expminus-chicago.cxx)
target_link_libraries(expminus-chicago Threads::Threads)
