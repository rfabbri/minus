# This is minus/cmd/CMakeLists.txt
find_package(Threads REQUIRED)

add_executable(expminus-chicago expminus-chicago.cxx expminus.cxx expminus.h)
target_link_libraries(expminus-chicago minus Threads::Threads)

#add_compile_options(-Wall -Weffc++ -pedantic -std=c++0x)
set_source_files_properties(expminus.cxx PROPERTIES COMPILE_FLAGS "-ffast-math -fcx-limited-range -march=native")
