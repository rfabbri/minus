# This is comp/tests/CMakeLists.txt

#--------------------------------------------------------------

set(minus_test_all_sources
  test_driver.cxx
  test_minus.cxx
  test_expminus.cxx
  expminus.cxx
  expminus.h
)

aux_source_directory(Templates minus_test_all_sources)
add_executable(minus_test_all ${minus_test_all_sources})

target_link_libraries(minus_test_all minus ${VXL_LIB_PREFIX}vnl ${VXL_LIB_PREFIX}vbl ${VXL_LIB_PREFIX}testlib)
add_test( NAME test_minus COMMAND minus_test_all test_minus )
add_test( NAME exptest_minus COMMAND minus_test_all exptest_minus )

add_executable(minus_test_include test_include.cxx)
target_link_libraries(minus_test_include minus)

#--------------------------------------------------------------
  

add_executable(expminus-chicago expminus-chicago.cxx expminus.cxx expminus.h)
target_link_libraries(expminus-chicago minus)
  
#add_compile_options(-Wall -Weffc++ -pedantic -std=c++0x)
set_source_files_properties(expminus.cxx PROPERTIES COMPILE_FLAGS "-ffast-math -fcx-limited-range -march=native")
