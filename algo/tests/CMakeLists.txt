# This is minus/algo/tests/CMakeLists.txt

#--------------------------------------------------------------

set(minus_algo_test_all_sources
  test_driver.cxx
  #  test_expminus.cxx
  test_linear.cxx
)

aux_source_directory(Templates minus_algo_test_all_sources)
add_executable(minus_algo_test_all ${minus_algo_test_all_sources})

target_link_libraries(minus_algo_test_all minus ${VXL_LIB_PREFIX}vnl ${VXL_LIB_PREFIX}vbl ${VXL_LIB_PREFIX}testlib)
# add_test( NAME test_expminus COMMAND minus_algo_test_all test_expminus )
add_test( NAME test_linear COMMAND minus_algo_test_all test_linear )

add_executable(minus_algo_test_include test_include.cxx)
target_link_libraries(minus_algo_test_include minus)
