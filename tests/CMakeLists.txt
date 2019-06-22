# This is minus/tests/CMakeLists.txt

#--------------------------------------------------------------

set(minus_test_all_sources
  test_driver.cxx
  test_minus.cxx
)

aux_source_directory(Templates minus_test_all_sources)
add_executable(minus_test_all ${minus_test_all_sources})

target_link_libraries(minus_test_all minus ${VXL_LIB_PREFIX}vnl ${VXL_LIB_PREFIX}vbl ${VXL_LIB_PREFIX}testlib)
# add_test( NAME test_minus COMMAND minus_test_all test_minus )

add_executable(minus_test_include test_include.cxx)
target_link_libraries(minus_test_include minus)
