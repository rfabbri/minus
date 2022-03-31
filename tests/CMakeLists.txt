# This is minus/tests/CMakeLists.txt

#--------------------------------------------------------------

set(minus_test_all_sources
  test-driver.cxx
  test-minus.cxx
  test-minus-cleveland.cxx
  test-internals.cxx
  test-internals-cleveland.cxx
)

aux_source_directory(Templates minus_test_all_sources)
add_executable(minus_test_all ${minus_test_all_sources})

target_link_libraries(minus_test_all minus ${VXL_LIB_PREFIX}vnl ${VXL_LIB_PREFIX}vbl ${VXL_LIB_PREFIX}testlib)
add_test( NAME test-minus COMMAND minus_test_all test_minus )
add_test( NAME test_minus_cleveland COMMAND minus_test_all test_minus_cleveland )
add_test( NAME test-internals COMMAND minus_test_all test_internals )
add_test( NAME test-internals_cleveland COMMAND minus_test_all test_internals_cleveland )

add_executable(minus_test_include test-include.cxx)
target_link_libraries(minus_test_include minus)
