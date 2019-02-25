# This is comp/tests/CMakeLists.txt

add_executable(comps_test_all
  # Driver source and utilities
  test_driver.cxx
  test_comps.cxx
)


target_link_libraries(comps_test_all ${VXL_LIB_PREFIX}vnl ${VXL_LIB_PREFIX}testlib)
add_test( NAME test_comps COMMAND comps_test_all test_comps )

add_executable(comps_test_include test_include.cxx)
target_link_libraries(comps_test_include comps)
