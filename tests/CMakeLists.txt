# This is comp/tests/CMakeLists.txt


set(comps_test_all_sources
  test_driver.cxx
  test_comps.cxx
)
aux_source_directory(Templates comps_test_all_sources)

add_executable(comps_test_all ${comps_test_all_sources})


target_link_libraries(comps_test_all comps ${VXL_LIB_PREFIX}vnl ${VXL_LIB_PREFIX}testlib)
add_test( NAME test_comps COMMAND comps_test_all test_comps )

add_executable(comps_test_include test_include.cxx)
target_link_libraries(comps_test_include comps)

