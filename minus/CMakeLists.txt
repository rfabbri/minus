# optional, header-only library can be used without this library
# library is useful when the function will be used in many places
# to avoid having to include Minus many times and slowdown the compiler
set (minus_sources
  minus.hxx minus.h chicago14a.hxx
)
aux_source_directory(Templates minus_sources)

# this library can be used if you want to precompile instantiations
# but header-only will also work, see README.md

include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_library( minus ${minus_sources} )
set_target_properties(minus PROPERTIES LINKER_LANGUAGE CXX)
