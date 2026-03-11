# optional, header-only library can be used without this library
# library is useful when the function will be used in many places
# to avoid having to include Minus many times and slowdown the compiler
set(minus_sources
  minus.hxx minus.h 
#
  chicago14a.h chicago-default.h chicago14a-default-data.h chicago14a-io.h
  chicago14a.hxx chicago14a-HxH.hxx chicago14a-Hxt.hxx chicago14a-default-data.hxx chicago14a-lsolve.hxx

  
  linecircle2a.h linecircle-default.h linecircle2a-default-data.h linecircle2a-io.h
  linecircle2a.hxx linecircle2a-HxH.hxx linecircle2a-Hxt.hxx linecircle2a-default-data.hxx linecircle2a-lsolve.hxx

  
  internal-util.h internal-util.hxx debug-util.h 

#
#  Work in progress: builds but needs more testing for production
#  cleveland/clevelandveland14a.h cleveland/clevelandveland-default.h cleveland/clevelandveland14a-default-data.h cleveland/clevelandveland14a-io.h
#  cleveland/clevelandveland14a.hxx cleveland/clevelandveland14a-default-data.hxx
#
  problem-defs.h
)
aux_source_directory(Templates minus_sources)

# this library can be used if you want to precompile instantiations
# but header-only will also work, see README.md

include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_library( minus ${minus_sources} )
set_target_properties(minus PROPERTIES LINKER_LANGUAGE CXX)
