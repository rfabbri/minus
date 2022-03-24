# MINUS
#

cmake_minimum_required(VERSION 3.10.2 FATAL_ERROR)
set(CMAKE_CXX_STANDARD "11")
if(NOT CMAKE_CXX_STANDARD_REQUIRED)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

#------------------------------------------------------------------------------
# This command must be the last command in this file
#if(NOT MINUS_NO_EXPORT)
#  include(${CMAKE_CURRENT_LIST_DIR}/config/cmake/export/MINUSCreateProject.cmake)
#endif()

## TODO: use add_compile_options(-Wall) etc
# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall" CACHE STRING "" FORCE)

# Crucial flags for performance
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native" CACHE STRING "" FORCE)

# add_compile_options(
#     "$<$<CONFIG:RELEASE>:-march=native -ffast-math>"
#     "$<$<CONFIG:DEBUG>:-Wall>"
# )

# This one is the most important one - order of magnitude faster
# See this: https://gcc.gnu.org/ml/gcc/2001-07/msg02150.html
# if GCC 8
#   -march=native -freciprocal-math -fno-trapping-math -fno-math-errno -ffinite-math-only -fno-signaling-nans -fno-rounding-math -fcx-limited-range -fexcess-precision=fast
#   problem with gcc 8 is -fno-signed-zeros
# else
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math" CACHE STRING "" FORCE)

# --fast-math breaks into these:
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-math-errno -ffinite-math-only -fno-rounding-math -fno-signaling-nans -fcx-limited-range -fexcess-precision=fast")

# This is the --fast-math flag that does a big job but sometimes gives wrong
# results in GCC8: -funsafe-math-optimizations
# of the flags implied by this one, the most impact comes from
# -fno-signed-zeros, which is also what breaks GCC8 build
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-signed-zeros -fno-trapping-math -fassociative-math -freciprocal-math")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")
#------------------------------------------------------------------------------
project(MINUS #Project name must be all caps to have properly generated MINUS_VERSION_* variables
    VERSION 1
    DESCRIPTION "Continuation-based Minimal Problem Solver"
    LANGUAGES CXX C)
include(CMakeDependentOption)
include(${CMAKE_CURRENT_SOURCE_DIR}/config/cmake/MINUSInitializeBuildType.cmake)
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-Wno-undefined-var-template HAS_NO_UNDEFINED_VAR_TEMPLATE)
if( HAS_NO_UNDEFINED_VAR_TEMPLATE )
  add_definitions( -Wno-undefined-var-template )
endif()

#-----------------------------------------------------------------------------
if(NOT COMMAND SETIFEMPTY)
  macro(SETIFEMPTY)
    set(KEY ${ARGV0})
    set(VALUE ${ARGV1})
    if(NOT ${KEY})
      set(${ARGV})
    endif()
  endmacro()
endif()

#-----------------------------------------------------------------------------
SETIFEMPTY(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

#-----------------------------------------------------------------------------
SETIFEMPTY(CMAKE_INSTALL_LIBRARY_DESTINATION lib)
SETIFEMPTY(CMAKE_INSTALL_ARCHIVE_DESTINATION lib)
SETIFEMPTY(CMAKE_INSTALL_RUNTIME_DESTINATION bin)

# Allow external project to override the export target
# if(NOT COMP_NO_EXPORT)
#   SETIFEMPTY(COMP_INSTALL_EXPORT_NAME COMPTargets)
# endif()

SETIFEMPTY(MINUS_INSTALL_RUNTIME_DIR bin)
SETIFEMPTY(MINUS_INSTALL_LIBRARY_DIR lib)
SETIFEMPTY(MINUS_INSTALL_ARCHIVE_DIR lib)
SETIFEMPTY(MINUS_INSTALL_INCLUDE_DIR include/)
if(NOT MINUS_LIB_PREFIX)
  set( MINUS_LIB_PREFIX "") # This is typically empty
endif()

# CMake support directory.
set(MINUS_ROOT_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})
set(MINUS_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR}/config/cmake/Modules)

include(${MINUS_CMAKE_DIR}/MINUSStandardOptions.cmake)
include( GenerateExportHeader )
#include(${CMAKE_CURRENT_LIST_DIR}/config/cmake/doxygen/doxygen.cmake)
# Only include these if we have vxd-specific cmake utilities. We currently use VXL's:
#include(${CMAKE_CURRENT_LIST_DIR}/config/cmake/config/MINUS_utils.cmake)

# Location of VXD's FindXXX.cmake CMake modules.
# This is identical to MINUS_CMAKE_DIR.  Perhaps we should eliminate MODULE_PATH?
set( MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/config/cmake/Modules CACHE STATIC "MINUS module path" )

# Options to add extra compiler and linker flags
#
# These options allow you to specify additional flags without
# affecting the default flags for a particular platform or build type.
# This is especially useful for adding extra warning flags.
set( MINUS_EXTRA_CMAKE_C_FLAGS             CACHE STRING "Extra flags appended to CMAKE_C_FLAGS" )
set( MINUS_EXTRA_CMAKE_CXX_FLAGS           CACHE STRING "Extra flags appended to CMAKE_CXX_FLAGS" )
set( MINUS_EXTRA_CMAKE_EXE_LINKER_FLAGS    CACHE STRING "Extra flags appended to CMAKE_EXE_LINKER_FLAGS" )
set( MINUS_EXTRA_CMAKE_MODULE_LINKER_FLAGS CACHE STRING "Extra flags appended to CMAKE_MODULE_LINKER_FLAGS" )
set( MINUS_EXTRA_CMAKE_SHARED_LINKER_FLAGS CACHE STRING "Extra flags appended to CMAKE_SHARED_LINKER_FLAGS" )
set( CMAKE_C_FLAGS             "${CMAKE_C_FLAGS} ${MINUS_EXTRA_CMAKE_C_FLAGS}" )
set( CMAKE_CXX_FLAGS           "${CMAKE_CXX_FLAGS} ${MINUS_EXTRA_CMAKE_CXX_FLAGS}" )
set( CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} ${MINUS_EXTRA_CMAKE_EXE_LINKER_FLAGS}" )
set( CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${MINUS_EXTRA_CMAKE_MODULE_LINKER_FLAGS}" )
set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${MINUS_EXTRA_CMAKE_SHARED_LINKER_FLAGS}" )

#-------------------------------------------------------------------
#-- BUILD CONFIG OPTIONS

# In order to be able to link vxl libraries into shared libraries on 64 bit linux, the -fPIC
# compiler flag must be added.  Only do this if we are on a x86_64 *nix platform, we're building
# static libraries, and the user has not explicitly requested position dependent code.
if(UNIX)
  if(NOT BUILD_SHARED_LIBS AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
    option(VXL_BUILD_POSITION_DEPENDENT_CODE "Generate position dependent code (i.e. code cannot be used in shared library)" OFF)
    mark_as_advanced(VXL_BUILD_POSITION_DEPENDENT_CODE)
    if(NOT VXL_BUILD_POSITION_DEPENDENT_CODE)
      message(STATUS "Adding -fPIC compiler flag to generate position independent code.")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
    endif()
  endif()
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # require at least gcc 4.8
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
      message(WARNING "***********************************")
      message(WARNING "GCC must be at least 5.0 for speed!")
      message(WARNING "*********************************")
    elseif(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 7.9)
      message(WARNING "******************************************\nGCC must be less than 8.0 for reliability!\n******************************************\n")
    endif()
elseif (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
else()
    message(WARNING "You are using an unsupported compiler! Compilation has only been tested with Clang and GCC.")
endif()


# Some types of path names can cause havoc with regular expressions,
# so avoid those.
if( ${PROJECT_BINARY_DIR} MATCHES  \\+ )
  message(SEND_ERROR "You cannot have a + in your binary path")
endif()
if( ${CMAKE_CURRENT_LIST_DIR} MATCHES  \\+ )
  message(SEND_ERROR "You cannot have a + in your source path")
endif()

# include CMakeListsLocal.txt from source directory if it exists
# also include it from the binary dir if different from source dir
if( ${PROJECT_BINARY_DIR} MATCHES ${CMAKE_CURRENT_LIST_DIR} )
  include( ${CMAKE_CURRENT_LIST_DIR}/CMakeListsLocal.txt OPTIONAL )
else()
  include( ${CMAKE_CURRENT_LIST_DIR}/CMakeListsLocal.txt OPTIONAL )
  include( ${PROJECT_BINARY_DIR}/CMakeListsLocal.txt OPTIONAL )
endif()

set(MINUS_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR})    #Source includes
set(MINUS_BINARY_INCLUDE_DIR ${PROJECT_BINARY_DIR}) #Generated includes
set(MINUS_INCLUDE_DIRS ${MINUS_BINARY_INCLUDE_DIR} ${MINUS_INCLUDE_DIR})
set(MINUS_INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include/)

include_directories(${MINUS_INCLUDE_DIRS})


# Optinonally look for VXL if you want to built tests
#
option(USE_VXL "Use VXL libraries" NO)
if (USE_VXL)
  find_package(VXL PATHS "../vxl-bin")
  if (VXL_FOUND)
    # include the standard stuff, such as Dart test
    set(VXL_PROVIDE_STANDARD_OPTIONS 1)
    include(${VXL_CMAKE_DIR}/UseVXL.cmake)
  else ()
    message( STATUS "VXL not found - OK, only needed if building tests" )
  endif()
endif()


if(VXL_FOUND)
  # This is required on some platforms because various libraries
  # include OpenGL indirectly even when it is not used
  include_directories( ${VXL_VGUI_INCLUDE_DIR} )
endif()

option(USE_Boost "Use Boost libraries" NO)
if ( USE_Boost )
  set(Boost_USE_STATIC_LIBS   OFF)
  find_package( Boost 1.35.0 COMPONENTS iostreams)
  #if(Boost_FOUND)
  #    message("aww yeah")
  # else(Boost_FOUND)
  #  message("OH NO")
  #endif(Boost_FOUND)
  set( Boost_FOUND TRUE )
else ( USE_Boost )
  set( Boost_FOUND FALSE )
endif( USE_Boost )

# adds a test for Endianness and a global variable that should be useful for file format programming -MM
include(TestBigEndian)
TEST_BIG_ENDIAN(MY_BIG_ENDIAN)

# For use in client projects that use CMake and not installed version
set(MINUS_LIBRARY_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} CACHE STATIC "Where all the minus libraries are, for clients to use." )

#-------------------------------------------------------------------------------
# main lib
add_subdirectory(minus)
# commands / executables
add_subdirectory(cmd) 

set(BUILD_TESTING OFF CACHE BOOL "Build tests for minus")
set(VXL_BUILD_POSITION_DEPENDENT_CODE ON CACHE BOOL "Generate position dependent code - code cannot be used in shared library")
  
if(VXL_FOUND)  # For now we need vxl to build tests
  if( BUILD_TESTING )
    add_subdirectory(tests)
  endif()
  add_subdirectory(algo)
endif()


# Good to have / secondary -------
#
# vectorize - used by eigen casually, reportedly giving 1.3x speedup in eg
# matrix multiplication, but no gain seen here so far
# already enabled by -march=native mostly. If you enable one of these and the
# arch is does not support this well, the software may crash with invalid
# instructions
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mavx -mfma")

# Use naive and fast complex multiplication without checking for corner cases
# most compiler will have this as part of -ffast-math
# later on, we can use only this flag if equivalent gains to the more aggressive
# -fast-math
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fcx-limited-range")

# For comparison purposes
# find_package(BLAS)
# find_package(LAPACK)
# if(LAPACK_FOUND AND BLAS_FOUND)
#   message("FOUND YEAH" ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
#   target_link_libraries( bnld_algo bnld vnl_algo vil brent ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
# else()
#   message("FOUND NO")
#   target_link_libraries( bnld_algo bnld vnl_algo vil brent)
# endif()
