# COMPS
#

cmake_minimum_required(VERSION 3.10.2 FATAL_ERROR)
set(CMAKE_CXX_STANDARD "11")
if(NOT CMAKE_CXX_STANDARD_REQUIRED)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

project(COMPS #Project name must be all caps to have properly generated COMPS_VERSION_* variables
    VERSION 1
    DESCRIPTION "Continuation-based Minimal Problem Solver"
    LANGUAGES CXX C)
include(CMakeDependentOption)
include(${CMAKE_CURRENT_SOURCE_DIR}/config/cmake/VXDInitializeBuildType.cmake)
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

SETIFEMPTY(COMPS_INSTALL_RUNTIME_DIR bin)
SETIFEMPTY(COMPS_INSTALL_LIBRARY_DIR lib)
SETIFEMPTY(COMPS_INSTALL_ARCHIVE_DIR lib)
SETIFEMPTY(COMPS_INSTALL_INCLUDE_DIR include/vxd)
if(NOT COMPS_LIB_PREFIX)
  set( COMPS_LIB_PREFIX "") # This is typically empty
endif()

# CMake support directory.
set(COMPS_ROOT_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})
set(COMPS_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR}/config/cmake/Modules)

include(${COMPS_CMAKE_DIR}/COMPSStandardOptions.cmake)
include( GenerateExportHeader )
#include(${CMAKE_CURRENT_LIST_DIR}/config/cmake/doxygen/doxygen.cmake)
# Only include these if we have vxd-specific cmake utilities. We currently use VXL's:
#include(${CMAKE_CURRENT_LIST_DIR}/config/cmake/config/COMPS_utils.cmake)

# Location of VXD's FindXXX.cmake CMake modules.
# This is identical to COMPS_CMAKE_DIR.  Perhaps we should eliminate MODULE_PATH?
set( MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/config/cmake/Modules CACHE STATIC "COMPS module path" )

# Options to add extra compiler and linker flags
#
# These options allow you to specify additional flags without
# affecting the default flags for a particular platform or build type.
# This is especially useful for adding extra warning flags.
set( COMPS_EXTRA_CMAKE_C_FLAGS             CACHE STRING "Extra flags appended to CMAKE_C_FLAGS" )
set( COMPS_EXTRA_CMAKE_CXX_FLAGS           CACHE STRING "Extra flags appended to CMAKE_CXX_FLAGS" )
set( COMPS_EXTRA_CMAKE_EXE_LINKER_FLAGS    CACHE STRING "Extra flags appended to CMAKE_EXE_LINKER_FLAGS" )
set( COMPS_EXTRA_CMAKE_MODULE_LINKER_FLAGS CACHE STRING "Extra flags appended to CMAKE_MODULE_LINKER_FLAGS" )
set( COMPS_EXTRA_CMAKE_SHARED_LINKER_FLAGS CACHE STRING "Extra flags appended to CMAKE_SHARED_LINKER_FLAGS" )
set( CMAKE_C_FLAGS             "${CMAKE_C_FLAGS} ${COMPS_EXTRA_CMAKE_C_FLAGS}" )
set( CMAKE_CXX_FLAGS           "${CMAKE_CXX_FLAGS} ${COMPS_EXTRA_CMAKE_CXX_FLAGS}" )
set( CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} ${COMPS_EXTRA_CMAKE_EXE_LINKER_FLAGS}" )
set( CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${COMPS_EXTRA_CMAKE_MODULE_LINKER_FLAGS}" )
set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${COMPS_EXTRA_CMAKE_SHARED_LINKER_FLAGS}" )

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

add_subdirectory(cmd) 

if(VXL_FOUND)
  # This is required on some platforms because various libraries
  # include OpenGL indirectly even when it is not used
  include_directories( ${VXL_VGUI_INCLUDE_DIR} )
  # xxx  tests
  add_subdirectory(tests)
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
set(COMPS_LIBRARY_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} CACHE STATIC "Where all the vxd libraries are, for clients to use." )

# This command must be the last command in this file
if(NOT COMPS_NO_EXPORT)
  include(${CMAKE_CURRENT_LIST_DIR}/config/cmake/export/COMPSCreateProject.cmake)
endif()
