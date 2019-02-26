# This VXDCreateProject.cmake file handles the creation of files needed by
# other client projects that use VXD.  Nothing is built by this
# CMakeLists.txt file.  This CMakeLists.txt file must be processed by
# CMake after all the other CMakeLists.txt files in the VXD tree,
# which is why the add_subdirectory(config/cmake/export) command is at the end
# of the top level CMakeLists.txt file.

# Save library dependencies.
#set(MINUS_CMAKE_DOXYGEN_DIR  ${MINUS_ROOT_SOURCE_DIR}/config/cmake/doxygen)

get_property(VXDTargets_MODULES GLOBAL PROPERTY VXDTargets_MODULES)

set(MINUS_CONFIG_CMAKE_DIR "share/vxd/cmake")
if(${CMAKE_VERSION} VERSION_LESS 2.8.12)
   set(INTERFACE_LINK_OPTION "")
else()
   set(INTERFACE_LINK_OPTION "EXPORT_LINK_INTERFACE_LIBRARIES")
endif()

if(VXDTargets_MODULES)
  export(TARGETS
    ${VXDTargets_MODULES}
    APPEND
    FILE "${CMAKE_CURRENT_BINARY_DIR}/VXDTargets.cmake"
    ${INTERFACE_LINK_OPTION}
  )
  install(EXPORT ${MINUS_INSTALL_EXPORT_NAME} DESTINATION ${MINUS_CONFIG_CMAKE_DIR}
          COMPONENT Development)
endif()

# Create the VXDConfig.cmake file for the build tree.
configure_file(${MINUS_CMAKE_DIR}/VXDConfig.cmake.in
               ${PROJECT_BINARY_DIR}/VXDConfig.cmake @ONLY)

configure_file(${MINUS_CMAKE_DIR}/VXDConfig_export.cmake.in
               ${PROJECT_BINARY_DIR}/config/cmake/export/VXDConfig.cmake
               @ONLY)

install(FILES
  ${PROJECT_BINARY_DIR}/config/cmake/export/VXDConfig.cmake
  ${MINUS_CMAKE_DIR}/VXDStandardOptions.cmake
  ${MINUS_CMAKE_DIR}/UseVXD.cmake
  DESTINATION ${MINUS_CONFIG_CMAKE_DIR}
)
