cmake_minimum_required(VERSION 3.5)

project(libUnmix LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(LIBUNMIX_MAJOR_VERSION 1)
set(LIBUNMIX_MINOR_VERSION 0)
set(LIBUNMIX_PATCH_VERSION 0)
set(LIBUNMIX_VERSION
  ${LIBUNMIX_MAJOR_VERSION}.${LIBUNMIX_MINOR_VERSION}.${LIBUNMIX_PATCH_VERSION})

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_INCLUDE_DIR include CACHE PATH
  "Installation directory for header files")
set(DEF_INSTALL_CMAKE_DIR CMake)

set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${PROJECT_BINARY_DIR}/${${var}}")
  endif()
endforeach()

message("Install include dir : ${INSTALL_INCLUDE_DIR}")
message("Install lib dir : ${INSTALL_LIB_DIR}")

include_directories(
  "${PROJECT_SOURCE_DIR}/include")

add_subdirectory(libUnmix)

# Add all targets to the build-tree export set
export(TARGETS libUnmix
  FILE "${PROJECT_BINARY_DIR}/libUnmixTargets.cmake")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE libUnmix)

# Create the libUnmixConfig.cmake and libUnmixConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}")
configure_file(libUnmixConfig.cmake.in
  "${PROJECT_BINARY_DIR}/libUnmixConfig.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "/${LIBUNMIX_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(libUnmixConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/libUnmixConfig.cmake" @ONLY)
# ... for both
configure_file(libUnmixConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/libUnmixConfigVersion.cmake" @ONLY)

# Install the libUnmixConfig.cmake and libUnmixConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/libUnmixConfig.cmake"
  "${PROJECT_BINARY_DIR}/libUnmixConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)

# Install the export set for use with the install-tree
install(EXPORT libUnmixTargets DESTINATION
  "${INSTALL_CMAKE_DIR}" COMPONENT dev)

