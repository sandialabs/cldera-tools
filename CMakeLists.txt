# Require a decently recent CMake...
cmake_minimum_required (VERSION 3.18)

# PROJECT is mainly in C++, with C hooks for C/Fortran
project (CLDERA C CXX Fortran)

set (CLDERA_VERSION_MAJOR 0)
set (CLDERA_VERSION_MINOR 1)
set (CLDERA_VERSION_PATCH 0)

##############################
#        Build EKAT          #
##############################
set (EKAT_SRC_DIR ${CLDERA_SOURCE_DIR}/externals/ekat)
set (EKAT_BLD_DIR ${CLDERA_BINARY_DIR}/externals/ekat)

set (EKAT_ENABLE_TESTS FALSE)
set (EKAT_DISABLE_TPL_WARNINGS TRUE)

# Make CMake honor the options we set
set (CMAKE_POLICY_DEFAULT_CMP0077 NEW)
cmake_policy(SET CMP0077 NEW)

add_subdirectory(${EKAT_SRC_DIR} ${EKAT_BLD_DIR})

##############################
#       Build CLDERA         #
##############################

set (EKAT_CMAKE_PATH ${EKAT_SRC_DIR}/cmake CACHE INTERNAL "")
list(APPEND CMAKE_MODULE_PATH ${EKAT_CMAKE_PATH})
include(EkatMpiUtils)
DisableMpiCxxBindings()

# Add search paths for cmake modules
list (APPEND CMAKE_MODULE_PATH
      ${CMAKE_CURRENT_SOURCE_DIR}/cmake
      ${CMAKE_CURRENT_SOURCE_DIR}/externals/ekat/cmake)

# Get all CLDERA cmake config options
include (cldera_options)

#  Use EKAT's utils to set a default for CLDERA_DEBUG
include (EkatUtils)
IsDebugBuild(CLDERA_DEBUG_DEFAULT)
option (CLDERA_DEBUG "Turn on some debug features in CLDERA" ${CLDERA_DEBUG_DEFAULT})

# Libraries/tools
add_subdirectory (src)

# Tests
if (CLDERA_ENABLE_TESTS)
  include (CTest)

  # System tests
  add_subdirectory (tests)
endif()

###########################################
###   Package CLDERA as a CMake package   ###
###########################################

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

configure_package_config_file(
  cmake/ClderaConfig.cmake.in
  "${CLDERA_BINARY_DIR}/ClderaConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}/cmake)

write_basic_package_version_file(
  "${CLDERA_BINARY_DIR}/ClderaConfigVersion.cmake"
  VERSION ${CLDERA_VERSION_MAJOR}.${CLDERA_VERSION_MINOR}.${CLDERA_VERSION_PATCH}.
  COMPATIBILITY SameMajorVersion)


# Install the ClderaConfig*.cmake files
install(FILES
  "${CLDERA_BINARY_DIR}/ClderaConfig.cmake"
  "${CLDERA_BINARY_DIR}/ClderaConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cldera)

# Install cmake targets
install(EXPORT ClderaTargets
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cldera)

# Install the template cldera_profiling_config.yaml
install(FILES
  "${CLDERA_SOURCE_DIR}/scripts/cldera_profiling_config.yaml"
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR})
