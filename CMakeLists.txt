# Require a decently recent CMake...
cmake_minimum_required (VERSION 3.18)

# PROJECT is mainly in C++, with C hooks for C/Fortran
project (CLDERA C CXX)

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

# Add search paths for cmake modules
list (APPEND CMAKE_MODULE_PATH
      ${CMAKE_CURRENT_SOURCE_DIR}/cmake
      ${CMAKE_CURRENT_SOURCE_DIR}/externals/ekat/cmake)

# Get all CLDERA cmake config options
include (cldera_options)

#  Use EKAT's utils to set CLDERA_DEBUG
include (EkatUtils)
IsDebugBuild(CLDERA_DEBUG)

# Libraries/tools
add_subdirectory (src)

# Tests
if (CLDERA_ENABLE_TESTS)
  include (CTest)

  # System tests
  add_subdirectory (tests)
endif()
