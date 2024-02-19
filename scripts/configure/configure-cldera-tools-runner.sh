#!/bin/bash
# A script to configure the project. 
# This is primarily used used by the gitlab runner,
# although it should work universally.
# Usage:
#   ./configure-cldera-tools-runner SOURCE
# where SOURCE is the source directory for cldera-tools.
# This should only be invoked from the build directory, 
# as it sets the build directory to the current directory.

SRC_DIR="$@"
INSTALL_DIR="$(pwd)/install"
 
rm -rf CMakeFiles
rm -f  CMakeCache.txt
 
cmake \
  -D CMAKE_BUILD_TYPE:STRING=DEBUG            \
  -D CMAKE_C_COMPILER:STRING=mpicc            \
  -D CMAKE_CXX_COMPILER:STRING=mpicxx         \
  -D CMAKE_Fortran_COMPILER:STRING=mpifort    \
  -D CMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} \
  -D CLDERA_PNETCDF_PATH:PATH=${PNETCDF_ROOT} \
  -D CLDERA_ENABLE_TESTS:BOOL=ON              \
  -D CLDERA_ENABLE_PROFILING_TOOL:BOOL=ON     \
  -D CLDERA_TESTS_MAX_RANKS:STRING=4          \
  ${SRC_DIR}
