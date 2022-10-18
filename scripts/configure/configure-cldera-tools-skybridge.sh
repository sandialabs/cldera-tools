#!/bin/bash

SRC_DIR=/ascldap/users/gbharpe/hpc_home/CLDERA/cldera-tools-source
INSTALL_DIR=/ascldap/users/gbharpe/hpc_home/CLDERA/cldera-tools-build

 
rm -rf CMakeFiles
rm -f  CMakeCache.txt
 
cmake \
  -D CMAKE_BUILD_TYPE:STRING=DEBUG            \
  -D CMAKE_C_COMPILER:STRING=mpicc            \
  -D CMAKE_CXX_COMPILER:STRING=mpicxx         \
  -D CMAKE_Fortran_COMPILER:STRING=mpifort    \
  -D CMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}/release \
  -D CLDERA_PNETCDF_PATH:PATH=${SEMS_NETCDF_ROOT}  \
  -D CLDERA_ENABLE_TESTS:BOOL=ON              \
  -D CLDERA_ENABLE_PROFILING_TOOL:BOOL=ON     \
  -D CLDERA_TESTS_MAX_RANKS:STRING=4          \
  ${SRC_DIR}
