#!/bin/bash
# Load cldera-tools modules for Luca's workstation

source /projects/sems/modulefiles/utils/sems-modules-init.sh
module purge
module load sems-git sems-cmake/3.21.1 sems-python/3.9.0 sems-gcc/10.1.0 sems-openmpi/4.0.5 sems-netcdf-c/4.7.4 sems-netcdf-fortran sems-parallel-netcdf

# The configure-cldera-tools-runner.sh expects PNETCDF_ROOT
export PNETCDF_ROOT=${PARALLEL_NETCDF_ROOT}
