#!/bin/bash
# Load cldera-tools modules for Graham's SEMS-based SRN workstation vaporwave

module purge
source /projects/sems/modulefiles/utils/sems-archive-modules-init.sh
module load acme-env
module load sems-archive-git
module load sems-archive-cmake/3.19.1
module load acme-gcc/8.1.0
module load acme-openmpi/4.1.4
module load acme-netcdf/4.7.4/acme

# The configure-cldera-tools-runner.sh expects PNETCDF_ROOT
export PNETCDF_ROOT=${SEMS_NETCDF_ROOT}
