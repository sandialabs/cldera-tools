# This file is for user convenience only and is not used by the model
# Changes to this file will be ignored and overwritten
# Changes to the environment should be made in env_mach_specific.xml
# Run ./case.setup --reset to regenerate this file
source /projects/sems/modulefiles/utils/sems-archive-modules-init.sh
source /usr/share/lmod/lmod/init/sh
module purge 
module load sems-env acme-env sems-git sems-cmake/3.19.1 gnu/6.3.1 sems-intel/17.0.0 sems-openmpi/1.10.5 acme-netcdf/4.7.4/acme
export NETCDFROOT=/projects/sems/install/toss3/acme/tpl/netcdf/4.7.4/intel/17.0.0/openmpi/1.10.5/acme
export NETCDF_INCLUDES=/projects/sems/install/toss3/acme/tpl/netcdf/4.7.4/intel/17.0.0/openmpi/1.10.5/acme/include
export NETCDF_LIBS=/projects/sems/install/toss3/acme/tpl/netcdf/4.7.4/intel/17.0.0/openmpi/1.10.5/acme/lib
export OMP_STACKSIZE=64M
export PNETCDFROOT=/projects/sems/install/toss3/acme/tpl/netcdf/4.7.4/intel/17.0.0/openmpi/1.10.5/acme
