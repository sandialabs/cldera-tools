# This file is for user convenience only and is not used by the model
# Changes to this file will be ignored and overwritten
# Changes to the environment should be made in env_mach_specific.xml
# Run ./case.setup --reset to regenerate this file
source /projects/sems/install/rhel7-x86_64/sems/v2/lmod/lmod/8.3/gcc/10.1.0/zbzzu7k/lmod/lmod/init/sh
module purge 
module load sems-archive-env acme-env sems-archive-git sems-archive-cmake/3.19.1 acme-gcc/8.1.0 acme-openmpi/4.1.4 acme-netcdf/4.7.4/acme
export CLDERA_PATH=/sems-data-store/ACME/cldera/cldera-tools/install
export NETCDFROOT=/projects/sems/install/rhel7-x86_64/acme/tpl/netcdf/4.7.4/gcc/8.1.0/openmpi/4.1.4/acme
export OMP_STACKSIZE=64M
export OMP_PROC_BIND=spread
export OMP_PLACES=threads