#!/bin/bash
# Simple script for running CLDERA-E3SM with cldera-tools

# Settings
readonly CASE_NAME="v2.ULR.FIDEAL-CLDERA.eruption"
readonly COMPSET="FIDEAL"
readonly RESOLUTION="ne16pg2_ne16pg2"
#ne4pg2_oQU480
readonly MACHINE="mappy"
readonly CODE_ROOT="${HOME}/Programming/CLDERA-E3SM"
readonly CASE_ROOT="${HOME}/Programming/CLDERA-ULRFIDEAL"
readonly STOP_OPTION="ndays"
readonly STOP_N="5"
readonly CLDERA_PROFILING_CONFIG="${HOME}/Programming/cldera_profiling_config.yaml"
readonly CLDERA_PROFILING_PNETCDF="${HOME}/Programming/cldera_HSW_Tstats.h2.nc"
#/sems-data-store/ACME/cldera/data/E3SM/model_input/emissions/volc/zero-presc-volc-mmr/BMW_volcanic_1850-3009_all_zero.nc'


# Derived settings
readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_RUN_DIR=${CASE_ROOT}/run

###################
### Create case ###
###################
# Joe's case create
printf "\n\n========== CREATING CASE ==========\n"
${CODE_ROOT}/cime/scripts/create_newcase \
  --compset ${COMPSET} --res ${RESOLUTION} --case ${CASE_NAME} \
  --output-root ${CASE_ROOT} --machine ${MACHINE} --script-root ${CASE_SCRIPTS_DIR}

# My old one
#${CODE_ROOT}/cime/scripts/create_newcase \
#  --compset ${COMPSET} \
#  --case ${CASE_NAME} \
#  --res ${RESOLUTION} \
#  --machine ${MACHINE} \
#  --output-root ${CASE_ROOT} \
#  --script-root ${CASE_SCRIPTS_DIR}

##################
### Setup case ###
##################
cd ${CASE_SCRIPTS_DIR}
./xmlchange DEBUG=FALSE,DOUT_S=FALSE,STOP_OPTION=${STOP_OPTION},STOP_N=${STOP_N}
./xmlchange JOB_WALLCLOCK_TIME=$WALLCLOCK
./xmlchange SAVE_TIMING=TRUE
./xmlchange JOB_QUEUE=$QUEUE
./xmlchange RUN_STARTDATE="1991-06-14"
./xmlchange EXEROOT=${CASE_BUILD_DIR}
./xmlchange RUNDIR=${CASE_RUN_DIR}
./xmlchange --append CAM_CONFIG_OPTS="-cldera_profiling"

# Here we update the ATM model's "config options", which are a set of options different from
# the case configration options set directly by xmlchange, and different from namelist settings.
# See qualitative description here (E3SM will not agree with all details on this page):
# https://ncar.github.io/CAM/doc/build/html/users_guide/customizing-compsets.html
# 
# 'cldera_sai_trcs' is a flag which activates the idealized injection. If this is ommited, 
# the tracer fields will never be initialized
./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS --val "-cldera_sai_trcs "

printf "\n\n========== POPULATING NAMELIST SETTINGS ==========\n"
cat << EOF >> ./user_nl_eam
! --------------------------------------------------------------------------------------
!
! Joe Hollowed 3/15/23
! Namelist settings for HSW++ SAI runs. See namelist setting descriptions inline below.
! By default, this outputs 1 history file with daily averages.
!
! Namelist definitions and default values are stored at
! components/eam/bld/namelist_files/namelist_definition.xml
! components/eam/bld/namelist_files/namelist_defaults_eam.xml
!
! To see all namelist values after case creation, see
! {case_output_root}/run/atm_in
! This file is a concatenation of this contents of the present file (user_nl_eam), and
! all other defaults
!
! --------------------------------------------------------------------------------------

! ------ Output frequencies etc
empty_htapes     = .TRUE.          ! output only the varibales listed below
avgflag_pertape  = 'A'             ! hist file 1 is avg
NHTFRQ           = -24             ! output frequency every 24 hours
MFILT            = 720             ! allow 720 time samples per hist file 
NDENS            = 2               ! single-precision for each hist file
inithist         = 'ENDOFRUN'      ! output initial conditions at end of run

! ------ Output fields currently present:
! U: zonal wind in m/s
! V: meridional wind in m/s
! T: temperature in K
! OMEGA: vertical pressure velocity in Pa/s
! PS: surface pressure in Pa
! Z3: geopotential height in km
! SO2: SO2 mixing-ratio in kg/kg
! SULFATE: Sulfate mixing-ratio in kg/kg
! ASH: Ash mixing-ratio in kg/kg
! SAI_AOD: Dimensionless Aerosol Optical Depth of all species
!
! ------ Other optional fields that can be added:
! SAI_HEAT: Local tracer heating rate in K/day
! SAI_COOL: Surface heating by tracer SAI_AOD in K/day
! AIR_MASS: Mass of air in grid cells in kg (multiply through mixing ratios to get tracer masses)
! T1000: 2D temperature interpolated at 1000 hPa
! T050:  2D temperature interpolated at 50 hPa
! T025:  2D temperature interpolated at 25 hPa

fincl1           = 'U','V','T','OMEGA','PS','Z3','SO2','ASH','SULFATE','SAI_AOD'

! ------ Point to initial condition; currently 5-year HSW spinup
NCDATA="/ascldap/users/gbharpe/Programming/ens05_ic.nc"

ideal_phys_analytic_ic = .false.              ! don't let analytic ICs overwrite input from NCDATA
ideal_phys_option = 'held-suarez-williamson'  ! select HSW idealized forcing
pertlim = 0                                   ! turn off random T perturbations

! ------ Injection namelist settings
! set cldera_sai_stratHeating and cldera_sai_surfCooling both to .false. for passive tracers
cldera_sai_read_from_ic_file = .false. ! activate built-in analytic init of cldera sai tracers
cldera_sai_formSulfate = .true.        ! toggle sulfate formation
cldera_sai_stratHeating = .true.       ! toggle local heating
cldera_sai_surfCooling = .true.        ! toggle surface cooling

! injection localization
cldera_sai_t0 = 0              ! day at which to begin injection (defaults to zero)
cldera_sai_lat0 = 15.15        ! injection latitude 
cldera_sai_lon0 = 120.35       ! injection longitude

! SAI namelist settings related to tuning parameters are not set here, and will assume their 
! default values

EOF
#cat ./user_nl_eam


# My old one
#cd ${CASE_SCRIPTS_DIR}
#./xmlchange STOP_OPTION=${STOP_OPTION}
#./xmlchange STOP_N=${STOP_N}
#./xmlchange NTASKS=1 # run on only one proc
#./xmlchange run_exe='gdb ${EXEROOT}/e3sm.exe' # run gdb for debugging
#./xmlchange run_exe='${EXEROOT}/e3sm.exe'

printf "\n\n========== CASE SETUP ==========\n"
./case.setup --debug


##################
### Build case ###
##################
printf "\n\n========== BUILDING, SUBMITTING JOB ==========\n"
./case.build 2>&1 | tee ./log.case.buid

# My old one
#source .env_mach_specific.sh
#./xmlchange DEBUG=FALSE
#./case.build

###################
### Submit case ###
###################
./case.submit 2>&1 | tee ./log.case.submit
cp ${CLDERA_PROFILING_CONFIG} ${CASE_RUN_DIR}/cldera_profiling_config.yaml
cp ${CLDERA_PROFILING_PNETCDF} ${CASE_RUN_DIR}/cldera_HSW_Tstats.h2.nc
