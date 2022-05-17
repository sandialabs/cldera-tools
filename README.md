## CLDERA TOOLS

NOTE: This README is a copy of [the confluence page](https://sems-atlassian-son.sandia.gov/confluence/display/CSPT/CLDERA+profiling+tool) intended more to keep relevant information with the software. For the most up-to-date information, see the confluence page.

The cldera-tools repo contains a library whose goal is to "profile" quantities in E3SM. As of now (5/17/22), the library allows to grab fields pointers from E3SM, compute some basic global "statistics" on them (max, min, sum, avg) and save their values to an output YAML file. The names of the fields to track, as well as the statistics to track, are specified via an input YAML file. Right now, the location/name of the file is hardcoded: must be called `cldera_profiling_config.yaml` and must be located in the run folder of the CIME case. For example, consider the following input file

```
%YAML 1.0
---
Stats Output File: cldera_stats.yaml
Fields To Track: [SO2, H2SO4]

SO2:
  Compute Stats: [Max]
H2SO4:
  Compute Stats: [Min]
...
```

Running an E3SM case that uses CLDERA (e.g., `F20TR-CLDERA`) would generate the file `cldera_stats.yaml` in the run folder, whose content will look like the following:

```
%YAML 1.1
---
H2SO4:
  Min:
    Timestamps: [19910101.3600, 19910101.7200, 19910101.10800, 19910101.14400, 19910101.18000, 19910101.21600, 19910101.25200, 19910101.28800, 19910101.32400, 19910101.36000, 19910101.39600, 19910101.43200, 19910101.46800, 19910101.50400, 19910101.54000, 19910101.57600, 19910101.61200, 19910101.64800, 19910101.68400, 19910101.72000, 19910101.75600, 19910101.79200, 19910101.82800, 19910102.0]
    Values: [7.8533024858398e-27, 2.5668891775532e-26, 1.2441217677756e-27, 6.5550604409998e-26, 7.0864413176058e-25, 1.6922991276391e-23, 4.1074848003180e-23, 1.3154042557344e-22, 1.3350997795038e-22, 5.7456865021541e-23, 5.9686624506205e-23, 2.2416232474606e-23, 1.1999704857521e-24, 2.4560533800177e-29, 1.2228147799198e-23, 7.2223714970792e-22, 5.5037631035254e-25, 6.6304561540602e-26, 2.1858451717351e-25, 4.8826322149140e-25, 1.1758402879695e-22, 8.3975348678835e-23, 4.3280576584943e-21, 6.2645093521745e-23]
SO2:
  Max:
    Timestamps: [19910101.3600, 19910101.7200, 19910101.10800, 19910101.14400, 19910101.18000, 19910101.21600, 19910101.25200, 19910101.28800, 19910101.32400, 19910101.36000, 19910101.39600, 19910101.43200, 19910101.46800, 19910101.50400, 19910101.54000, 19910101.57600, 19910101.61200, 19910101.64800, 19910101.68400, 19910101.72000, 19910101.75600, 19910101.79200, 19910101.82800, 19910102.0]
    Values: [3.6087438125778e-08, 3.5894425708473e-08, 3.5721556136599e-08, 3.5483199956125e-08, 3.5253918091705e-08, 3.6059542506238e-08, 3.6836519271650e-08, 3.8432711743365e-08, 3.9137532803743e-08, 4.0995685281938e-08, 4.1824806810879e-08, 4.3649447241393e-08, 4.5322348063647e-08, 4.7589193008743e-08, 4.5886412540401e-08, 4.7434266908599e-08, 4.6544989213327e-08, 4.7089547608924e-08, 4.3591178610284e-08, 4.4747558792056e-08, 4.6106373923549e-08, 4.4921964309152e-08, 4.5909736907053e-08, 4.6846372141518e-08]
...
```

The output file contains a sublist for each field; for each field, there's a sublist for each requested statistic; for each statistic, CLDERA saves the timestamps where the corresponding values were obtained (timestamp is in the format YYYYMMDD.TOD, where TOD is the time-of-day in seconds).

Currently, one can specify in the input file any field that is a "tracer" in EAM and/or stored in "pbuf". This is ruling out state variables (u,dp,T), but it should be easy to add support for those. For each field, one can ask for more than one statistic.

## Where/how to run

Currently, cldera-tools is installed on mappy. The cldera-e3sm fork of E3SM has the location of this pre-built library, so that it can be linked at build time if a CLDERA compset is enabled. To add a cldera compset, you need to edit the `config_compset.xml` file in `eam/cime_config`. Currently, the only compset added is the following:

```
  <compset>
    <alias>F20TR-CLDERA</alias>
    <lname>20TR_EAM%CMIP6-CLDERA_ELM%SPBC_MPASSI%PRES_DOCN%DOM_SROF_SGLC_SWAV</lname>
  </compset>
```

Which runs active atm (EAM) along with active sea-ice (mpas-seaice), and data for the other components. So long as your compset name contains the string "CLDERA", the cldera-profiling tool will be enabled. If you want to enable it for existing compsets, you need to edit `config_component.xml` in `eam/cime_config`, adding `-cldera-profiling` to the `CAM_CONFIG_OPTS` xml entry.
Once you created the compset, you simply create a new case:

```
$ /path/to/cime/scripts/create_newcase --compset YOUR_COMPSET --case CASE_NAME --res GRID --machine mappy
```

Here, GRID is a grid that works with your compset. For testing the tool with the `F20TR-CLDERA` compset, I used `ne4pg2_oQU480`. But you probably will want to use higher resolutions (`query-config --grid` can help to pick the right grid for your compset). After that, it's the usual

```
$ cd CASE_NAME
$ ./case.setup
$ ./case.buid
$ ./case.submit
```

Before the submit step, makes sure you copy your `cldera_profiling_config.yaml` file in the run directory. You can find the path to the run dir by doing

```
$ ./xmlquery RUNDIR
```

from the case dir (xmlquery is a very useful tool to check your config). You can also change settings via ./xmlchange . Most notably, I think the default value for DEBUG is TRUE, so you might want to do ./xmlchange DEBUG=FALSE (you'll have to do it before the build phase, or CIME will inform you that you need to clean and rebuild).

## Build your own version of cldera-tools

If you start developing in cldera tools, you might be interested in building a local branch, and use it in E3SM. To do so, you can use the following configuration script (I'll call it `do-cmake.sh`)

```
SRC_DIR=/path/to/source/tree
INSTALL_DIR=/path/to/install/folder

rm -rf CMakeFiles
rm -f  CMakeCache.txt

cmake \
  \
  -D CMAKE_BUILD_TYPE:STRING=DEBUG            \
  -D CMAKE_C_COMPILER:STRING=mpicc            \
  -D CMAKE_CXX_COMPILER:STRING=mpicxx         \
  -D CMAKE_Fortran_COMPILER:STRING=mpifort    \
  -D CMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}/debug \
  \
  -D CLDERA_ENABLE_TESTS:BOOL=ON              \
  -D CLDERA_ENABLE_PROFILING_TOOL:BOOL=ON     \
  -D CLDERA_TESTS_MAX_RANKS:STRING=4          \
  \
  ${SRC_DIR}
```

Before running the script, be sure your environment matches that of a typical E3SM run, to avoid any potential conflict in library versions. To do so, create a case (as described above), and run case.setup . In the case folder, you should see a file called `.env_mach_specific.sh`. Source that file to get the E3SM-compatible environment. Notice that the software is installed in `${INSTALL_DIR}/debug`. If you build in release mode, install in `${INSTALL_DIR}/release`.

Then,

```
$ ./do-cmake.sh
$ make -j
$ make install
```

Once cldera tools is installed, go in the source tree of E3SM that you will build, and open the file `cime_config/machines/cmake_macros/$compiler_$mach.cmake`, where `$compiler` and `$cmake` are the names of the compiler and machine you are using (e.g., `gnu_mappy.cmake` when running on mappy), and edit the variable `CLDERA_PATH`  to point to the what you used for `INSTALL_DIR`  in the `do-cmake.sh` script. Notice that `CLDERA_PATH`  should be whatever comes before `debug` or `release`. CIME will then add `debug` or `release` depending on whether E3SM is built with DEBUG on or off.

NOTE: you need to make sure the submodules are up to date in the source folder. If you cloned the repo with `git clone --recursive`, they should already be. Otherwise, you need to issue `git submodule update --init --recursive`  in the source folder.
