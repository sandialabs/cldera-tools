# List of all options and cmake cache variables that CLDERA uses

option (CLDERA_ENABLE_TESTS "Whether to enable CLDERA tests" ON)
option (CLDERA_ENABLE_PROFILING_TOOL "Whether to build the cldera profiling tool" ON)

if (CLDERA_ENABLE_TESTS)
  # Cache vars used for testing
  set (CLDERA_TESTS_MAX_RANKS 1 CACHE STRING "Max number of ranks to use in testing")
  set (CLDERA_TESTS_MPI_EXEC_NAME mpiexec CACHE STRING "Command to be used to run MPI executables")
  set (CLDERA_TESTS_MPI_NP_FLAG -n CACHE STRING "Flag to be used to specify number of ranks")
endif()
