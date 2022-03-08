include (EkatCreateUnitTest)

# Test ProfilingSession
EkatCreateUnitTest (session session.cpp
  LIBS cldera-profiling ekat
  MPI_RANKS 1 ${CLDERA_TESTS_MAX_RANKS}
  MPI_EXEC_NAME ${CLDERA_TESTS_MPI_EXEC_NAME}
  MPI_NP_FLAG ${CLDERA_TESTS_MPI_NP_FLAG})

# Test ProfilingArchive
EkatCreateUnitTest (archive archive.cpp
  LIBS cldera-profiling ekat)

# Test functions to compute stats
EkatCreateUnitTest (compute_stats compute_stats.cpp
  LIBS cldera-profiling ekat)
