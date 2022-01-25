include (EkatCreateUnitTest)

# Create a test
EkatCreateUnitTest (profiling profiling.cpp
  LIBS cldera-profiling ekat
  MPI_RANKS 1 ${CLDERA_TESTS_MAX_RANKS}
  MPI_EXEC_NAME ${CLDERA_TESTS_MPI_EXEC_NAME}
  MPI_NP_FLAG ${CLDERA_TESTS_MPI_NP_FLAG})
