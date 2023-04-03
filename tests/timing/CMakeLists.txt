include (EkatCreateUnitTest)

# Test timing utils
EkatCreateUnitTest (timing cldera_timing_tests.cpp
  LIBS cldera-timing
  MPI_RANKS 1 ${CLDERA_TESTS_MAX_RANKS}
  MPI_EXEC_NAME ${CLDERA_TESTS_MPI_EXEC_NAME}
  MPI_NP_FLAG ${CLDERA_TESTS_MPI_NP_FLAG}
)
