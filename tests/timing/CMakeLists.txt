include (EkatCreateUnitTest)

# Test timing utils
EkatCreateUnitTest (timing cldera_timing_tests.cpp
  LIBS cldera-timing
  MPI_RANKS 1 ${CLDERA_TESTS_MAX_RANKS}
)
