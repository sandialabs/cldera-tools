include (EkatCreateUnitTest)

# Test ProfilingContext
EkatCreateUnitTest (control control.cpp
  LIBS cldera-control ekat
  MPI_RANKS 1 ${CLDERA_TESTS_MAX_RANKS}
)

# Set yaml files for tests that utilize input decks
set (CLDERA_CONTROL_TEST_YAMLS
  cldera_control_input.yaml)

# Copy the yaml files to the destination
FILE(COPY ${CLDERA_CONTROL_TEST_YAMLS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
