include (EkatCreateUnitTest)

# Test ProfilingSession
EkatCreateUnitTest (session session.cpp
  LIBS cldera-profiling ekat
  MPI_RANKS 1 ${CLDERA_TESTS_MAX_RANKS}
  MPI_EXEC_NAME ${CLDERA_TESTS_MPI_EXEC_NAME}
  MPI_NP_FLAG ${CLDERA_TESTS_MPI_NP_FLAG})

# Test Field
EkatCreateUnitTest (time_stamp time_stamp.cpp
  LIBS cldera-profiling ekat)

# Test Field
EkatCreateUnitTest (field field.cpp
  LIBS cldera-profiling ekat)

# Test ProfilingArchive
EkatCreateUnitTest (archive archive.cpp
  LIBS cldera-profiling ekat
  PROPERTIES FIXTURES_SETUP archive_output)

add_test (NAME archive_check
          COMMAND ncdump -v foo_global_max archive_tests.nc)
set_tests_properties(archive_check PROPERTIES
  PASS_REGULAR_EXPRESSION "foo_global_max = 3"
  FIXTURES_REQUIRED archive_output
)

# Test Graph
EkatCreateUnitTest (graph graph.cpp
  LIBS cldera-profiling ekat)

# Test functions to compute stats
EkatCreateUnitTest (compute_stats compute_stats.cpp
  LIBS cldera-profiling ekat)

# Test functions to compute stats (requires cldera-pnetcdf)
if (CLDERA_PNETCDF_PATH)
  EkatCreateUnitTest (compute_stats_pnetcdf compute_stats_pnetcdf.cpp
    LIBS cldera-profiling cldera-pnetcdf ekat)
endif()

# Test Stats Field Test
EkatCreateUnitTest (bounds_field_test bounds_field_test.cpp
  LIBS cldera-profiling ekat)

# Test Pathway
EkatCreateUnitTest (pathway pathway.cpp
  LIBS cldera-profiling ekat)

# Test ProfilingTestManager
EkatCreateUnitTest (test_manager test_manager.cpp
  LIBS cldera-profiling ekat)

# Set yaml files for tests that utilize input decks
set (CLDERA_TEST_YAMLS
  cldera_field_test_input.yaml
  cldera_graph_input.yaml
  cldera_pathway_input.yaml)

# Copy the yaml files to the destination
FILE(COPY ${CLDERA_TEST_YAMLS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
