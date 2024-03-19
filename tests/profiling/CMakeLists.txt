include (EkatCreateUnitTest)

# Test ProfilingContext
EkatCreateUnitTest (context context.cpp
  LIBS cldera-profiling ekat
  MPI_RANKS 1 ${CLDERA_TESTS_MAX_RANKS}
)

# Test Field
EkatCreateUnitTest (time_stamp time_stamp.cpp
  LIBS cldera-profiling ekat)

# Test Field
EkatCreateUnitTest (field field.cpp
  LIBS cldera-profiling ekat)

# Test ProfilingArchive
EkatCreateUnitTest (archive archive.cpp
  LIBS cldera-profiling ekat
  FIXTURES_SETUP archive_output)

add_test (NAME archive_check
          COMMAND ncdump -v foo_global_max archive_tests-2022-09-15-43000.nc)
set_tests_properties(archive_check PROPERTIES
  PASS_REGULAR_EXPRESSION "foo_global_max = 3"
  FIXTURES_REQUIRED archive_output
)

# Test subview utils
EkatCreateUnitTest (subview_utils subview_utils.cpp
  LIBS cldera-profiling ekat)

# Test Graph
EkatCreateUnitTest (graph graph.cpp
  LIBS cldera-profiling ekat)

# Test functions to compute stats
EkatCreateUnitTest (compute_stats compute_stats.cpp
  LIBS cldera-profiling ekat)

if (CLDERA_PNETCDF_PATH)
  EkatCreateUnitTest (compute_stats_pnetcdf compute_stats_pnetcdf.cpp
    LIBS cldera-profiling cldera-pnetcdf ekat)
endif()

EkatCreateUnitTest (bounds_field_test bounds_field_test.cpp
  LIBS cldera-profiling ekat)

EkatCreateUnitTest (vertical_contraction vertical_contraction.cpp
  LIBS cldera-profiling ekat)

EkatCreateUnitTest (piped_stats piped_stats.cpp
  LIBS cldera-profiling ekat)

EkatCreateUnitTest (masked_integral masked_integral.cpp
  LIBS cldera-profiling ekat
  MPI_RANKS 1 ${CLDERA_TESTS_MAX_RANKS}
)

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
