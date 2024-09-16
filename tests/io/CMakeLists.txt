include (EkatCreateUnitTest)

# Test pnetcdf IO
EkatCreateUnitTest (pnetcdf-tests cldera_pnetcdf_tests.cpp
  LIBS cldera-pnetcdf
  MPI_RANKS 1 ${CLDERA_TESTS_MAX_RANKS}
  FIXTURES_SETUP setup_nc_files)

set (compare_script ${CMAKE_CURRENT_SOURCE_DIR}/compare_nc.cmake)

# Check monolithic vs restarted output
foreach (RANK RANGE 1 ${CLDERA_TESTS_MAX_RANKS})
  set (SRC_FILE ./monolithic_np${RANK}.nc)
  set (TGT_FILE ./restarted_np${RANK}.nc)
  add_test (NAME cldera-pnetcdf-restart-check-np${RANK}
            COMMAND ${CMAKE_COMMAND} -P ${compare_script} ${SRC_FILE} ${TGT_FILE}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties (cldera-pnetcdf-restart-check-np${RANK} PROPERTIES
            FIXTURES_REQUIRED setup_nc_files)
endforeach()

# Compare simple output across number of ranks
if (CLDERA_TESTS_MAX_RANKS GREATER 1)
  foreach (RANK RANGE 2 ${CLDERA_TESTS_MAX_RANKS})
    set (SRC_FILE ./simple_test_np${RANK}.nc)
    set (TGT_FILE ./simple_test_np1.nc)
    add_test (NAME cldera-pnetcdf-np1-vs-np${RANK}
              COMMAND ${CMAKE_COMMAND} -P ${compare_script} ${SRC_FILE} ${TGT_FILE}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties (cldera-pnetcdf-np1-vs-np${RANK} PROPERTIES
              FIXTURES_REQUIRED setup_nc_files)
  endforeach()
endif()
