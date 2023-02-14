include (EkatCreateUnitTest)

# Test pnetcdf IO
EkatCreateUnitTest (pnetcdf-tests cldera_pnetcdf_tests.cpp
  LIBS cldera-pnetcdf
  MPI_RANKS 1 ${CLDERA_TESTS_MAX_RANKS}
  MPI_EXEC_NAME ${CLDERA_TESTS_MPI_EXEC_NAME}
  MPI_NP_FLAG ${CLDERA_TESTS_MPI_NP_FLAG}
  PROPERTIES FIXTURES_SETUP setup_nc_files)

if (CLDERA_TESTS_MAX_RANKS GREATER 1)
  set (scriptName ${CMAKE_CURRENT_SOURCE_DIR}/compare_nc.cmake)
  foreach (RANK RANGE 2 ${CLDERA_TESTS_MAX_RANKS})
    set (SRC_FILE ./test_np${RANK}.nc)
    set (TGT_FILE ./test_np1.nc)
    add_test (NAME cldera-pnetcdf-np1-vs-np${RANK}
              COMMAND ${CMAKE_COMMAND} -P ${scriptName} ${SRC_FILE} ${TGT_FILE}
              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties (cldera-pnetcdf-np1-vs-np${RANK} PROPERTIES
              FIXTURES_REQUIRED setup_nc_files)
  endforeach()
endif()
