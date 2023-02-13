# Invoke this script as
#  cmake -P CprncTest.cmake src_nc_file tgt_nc_file

if (NOT ${CMAKE_ARGC} EQUAL 5)
  message (FATAL_ERROR "CprncTest should be invoked with 2 arguments (src and tgt nc file).")
endif()

set (SRC_FILE ${CMAKE_ARGV3})
set (TGT_FILE ${CMAKE_ARGV4})

# Dump src file, but skip 1st line, which contains filename
set (CMD ncdump -v T,V,I ${SRC_FILE})
execute_process (
  COMMAND ${CMD}
  COMMAND tail -n +2
  RESULT_VARIABLE ncdump_src_result
  OUTPUT_VARIABLE ncdump_src_output
  ERROR_VARIABLE  ncdump_src_output
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

if (NOT ncdump_src_result EQUAL 0)
  string (CONCAT msg
          "Command\n"
          "  '${CMD}'\n"
          "returned '${ncdump_src_result}', and output:\n"
          "${ncdump_src_output}")
  message ("${msg}")
  message (FATAL_ERROR "Aborting.")
endif()

# Dump tgt file, but skip 1st line, which contains filename
set (CMD ncdump -v T,V,I ${TGT_FILE})
execute_process (
  COMMAND ${CMD}
  COMMAND tail -n +2
  RESULT_VARIABLE ncdump_tgt_result
  OUTPUT_VARIABLE ncdump_tgt_output
  ERROR_VARIABLE  ncdump_tgt_output
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

if (NOT ncdump_tgt_result EQUAL 0)
  string (CONCAT msg
          "Command\n"
          "  '${CMD}'\n"
          "returned '${ncdump_tgt_result}', and output:\n"
          "${ncdump_tgt_output}")
  message ("${msg}")
  message (FATAL_ERROR "Aborting.")
endif()

# Compare dumped outputs
if (NOT "${ncdump_src_output}" STREQUAL "${ncdump_tgt_output}")
  string (CONCAT msg
          "The two nc files appear to be different.\n"
          "  src file: ${SRC_FILE}\n"
          "  tgt file: ${TGT_FILE}\n")
  message ("${msg}")
  message ("${cprnc_output}")
  message (FATAL_ERROR "Aborting.")
endif()
