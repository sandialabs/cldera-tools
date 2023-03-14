include (EkatUtils)

# Generate cldera_config.h
EkatConfigFile (${CLDERA_SOURCE_DIR}/src/cldera_config.h.in
                ${CLDERA_BINARY_DIR}/src/cldera_config.h
                F90_FILE ${CLDERA_BINARY_DIR}/src/cldera_config.f)

# Generate cldera_version.h
execute_process(
  COMMAND git rev-parse HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE CLDERA_GIT_SHA
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/cldera_version.h.in
                ${CMAKE_CURRENT_BINARY_DIR}/cldera_version.h)


# Enable sub-packages

option (CLDERA_ENABLE_PROFILING_TOOL "Enable profiling source code" ON)
set (CLDERA_PNETCDF_PATH "" CACHE PATH "Path to pnetcdf installation")

if (CLDERA_PNETCDF_PATH)
  add_subdirectory(io)
endif()

if (CLDERA_ENABLE_PROFILING_TOOL)
  add_subdirectory(profiling)
endif()

include(GNUInstallDirs)
install (FILES ${CLDERA_BINARY_DIR}/src/cldera_config.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install (FILES ${CLDERA_BINARY_DIR}/src/cldera_version.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
