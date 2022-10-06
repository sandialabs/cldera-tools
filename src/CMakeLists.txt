include (EkatUtils)

# Generate cldera_config.h
EkatConfigFile (${CLDERA_SOURCE_DIR}/src/cldera_config.h.in
                ${CLDERA_BINARY_DIR}/src/cldera_config.h
                F90_FILE ${CLDERA_BINARY_DIR}/src/cldera_config.f)

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
