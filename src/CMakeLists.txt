# Generate cldera_config.h
configure_file (${CLDERA_SOURCE_DIR}/src/cldera_config.h.in
                ${CLDERA_BINARY_DIR}/src/cldera_config.h)

if (CLDERA_ENABLE_PROFILING_TOOL)
  add_subdirectory(profiling)
endif()

include(GNUInstallDirs)
install (FILES ${CLDERA_BINARY_DIR}/src/cldera_config.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
