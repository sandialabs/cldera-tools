include(GNUInstallDirs)

add_library (cldera-timing
    cldera_timing_session.cpp
)
target_link_libraries (cldera-timing PUBLIC ekat)

target_include_directories (cldera-timing
  PUBLIC
  "$<BUILD_INTERFACE:${CLDERA_SOURCE_DIR}/src>"
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set (CLDERA_TIMING_HEADERS
    cldera_timer_history.cpp
    cldera_timing_session.hpp
)
set_target_properties (cldera-timing PROPERTIES
  PUBLIC_HEADER "${CLDERA_PROFILING_HEADERS}")

install (TARGETS cldera-timing
  EXPORT        ClderaTargets
  LIBRARY       DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE       DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
