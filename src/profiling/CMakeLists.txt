include(GNUInstallDirs)

add_library (cldera-profiling
    cldera_field.cpp
    cldera_profiling_session.cpp
    cldera_profiling_archive.cpp
    cldera_profiling_interface.cpp
    stats/cldera_compute_max.cpp
    stats/cldera_compute_min.cpp
    stats/cldera_compute_avg.cpp
    stats/cldera_compute_sum.cpp
)

target_link_libraries (cldera-profiling PUBLIC ekat)

target_include_directories (cldera-profiling
  PUBLIC
  "$<BUILD_INTERFACE:${CLDERA_SOURCE_DIR}/src;${CLDERA_BINARY_DIR}/src>"
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set (CLDERA_PROFILING_HEADERS
  cldera_profiling_session.hpp
  cldera_profiling_archive.hpp
  cldera_profiling_interface.hpp
  cldera_profiling_types.hpp
  stats/cldera_compute_stat.hpp
)

install (TARGETS cldera-profiling
  EXPORT        ClderaTargets
  LIBRARY       DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE       DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (FILES ${CLDERA_PROFILING_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
