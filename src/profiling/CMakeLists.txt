include(GNUInstallDirs)

add_library (cldera-profiling
    cldera_bounds_field_test.cpp
    cldera_field.cpp
    cldera_graph_vertex.cpp
    cldera_graph.cpp
    cldera_graph_factory.cpp
    cldera_profiling_session.cpp
    cldera_profiling_archive.cpp
    cldera_profiling_interface.cpp
    cldera_profiling_test_manager.cpp
    cldera_interface_mod.F90
    cldera_interface_f2c_mod.F90
    stats/cldera_compute_max.cpp
    stats/cldera_compute_min.cpp
    stats/cldera_compute_avg.cpp
    stats/cldera_compute_sum.cpp
)
set (MODULES_DIR ${CMAKE_CURRENT_BINARY_DIR}/profiling_modules)
set_target_properties(cldera-profiling PROPERTIES
  Fortran_MODULE_DIRECTORY ${MODULES_DIR})

target_link_libraries (cldera-profiling PUBLIC ekat)

target_include_directories (cldera-profiling
  PUBLIC
  "$<BUILD_INTERFACE:${CLDERA_SOURCE_DIR}/src;${CLDERA_BINARY_DIR}/src;${MODULES_DIR}>"
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set (CLDERA_PROFILING_HEADERS
  cldera_profiling_session.hpp
  cldera_profiling_archive.hpp
  cldera_profiling_interface.hpp
  cldera_profiling_test_manager.hpp
  cldera_profiling_types.hpp
  stats/cldera_compute_stat.hpp
)

install (TARGETS cldera-profiling
  EXPORT        ClderaTargets
  LIBRARY       DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE       DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (FILES ${CLDERA_PROFILING_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install (DIRECTORY ${MODULES_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
