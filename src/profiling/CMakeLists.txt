include(GNUInstallDirs)

if (NOT TARGET cldera-pnetcdf)
  string (CONCAT msg
    " +----------------------- ERROR -------------------------+\n"
    " | Error! Cldera profiling tool requires cldera-pnetcdf. |\n"
    " | Please, reconfigure with CLDERA_PNETCDF_PATH pointing |\n"
    " | to a valid installation of pnetcdf.                   |\n"
    " +-------------------------------------------------------+\n")
  message ("${msg}")
  message (FATAL_ERROR "Aborting.")
endif()

add_library (cldera-profiling
    cldera_bounds_field_test.cpp
    cldera_data_type.cpp
    cldera_max_field_test.cpp
    cldera_min_field_test.cpp
    cldera_field_test_factory.cpp
    cldera_field.cpp
    cldera_graph_vertex.cpp
    cldera_graph.cpp
    cldera_graph_factory.cpp
    cldera_pathway.cpp
    cldera_pathway_factory.cpp
    cldera_profiling_session.cpp
    cldera_profiling_archive.cpp
    cldera_profiling_interface.cpp
    cldera_profiling_test_manager.cpp
    cldera_interface_mod.F90
    cldera_interface_f2c_mod.F90
    stats/cldera_field_stat_factory.cpp
)
set (MODULES_DIR ${CMAKE_CURRENT_BINARY_DIR}/profiling_modules)
set_target_properties(cldera-profiling PROPERTIES
  Fortran_MODULE_DIRECTORY ${MODULES_DIR})

target_link_libraries (cldera-profiling PUBLIC cldera-pnetcdf ekat)

target_include_directories (cldera-profiling
  PUBLIC
  "$<BUILD_INTERFACE:${CLDERA_SOURCE_DIR}/src;${CLDERA_BINARY_DIR}/src;${MODULES_DIR}>"
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set (CLDERA_PROFILING_HEADERS
  cldera_field.hpp
  cldera_field_layout.hpp
  cldera_profiling_session.hpp
  cldera_profiling_archive.hpp
  cldera_profiling_interface.hpp
  cldera_profiling_test_manager.hpp
  cldera_profiling_types.hpp
  stats/cldera_field_stat.hpp
  stats/cldera_field_global_max.hpp
  stats/cldera_field_global_min.hpp
  stats/cldera_field_global_sum.hpp
  stats/cldera_field_global_avg.hpp
  stats/cldera_field_stat_factory.hpp
)
set_target_properties (cldera-profiling PROPERTIES
  PUBLIC_HEADER "${CLDERA_PROFILING_HEADERS}")

install (TARGETS cldera-profiling
  EXPORT        ClderaTargets
  LIBRARY       DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE       DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install (DIRECTORY ${MODULES_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
