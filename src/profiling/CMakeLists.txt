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
    cldera_field_layout.cpp
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
    stats/cldera_register_stats.cpp
    stats/cldera_field_identity.cpp
    stats/cldera_field_stat.cpp
    stats/cldera_field_masked_integral.cpp
    stats/cldera_field_bounded_masked_integral.cpp
    stats/cldera_field_bounded.cpp
    stats/cldera_field_sum_along_columns.cpp
    stats/cldera_field_zonal_mean.cpp
)
set (MODULES_DIR ${CMAKE_CURRENT_BINARY_DIR}/profiling_modules)
set_target_properties(cldera-profiling PROPERTIES
  Fortran_MODULE_DIRECTORY ${MODULES_DIR})

target_link_libraries (cldera-profiling PUBLIC cldera-timing cldera-pnetcdf ekat)

target_include_directories (cldera-profiling
  PUBLIC
  "$<BUILD_INTERFACE:${CLDERA_SOURCE_DIR}/src;${CLDERA_BINARY_DIR}/src;${MODULES_DIR}>"
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set (CLDERA_PROFILING_HEADERS
  cldera_bounds_field_test.hpp
  cldera_data_type.hpp
  cldera_field.hpp
  cldera_field_layout.hpp
  cldera_field_test.hpp
  cldera_field_test_factory.hpp
  cldera_graph.hpp
  cldera_graph_factory.hpp
  cldera_graph_vertex.hpp
  cldera_max_field_test.hpp
  cldera_min_field_test.hpp
  cldera_mpi_timing_wrappers.hpp
  cldera_pathway.hpp
  cldera_pathway_factory.hpp
  cldera_profiling_archive.hpp
  cldera_profiling_interface.hpp
  cldera_profiling_session.hpp
  cldera_profiling_test_manager.hpp
  cldera_profiling_types.hpp
  cldera_time_stamp.hpp
  stats/cldera_field_avg_along_columns.hpp
  stats/cldera_field_bounded.hpp
  stats/cldera_field_bounding_box.hpp
  stats/cldera_field_global_avg.hpp
  stats/cldera_field_global_max.hpp
  stats/cldera_field_global_min.hpp
  stats/cldera_field_global_sum.hpp
  stats/cldera_field_identity.hpp
  stats/cldera_field_bounded_masked_integral.hpp
  stats/cldera_field_masked_integral.hpp
  stats/cldera_field_max_along_columns.hpp
  stats/cldera_field_min_along_columns.hpp
  stats/cldera_field_pnetcdf_reference.hpp
  stats/cldera_field_stat.hpp
  stats/cldera_field_stat_along_axis.hpp
  stats/cldera_field_stat_pipe.hpp
  stats/cldera_field_stat_utils.hpp
  stats/cldera_field_sum_along_columns.hpp
  stats/cldera_field_vertical_contraction.hpp
  stats/cldera_field_zonal_mean.hpp
  stats/cldera_register_stats.hpp
  utils/cldera_subview_utils.hpp
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
