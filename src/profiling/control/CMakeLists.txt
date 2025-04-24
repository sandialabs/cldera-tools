include(GNUInstallDirs)

add_library (cldera-control
  cldera_control.cpp  
  cldera_register_controls.cpp  
)
set (MODULES_DIR ${CMAKE_CURRENT_BINARY_DIR}/control_modules)
set_target_properties(cldera-control PROPERTIES
  Fortran_MODULE_DIRECTORY ${MODULES_DIR})

target_link_libraries (cldera-control PUBLIC cldera-profiling cldera-timing cldera-pnetcdf ekat)

target_include_directories (cldera-control
  PUBLIC
  "$<BUILD_INTERFACE:${CLDERA_SOURCE_DIR}/src;${CLDERA_BINARY_DIR}/src;${MODULES_DIR}>"
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set (CLDERA_CONTROL_HEADERS
  cldera_control.hpp
  cldera_register_controls.hpp
)
set_target_properties (cldera-control PROPERTIES
  PUBLIC_HEADER "${CLDERA_CONTROL_HEADERS}")

install (TARGETS cldera-control
  EXPORT        ClderaTargets
  LIBRARY       DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE       DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install (DIRECTORY ${MODULES_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
