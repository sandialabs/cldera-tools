find_package (PNETCDF
  NAMES pnetcdf
  PATHS ${CLDERA_PNETCDF_PATH}
  PATH_SUFFIXES lib lib64)

if (NOT PNETCDF_FOUND)
  # No cmake config file found. Try to locate pnetcdf manually
  find_library (PNETCDF_LIBRARIES
    NAMES pnetcdf REQUIRED
    PATHS ${CLDERA_PNETCDF_PATH}
    PATH_SUFFIXES lib lib64)
  find_path (PNETCDF_INCLUDE_DIR
    NAMES pnetcdf.h REQUIRED
    PATHS ${CLDERA_PNETCDF_PATH}
    PATH_SUFFIXES include)
endif()

add_library (cldera-pnetcdf
  cldera_pnetcdf.cpp)

target_link_libraries (cldera-pnetcdf PUBLIC ekat)
target_include_directories (cldera-pnetcdf
  PUBLIC
  "$<BUILD_INTERFACE:${CLDERA_SOURCE_DIR}/src;${CLDERA_BINARY_DIR}/src>"
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Pnetcdf is a private implementation detail, never exposed via headers.
target_link_libraries (cldera-pnetcdf PRIVATE ${PNETCDF_LIBRARIES})
target_include_directories (cldera-pnetcdf PRIVATE ${PNETCDF_INCLUDE_DIR})
