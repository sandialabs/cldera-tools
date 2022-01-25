add_library (cldera-profiling)

target_sources (cldera-profiling
  PUBLIC
    cldera_profiling_session.cpp
)

target_link_libraries (cldera-profiling PUBLIC ekat)
target_include_directories (cldera-profiling
  PUBLIC
    ${CLDERA_SOURCE_DIR}/src
    ${CLDERA_BINARY_DIR}/src) # For cldera_config.h
