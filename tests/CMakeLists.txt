add_subdirectory (timing)

if (CLDERA_ENABLE_PROFILING_TOOL)
  add_subdirectory(profiling)
endif()

if (CLDERA_PNETCDF_PATH)
  add_subdirectory(io)
endif()
