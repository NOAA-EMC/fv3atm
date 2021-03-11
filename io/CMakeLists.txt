list(APPEND _io_defs_private use_WRTCOMP)

if(INLINE_POST)
  set(POST_SRC inline_post.F90 post_nems_routines.F90 post_gfs.F90 post_regional.F90)
else()
  set(POST_SRC inline_post_stub.F90)
  list(APPEND _io_defs_private NO_INLINE_POST)
endif()

if(NOT PARALLEL_NETCDF)
  list(APPEND _io_defs_private NO_PARALLEL_NETCDF)
endif()

list(APPEND _io_defs_private CCPP)

add_library(
    io

    ffsync.F90
    FV3GFS_io.F90
    module_write_nemsio.F90
    module_write_netcdf.F90
    module_write_netcdf_parallel.F90
    module_fv3_io_def.F90
    module_write_internal_state.F90
    module_wrt_grid_comp.F90
    ${POST_SRC}
)

set_target_properties(io PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_compile_definitions(io PRIVATE "${_io_defs_private}")
target_include_directories(io PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)

target_link_libraries(io PRIVATE fms
                                 ccppdriver)

if(INLINE_POST)
  target_link_libraries(io PRIVATE upp::upp)
endif()

target_include_directories(io PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/framework/src
                                      ${CMAKE_BINARY_DIR}/FV3/ccpp/physics
                                      ${CMAKE_BINARY_DIR}/FV3/ccpp/driver/mod)

target_link_libraries(io PRIVATE nemsio::nemsio
                                 esmf)
