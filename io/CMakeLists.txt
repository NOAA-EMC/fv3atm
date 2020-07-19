if(INLINE_POST)
  set(POST_SRC post_gfs.F90 post_nems_routines.F90)
else()
  set(POST_SRC post_gfs_stub.F90)
  add_definitions(-DNO_INLINE_POST)
endif()

if(NOT PARALLEL_NETCDF)
  add_definitions(-DNO_PARALLEL_NETCDF)
endif()

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

target_link_libraries(io PRIVATE fms)
target_link_libraries(io PRIVATE gfsphysics)
target_link_libraries(io PRIVATE ipd)

if(INLINE_POST)
  target_link_libraries(io PRIVATE nceppost::nceppost)
endif()

if(CCPP)
  target_include_directories(io PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/framework/src)
  target_include_directories(io PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/physics)
endif()

target_link_libraries(io PRIVATE nemsio::nemsio)
target_link_libraries(io PRIVATE esmf)

set_target_properties(io PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(io PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
