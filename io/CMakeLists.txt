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

set_source_files_properties(post_gfs.F90 post_nems_routines.F90
                            PROPERTIES
                            COMPILE_FLAGS "-I${POST_INC}")

target_link_libraries(io fms)
target_link_libraries(io gfsphysics)
target_link_libraries(io ipd)

if(CCPP)
  target_include_directories(io PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/framework/src)
  target_include_directories(io PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/physics)
endif()
target_include_directories(io PRIVATE ${ESMF_MOD})
target_include_directories(io PRIVATE ${NEMSIO_INC})

set_target_properties(io PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(io PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/mod)
