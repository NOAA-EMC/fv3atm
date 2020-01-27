if(INLINE_POST)
  message("inline post is set")
  set(POST_SRC post_gfs.F90 post_nems_routines.F90)
else()
  message("inline post is NOT set")
  set(POST_SRC post_gfs_stub.F90)
  add_definitions(-DNO_INLINE_POST)
endif()

add_library(
    io

    ffsync.F90
    FV3GFS_io.F90
    module_write_nemsio.F90
    module_write_netcdf.F90
    module_fv3_io_def.F90
    module_write_internal_state.F90
    module_wrt_grid_comp.F90
    ${POST_SRC}
)

message("HEY post_inc is ${POST_INC}, and src is ${POST_SRC}")
set_source_files_properties(post_gfs.F90 post_nems_routines.F90
                            PROPERTIES
                            COMPILE_FLAGS "-I${POST_INC}")

target_link_libraries(io fms)
target_link_libraries(io ${PHYS}physics)
target_link_libraries(io ipd)

target_include_directories(io PRIVATE ${ESMF_F90COMPILEPATHS})
target_include_directories(io PRIVATE ${NEMSIO_INC})
target_include_directories(io PRIVATE ${MPI_Fortran_INCLUDE_PATH})
target_include_directories(io PRIVATE ${NCEPLIBS_INSTALL_DIR}/include)
target_include_directories(io PRIVATE ${NCEPLIBS_INSTALL_DIR}/mod)

set_target_properties(io PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(io PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(io PUBLIC ${MPI_Fortran_INCLUDE_PATH})
