
###############################################################################
### CCPP
###############################################################################

add_subdirectory(ccpp)

###############################################################################
### fv3dycore
###############################################################################
list(APPEND _fv3dycore_srcs
  atmos_cubed_sphere/model/a2b_edge.F90
  atmos_cubed_sphere/model/multi_gases.F90
  atmos_cubed_sphere/model/boundary.F90
  atmos_cubed_sphere/model/dyn_core.F90
  atmos_cubed_sphere/model/fv_arrays.F90
  atmos_cubed_sphere/model/fv_control.F90
  atmos_cubed_sphere/model/fv_dynamics.F90
  atmos_cubed_sphere/model/fv_fill.F90
  atmos_cubed_sphere/model/fv_grid_utils.F90
  atmos_cubed_sphere/model/fv_mapz.F90
  atmos_cubed_sphere/model/fv_nesting.F90
  atmos_cubed_sphere/model/fv_regional_bc.F90
  atmos_cubed_sphere/model/fv_sg.F90
  atmos_cubed_sphere/model/fv_tracer2d.F90
  atmos_cubed_sphere/model/fv_update_phys.F90
  atmos_cubed_sphere/model/sw_core.F90
  atmos_cubed_sphere/model/tp_core.F90
  atmos_cubed_sphere/model/nh_core.F90
  atmos_cubed_sphere/model/nh_utils.F90
  atmos_cubed_sphere/tools/coarse_grained_diagnostics.F90
  atmos_cubed_sphere/tools/coarse_grained_restart_files.F90
  atmos_cubed_sphere/tools/coarse_graining.F90
  atmos_cubed_sphere/tools/external_ic.F90
  atmos_cubed_sphere/tools/external_sst.F90
  atmos_cubed_sphere/tools/fv_diag_column.F90
  atmos_cubed_sphere/tools/fv_diagnostics.F90
  atmos_cubed_sphere/tools/fv_eta.F90
  atmos_cubed_sphere/tools/fv_grid_tools.F90
  atmos_cubed_sphere/tools/fv_io.F90
  atmos_cubed_sphere/tools/fv_mp_mod.F90
  atmos_cubed_sphere/tools/fv_nudge.F90
  atmos_cubed_sphere/tools/fv_treat_da_inc.F90
  atmos_cubed_sphere/tools/fv_iau_mod.F90
  atmos_cubed_sphere/tools/fv_restart.F90
  atmos_cubed_sphere/tools/fv_surf_map.F90
  atmos_cubed_sphere/tools/fv_timing.F90
  atmos_cubed_sphere//tools/init_hydro.F90
  atmos_cubed_sphere/tools/sim_nc_mod.F90
  atmos_cubed_sphere/tools/sorted_index.F90
  atmos_cubed_sphere/tools/test_cases.F90
  atmos_cubed_sphere/driver/fvGFS/DYCORE_typedefs.F90
  atmos_cubed_sphere/driver/fvGFS/fv_nggps_diag.F90
  atmos_cubed_sphere/driver/fvGFS/atmosphere.F90)

add_library(fv3dycore ${_fv3dycore_srcs})

list(APPEND _fv3dycore_defs_private SPMD
                                    use_WRTCOMP
                                    GFS_PHYS
                                    GFS_TYPES
                                    USE_GFSL63
                                    MOIST_CAPPA
                                    USE_COND)

if(MULTI_GASES)
  list(APPEND _fv3dycore_defs_private MULTI_GASES)
endif()

if(32BIT)
  list(APPEND _fv3dycore_defs_private OVERLOAD_R4
                                      OVERLOAD_R8)
endif()

if(OPENMP)
  list(APPEND _fv3dycore_defs_private OPENMP)
endif()

set_property(SOURCE atmos_cubed_sphere/model/nh_utils.F90 APPEND_STRING PROPERTY COMPILE_FLAGS "${FAST}")
set_property(SOURCE atmos_cubed_sphere/model/fv_mapz.F90  APPEND_STRING PROPERTY COMPILE_FLAGS "${FAST}")

set_target_properties(fv3dycore PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(fv3dycore INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                               $<INSTALL_INTERFACE:mod>)

target_compile_definitions(fv3dycore PRIVATE "${_fv3dycore_defs_private}")

#  So much for being consistent:
# atmos_cubed_sphere/tools/fv_diagnostics.F90: #include <fv_diagnostics.h>
# atmos_cubed_sphere/tools/fv_eta.F90:         #include <tools/fv_eta.h>
target_include_directories(fv3dycore PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/atmos_cubed_sphere
                                             ${CMAKE_CURRENT_SOURCE_DIR}/atmos_cubed_sphere/tools)

target_link_libraries(fv3dycore PUBLIC fms
                                       fv3ccpp
                                       esmf)
if(OPENMP)
  target_link_libraries(fv3dycore PUBLIC OpenMP::OpenMP_Fortran)
endif()

# This should not be necessary once framework and physics targets define BUILD_INTERFACE
target_include_directories(fv3dycore PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/ccpp/framework/src
                                             ${CMAKE_CURRENT_BINARY_DIR}/ccpp/physics)

###############################################################################
### fv3atm
###############################################################################
if(INLINE_POST)
  set(POST_SRC io/inline_post.F90 io/post_nems_routines.F90 io/post_gfs.F90 io/post_regional.F90)
else()
  set(POST_SRC io/inline_post_stub.F90)
  list(APPEND _fv3atm_defs_private NO_INLINE_POST)
endif()

if(NOT PARALLEL_NETCDF)
  list(APPEND _fv3atm_defs_private NO_PARALLEL_NETCDF)
endif()

add_library(fv3atm
  atmos_model.F90
  fv3_cap.F90
  module_fv3_config.F90
  module_fcst_grid_comp.F90
  stochastic_physics/stochastic_physics_wrapper.F90
  cpl/module_block_data.F90
  cpl/module_cplfields.F90
  cpl/module_cap_cpl.F90
  io/ffsync.F90
  io/FV3GFS_io.F90
  io/module_write_nemsio.F90
  io/module_write_netcdf.F90
  io/module_write_netcdf_parallel.F90
  io/module_fv3_io_def.F90
  io/module_write_internal_state.F90
  io/module_wrt_grid_comp.F90
  ${POST_SRC}
)

list(APPEND _fv3atm_defs_private GFS_PHYS
                                 INTERNAL_FILE_NML
                                 use_WRTCOMP
                                 ESMF_VERSION_MAJOR=${ESMF_VERSION_MAJOR})

target_compile_definitions(fv3atm PRIVATE "${_fv3atm_defs_private}")

set_target_properties(fv3atm PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(fv3atm INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                            $<INSTALL_INTERFACE:mod>)

# This should not be necessary once framework and physics targets define BUILD_INTERFACE
target_include_directories(fv3atm PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/ccpp/physics)

target_link_libraries(fv3atm PUBLIC fv3dycore
                                    fv3ccpp
                                    stochastic_physics
                                    fms)

target_link_libraries(fv3atm PUBLIC nemsio::nemsio
                                    w3emc::w3emc_d
                                    w3nco::w3nco_d
                                    sp::sp_d
                                    bacio::bacio_4
                                    esmf)

if(INLINE_POST)
  target_link_libraries(fv3atm PUBLIC upp::upp)
endif()

if(OPENMP)
  target_link_libraries(fv3atm PUBLIC OpenMP::OpenMP_Fortran)
endif()

###############################################################################
### Install
###############################################################################
install(
  TARGETS fv3atm fv3dycore fv3ccpp
  EXPORT fv3atm-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX})

install(EXPORT fv3atm-config
  DESTINATION lib/cmake)
