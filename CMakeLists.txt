# Call to CCPP code generator
if(DEBUG)
  set(_ccpp_debug_arg "--debug")
endif()
if(DEFINED CCPP_SUITES)
  set(_ccpp_suites_arg "--suites=${CCPP_SUITES}")
  message("Calling CCPP code generator (ccpp_prebuild.py) for suites ${_ccpp_suites_arg} ...")
else()
  message("Calling CCPP code generator (ccpp_prebuild.py) for all available suites ...")
endif()
execute_process(COMMAND ${Python_EXECUTABLE}
                        "ccpp/framework/scripts/ccpp_prebuild.py"
                        "--config=ccpp/config/ccpp_prebuild_config.py"
                        "--builddir=${CMAKE_CURRENT_BINARY_DIR}" ${_ccpp_suites_arg} ${_ccpp_debug_arg}
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/ccpp_prebuild.out
                ERROR_FILE  ${CMAKE_CURRENT_BINARY_DIR}/ccpp_prebuild.err
                RESULT_VARIABLE RC)
# Check return code from ccpp_prebuild.py
if(NOT RC EQUAL 0)
  message(FATAL_ERROR "An error occured while running ccpp_prebuild.py, check ${CMAKE_CURRENT_BINARY_DIR}/ccpp_prebuild.{out,err}")
endif()

add_subdirectory(cpl)
add_subdirectory(io)

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
  atmos_cubed_sphere/tools/external_ic.F90
  atmos_cubed_sphere/tools/external_sst.F90
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
  set(DYN32 ON CACHE BOOL "Enable support for 32bit fast physics in CCPP")
  list(APPEND _fv3dycore_defs_private OVERLOAD_R4
                                      OVERLOAD_R8)
else()
  set(DYN32 OFF CACHE BOOL "Disable support for 32bit fast physics in CCPP")
endif()

list(APPEND _fv3dycore_defs_private CCPP)

if(OpenMP_Fortran_FOUND)
  list(APPEND _fv3dycore_defs_private OPENMP)
endif()

set_property(SOURCE atmos_cubed_sphere/model/nh_utils.F90 APPEND_STRING PROPERTY COMPILE_FLAGS "${FAST}")
set_property(SOURCE atmos_cubed_sphere/model/fv_mapz.F90  APPEND_STRING PROPERTY COMPILE_FLAGS "${FAST}")

set_target_properties(fv3dycore PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)

target_compile_definitions(fv3dycore PRIVATE "${_fv3dycore_defs_private}")
target_include_directories(fv3dycore PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/atmos_cubed_sphere
                                             ${CMAKE_CURRENT_BINARY_DIR}/ccpp/driver/mod)
target_include_directories(fv3dycore INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                               $<INSTALL_INTERFACE:mod>)

target_link_libraries(fv3dycore PUBLIC fms
                                       ccppphys
                                       ccppdata
                                       esmf)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(fv3dycore PUBLIC OpenMP::OpenMP_Fortran)
endif()

###############################################################################
### CCPP
###############################################################################

add_subdirectory(ccpp)
add_subdirectory(ccpp/data)
add_subdirectory(ccpp/driver)
add_dependencies(ccppphys   ccpp)
add_dependencies(ccppdata   ccpp ccppphys)
add_dependencies(ccppdriver ccpp ccppphys ccppdata)
add_dependencies(fv3dycore  ccppdriver ccpp ccppphys ccppdata)
target_include_directories(fv3dycore PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/ccpp/framework/src
                                             ${CMAKE_CURRENT_BINARY_DIR}/ccpp/physics
                                             ${CMAKE_CURRENT_BINARY_DIR}/ccpp/driver)
target_link_libraries(ccppphys PUBLIC sp::sp_d
                                      w3nco::w3nco_d)

###############################################################################
### stochastic_physics
###############################################################################
# When specifying an out-of-tree source a binary directory must be explicitly specified.
#add_subdirectory(../stochastic_physics stochastic_physics)
add_subdirectory(stochastic_physics)

###############################################################################
### fv3atm
###############################################################################
add_library(fv3atm
  atmos_model.F90
  fv3_cap.F90
  module_fv3_config.F90
  module_fcst_grid_comp.F90
  time_utils.F90)

list(APPEND _fv3atm_defs_private GFS_PHYS
                                 INTERNAL_FILE_NML
                                 ESMF_VERSION_MAJOR=${ESMF_VERSION_MAJOR})

set_target_properties(fv3atm PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(fv3atm INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                            $<INSTALL_INTERFACE:mod>)

list(APPEND _fv3atm_defs_private CCPP)
target_include_directories(fv3atm PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/ccpp/physics
                                          ${CMAKE_CURRENT_BINARY_DIR}/ccpp/driver/mod)
set(CCPP_LIBRARIES ccppdriver ccppdata ccppphys ccpp)
add_dependencies(fv3atm ccppdriver ccppdata ccppphys ccpp)
target_link_libraries(fv3atm PUBLIC ccppdriver ccppdata ccppphys ccpp)

target_include_directories(fv3atm PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/stochastic_physics)
target_compile_definitions(fv3atm PRIVATE "${_fv3atm_defs_private}")

target_link_libraries(fv3atm PUBLIC fv3dycore
                                    cpl
                                    io
                                    stochastic_physics
                                    stochastic_physics_wrapper)

target_link_libraries(fv3atm PUBLIC nemsio::nemsio
                                    w3emc::w3emc_d
                                    w3nco::w3nco_d
                                    sp::sp_d
                                    bacio::bacio_4
                                    esmf)

if(INLINE_POST)
  target_link_libraries(fv3atm PUBLIC upp::upp)
endif()

if(OpenMP_Fortran_FOUND)
  target_link_libraries(fv3atm PUBLIC OpenMP::OpenMP_Fortran)
endif()

###############################################################################
### Install
###############################################################################
install(
  TARGETS fv3atm fv3dycore io ${CCPP_LIBRARIES} cpl stochastic_physics stochastic_physics_wrapper
  EXPORT fv3atm-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX})

install(EXPORT fv3atm-config
  DESTINATION lib/cmake)
