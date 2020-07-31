if(CCPP)
    add_definitions(-DCCPP)
endif()

add_subdirectory(cpl)
add_subdirectory(gfsphysics)
add_subdirectory(ipd)
add_subdirectory(io)

###############################################################################
### fv3dycore
###############################################################################
if(NOT CCPP)
    set(FAST_PHYSICS_SRCS_F90 atmos_cubed_sphere/model/fv_cmp.F90)
endif()

add_library(
    fv3dycore

    atmos_cubed_sphere/model/a2b_edge.F90
    atmos_cubed_sphere/model/multi_gases.F90
    atmos_cubed_sphere/model/boundary.F90
    atmos_cubed_sphere/model/dyn_core.F90
    atmos_cubed_sphere/model/fv_arrays.F90
    ${FAST_PHYSICS_SRCS_F90}
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
    atmos_cubed_sphere/driver/fvGFS/atmosphere.F90
)

set_property(SOURCE atmos_cubed_sphere/model/nh_utils.F90 APPEND_STRING PROPERTY COMPILE_FLAGS "${FAST}")
set_property(SOURCE atmos_cubed_sphere/model/fv_mapz.F90  APPEND_STRING PROPERTY COMPILE_FLAGS "${FAST}")

set_target_properties(fv3dycore PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(fv3dycore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)

target_link_libraries(fv3dycore PRIVATE fms)
target_link_libraries(fv3dycore PRIVATE gfsphysics)
target_link_libraries(fv3dycore PRIVATE ipd)
target_link_libraries(fv3dycore PRIVATE esmf)

target_include_directories(fv3dycore PRIVATE atmos_cubed_sphere ${ESMF_MOD})

###############################################################################
### ccpp
###############################################################################
if(CCPP)
    add_subdirectory(ccpp)
    add_subdirectory(ccpp/driver)

    add_dependencies(gfsphysics ccpp ccppphys)
    add_dependencies(ccppdriver ccpp ccppphys)
    add_dependencies(ccppphys ccpp)

    target_include_directories(fv3dycore PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/framework/src)
    target_include_directories(fv3dycore PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/driver)

    target_link_libraries(ccppphys PRIVATE sp::sp_d)
    target_link_libraries(ccppphys PRIVATE w3nco::w3nco_d)
endif()

###############################################################################
### stochastic_physics
###############################################################################
# When specifying an out-of-tree source a binary directory must be explicitly specified.
#add_subdirectory(../stochastic_physics stochastic_physics)
add_subdirectory(stochastic_physics)

###############################################################################
### fv3atm
###############################################################################
add_library(
    fv3atm

    atmos_model.F90
    fv3_cap.F90
    module_fv3_config.F90
    module_fcst_grid_comp.F90
    time_utils.F90
)

set_target_properties(fv3atm PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(fv3atm PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)

if(CCPP)
    target_include_directories(fv3atm PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/ccpp/driver/mod)
    set(CCPP_LIBRARIES ccppdriver ccppphys ccpp)
    add_dependencies(fv3atm ccppdriver ccppphys ccpp)
    target_link_libraries(fv3atm PUBLIC ccppdriver ccppphys ccpp)
endif()
target_include_directories(fv3atm PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/stochastic_physics)
target_compile_definitions(fv3atm PRIVATE -DESMF_VERSION_MAJOR=${ESMF_VERSION_MAJOR})

target_link_libraries(fv3atm PUBLIC fms)
target_link_libraries(fv3atm PUBLIC cpl)
target_link_libraries(fv3atm PUBLIC gfsphysics)
target_link_libraries(fv3atm PUBLIC ipd)
target_link_libraries(fv3atm PUBLIC io)
target_link_libraries(fv3atm PUBLIC fv3dycore)
target_link_libraries(fv3atm PUBLIC stochastic_physics)
target_link_libraries(fv3atm PUBLIC stochastic_physics_wrapper)

target_link_libraries(fv3atm PUBLIC
                      nemsio::nemsio w3emc::w3emc_d w3nco::w3nco_d sp::sp_d bacio::bacio_4)
if(INLINE_POST)
  target_link_libraries(fv3atm PUBLIC nceppost::nceppost)
endif()

target_link_libraries(fv3atm PUBLIC esmf)

###############################################################################
### Install
###############################################################################
install(
  TARGETS fv3atm fv3dycore io ipd gfsphysics ${CCPP_LIBRARIES} cpl stochastic_physics stochastic_physics_wrapper
  EXPORT fv3atm-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX})

install(EXPORT fv3atm-config
  DESTINATION lib/cmake
)

