
if(32BIT)
remove_definitions(-DOVERLOAD_R4)
remove_definitions(-DOVERLOAD_R8)
message ("Force 64 bits in stochastic_physics")
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    if(REPRO)
        string (REPLACE "-i4 -real-size 32" "-i4 -real-size 64" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
    else()
        string (REPLACE "-i4 -real-size 32" "-i4 -real-size 64 -no-prec-div -no-prec-sqrt" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
    endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
endif()
endif()

add_library(
    stochastic_physics

    ../../stochastic_physics/plumes.f90

    ../../stochastic_physics/stochy_gg_def.f
    ../../stochastic_physics/stochy_resol_def.f
    ../../stochastic_physics/stochy_layout_lag.f
    ../../stochastic_physics/four_to_grid_stochy.F
    ../../stochastic_physics/glats_stochy.f
    ../../stochastic_physics/sumfln_stochy.f
    ../../stochastic_physics/gozrineo_stochy.f
    ../../stochastic_physics/num_parthds_stochy.f
    ../../stochastic_physics/get_ls_node_stochy.f
    ../../stochastic_physics/get_lats_node_a_stochy.f
    ../../stochastic_physics/setlats_a_stochy.f
    ../../stochastic_physics/setlats_lag_stochy.f
    ../../stochastic_physics/epslon_stochy.f
    ../../stochastic_physics/getcon_lag_stochy.f
    ../../stochastic_physics/pln2eo_stochy.f
    ../../stochastic_physics/dozeuv_stochy.f
    ../../stochastic_physics/dezouv_stochy.f

    ../../stochastic_physics/spectral_layout.F90
    ../../stochastic_physics/getcon_spectral.F90
    ../../stochastic_physics/stochy_namelist_def.F90
    ../../stochastic_physics/compns_stochy.F90
    ../../stochastic_physics/stochy_internal_state_mod.F90
    ../../stochastic_physics/stochastic_physics.F90
    ../../stochastic_physics/stochy_patterngenerator.F90
    ../../stochastic_physics/stochy_data_mod.F90
    ../../stochastic_physics/get_stochy_pattern.F90
    ../../stochastic_physics/initialize_spectral_mod.F90
    ../../stochastic_physics/cellular_automata_global.F90
    ../../stochastic_physics/cellular_automata_sgs.F90
    ../../stochastic_physics/update_ca.F90
)

if(CCPP)
  target_include_directories(stochastic_physics PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/framework/src)
  target_include_directories(stochastic_physics PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/physics)
endif()
target_link_libraries(stochastic_physics fms)
target_link_libraries(stochastic_physics gfsphysics)
target_link_libraries(stochastic_physics fv3core)
