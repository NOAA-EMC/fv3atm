
if(32BIT)
remove_definitions(-DOVERLOAD_R4)
remove_definitions(-DOVERLOAD_R8)
message ("Force 64 bits in namphysics")
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
    namphysics

    physics/grrad_nmmb.f
    physics/gscond.f
    physics/iounitdef.f
    physics/mersenne_twister.f
    physics/mfdeepcnv.f
    physics/mfpbl.f
    physics/mfshalcnv.f
    physics/moninedmf.f
    physics/moninq.f
    physics/namelist_soilveg.f
    physics/ozne_def.f
    physics/physparam.f
    physics/precpd_nmmb.f
    physics/rad_initialize_nmmb.f
    physics/radiation_aerosols_nmmb.f
    physics/radiation_astronomy_nmmb.f
    physics/radiation_clouds_nmmb.f
    physics/radiation_gases_nmmb.f
    physics/radiation_surface_nmmb.f
    physics/radlw_datatb.f
    physics/radlw_main_nmmb.f
    physics/radlw_param.f
    physics/radsw_datatb.f
    physics/radsw_main_nmmb.f
    physics/radsw_param.f
    physics/sascnvn.f
    physics/h2o_def.f
    physics/shalcnv.f

    physics/machine.F

    physics/funcphys.f90
    physics/physcons.f90
    physics/noahmp_tables.f90

    physics/module_BL_GFSPBL.F90
    physics/module_BL_GFSPBLEDMFHUR.F90
    physics/module_BL_GFSPBLHUR.F90
    physics/module_BL_MYJPBL.F90
    physics/module_CONSTANTS.F90
    physics/module_CONTROL.F90
    physics/module_CONVECTION.F90
    physics/module_CU_BMJ.F90
    physics/module_CU_SAS.F90
    physics/module_CU_SASHUR.F90
    physics/module_CU_SCALE.F90
    physics/module_DIAGNOSE.F90
    physics/module_GWD.F90
    physics/module_KINDS.F90
    physics/module_LS_LISS.F90
    physics/module_LS_NOAHLSM.F90
    physics/module_MICROPHYSICS.F90
    physics/module_MP_ETANEW.F90
    physics/module_MP_FER_HIRES.F90
    physics/module_MP_GFS.F90
    physics/module_MP_WSM6.F90
    physics/module_RADIATION.F90
    physics/module_RA_GFDL.F90
    physics/module_RA_RRTM.F90
    physics/module_SF_GFDL.F90
    physics/module_SF_JSFC.F90
    physics/module_SF_URBAN.F90
    physics/module_TURBULENCE.F90
    physics/module_mp_radar.F90
    physics/module_mp_thompson.F90
    physics/GFDL_parse_tracers.F90
    physics/gfdl_cloud_microphys.F90

    NAM_layer/module_SOLVER_GRID_COMP.F90
    NAM_layer/NAM_abstraction_layer.F90
    NAM_layer/NAM_diagnostics.F90
    NAM_layer/NAM_driver.F90
    NAM_layer/NAM_radiation_driver.F90
    NAM_layer/NAM_restart.F90
    NAM_layer/NAM_typedefs.F90
)

set_property(SOURCE physics/grrad_nmmb.f               APPEND_STRING PROPERTY COMPILE_FLAGS "-free -r8")
set_property(SOURCE physics/radiation_astronomy_nmmb.f APPEND_STRING PROPERTY COMPILE_FLAGS "-free -r8")
set_property(SOURCE physics/radiation_clouds_nmmb.f    APPEND_STRING PROPERTY COMPILE_FLAGS "-free -r8")
set_property(SOURCE physics/radsw_main_nmmb.f          APPEND_STRING PROPERTY COMPILE_FLAGS "-free -r8")

set_property(SOURCE NAM_layer/NAM_diagnostics.F90      APPEND_STRING PROPERTY COMPILE_FLAGS "-O0")

target_include_directories(namphysics PRIVATE ${ESMF_MOD})

set_target_properties(namphysics PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(namphysics PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/mod)
