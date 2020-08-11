
if(32BIT)
remove_definitions(-DOVERLOAD_R4)
remove_definitions(-DOVERLOAD_R8)
message ("Force 64 bits in gfsphysics")
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

set(CCPP_SOURCES
    physics/mersenne_twister.f
    physics/namelist_soilveg.f
    physics/physparam.f
    physics/set_soilveg.f

    physics/noahmp_tables.f90

    physics/GFDL_parse_tracers.F90
    physics/physcons.F90

    CCPP_layer/CCPP_data.F90
    ${CMAKE_BINARY_DIR}/FV3/ccpp/physics/ccpp_static_api.F90

    GFS_layer/GFS_abstraction_layer.F90
    GFS_layer/GFS_diagnostics.F90
    GFS_layer/GFS_driver.F90
    GFS_layer/GFS_restart.F90
)

set(IPD_SOURCES

    physics/cnvc90.f
    physics/co2hc.f
    physics/date_def.f
    physics/dcyc2.f
    physics/dcyc2.pre.rad.f
    physics/efield.f
    physics/get_prs.f
    physics/gocart_tracer_config_stub.f
    physics/gscond.f
    physics/gscondp.f
    physics/gwdc.f
    physics/gwdps.f
    physics/ugwp_driver_v0.f
    physics/cires_orowam2017.f
    physics/h2o_def.f
    physics/h2oc.f
    physics/h2ohdc.f
    physics/h2ophys.f
    physics/ideaca.f
    physics/idea_co2.f
    physics/idea_composition.f
    physics/idea_dissipation.f
    physics/idea_h2o.f
    physics/idea_ion.f
    physics/idea_o2_o3.f
    physics/idea_phys.f
    physics/idea_solar_heating.f
    physics/idea_tracer.f
    physics/iounitdef.f
    physics/lrgsclr.f
    physics/mersenne_twister.f
    physics/mfpbl.f
    physics/mfpblt.f
    physics/mfpbltq.f
    physics/mfscu.f
    physics/mfscuq.f
    physics/module_bfmicrophysics.f
    physics/moninedmf.f
    physics/moninedmf_hafs.f
    physics/moninp.f
    physics/moninp1.f
    physics/moninq.f
    physics/moninq1.f
    physics/moninshoc.f
    physics/mstadb.f
    physics/mstadbtn.f
    physics/mstadbtn2.f
    physics/mstcnv.f
    physics/namelist_soilveg.f
    physics/ozne_def.f
    physics/iccn_def.f
    physics/aerclm_def.f
    physics/ozphys.f
    physics/ozphys_2015.f
    physics/physparam.f
    physics/precpd.f
    physics/precpd_shoc.f
    physics/precpdp.f
    physics/precpd_shoc.f
    physics/progt2.f
    physics/progtm_module.f
    physics/rad_initialize.f
    physics/radiation_aerosols.f
    physics/radiation_astronomy.f
    physics/radiation_clouds.f
    physics/radiation_gases.f
    physics/radiation_surface.f
    physics/radlw_datatb.f
    physics/radlw_main.f
    physics/radlw_param.f
    physics/radsw_datatb.f
    physics/radsw_main.f
    physics/radsw_param.f
    physics/rascnvv2.f
    physics/rayleigh_damp.f
    physics/rayleigh_damp_mesopause.f
    physics/samfaerosols.f
    physics/samfdeepcnv.f
    physics/samfshalcnv.f
    physics/sascnv.f
    physics/sascnvn.f
    physics/satmedmfvdif.f
    physics/satmedmfvdifq.f
    physics/set_soilveg.f
    physics/sfc_cice.f
    physics/sfc_diag.f
    physics/sfc_diff.f
    physics/sfc_drv.f
    physics/sfc_noahmp_drv.f
    physics/sfc_nst.f
    physics/sfc_ocean.f
    physics/sfc_sice.f
    physics/sflx.f
    physics/shalcnv.f
    physics/shalcv.f
    physics/shalcv_1lyr.f
    physics/shalcv_fixdp.f
    physics/shalcv_opr.f
    physics/tracer_const_h.f
    physics/tridi2t3.f

    physics/calpreciptype.f90
    physics/funcphys.f90
    physics/gcm_shoc.f90
    physics/get_prs_fv3.f90
    physics/h2ointerp.f90
    physics/module_nst_model.f90
    physics/module_nst_parameters.f90
    physics/module_nst_water_prop.f90
    physics/ozinterp.f90
    physics/module_wrf_utl.f90
    physics/noahmp_tables.f90
    physics/module_sf_noahmplsm.f90
    physics/module_sf_noahmp_glacier.f90
    physics/iccninterp.f90
    physics/aerinterp.f90
    physics/wam_f107_kp_mod.f90

    physics/aer_cloud.F
    physics/cldmacro.F
    physics/cldwat2m_micro.F
    physics/gfs_phy_tracer_config.F
    physics/machine.F
    physics/num_parthds.F
    physics/sfcsub.F
    physics/wv_saturation.F

    physics/GFDL_parse_tracers.F90
    physics/gcycle.F90
    physics/cires_ugwp_initialize.F90
    physics/cires_ugwp_module.F90
    physics/cires_ugwp_utils.F90
    physics/cires_ugwp_triggers.F90
    physics/cires_ugwp_solvers.F90
    physics/cires_vert_lsatdis.F90
    physics/cires_vert_orodis.F90
    physics/cires_vert_wmsdis.F90
    physics/gfdl_cloud_microphys.F90
    physics/micro_mg_utils.F90
    physics/micro_mg2_0.F90
    physics/micro_mg3_0.F90
    physics/m_micro_driver.F90
    physics/cs_conv.F90
    physics/module_mp_radar.F90
    physics/module_mp_thompson_gfs.F90
    physics/module_mp_wsm6_fv3.F90
    physics/physcons.F90
    physics/surface_perturbation.F90

    GFS_layer/GFS_abstraction_layer.F90
    GFS_layer/GFS_diagnostics.F90
    GFS_layer/GFS_driver.F90
    GFS_layer/GFS_physics_driver.F90
    GFS_layer/GFS_radiation_driver.F90
    GFS_layer/GFS_restart.F90
    GFS_layer/GFS_typedefs.F90
)

if(CCPP)
  add_definitions(-DCCPP)
  if(DYN32)
    add_definitions(-DOVERLOAD_R4)
  endif()
  add_library(gfsphysics ${CCPP_SOURCES})
  target_include_directories(gfsphysics PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/framework/src)
  target_include_directories(gfsphysics PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/physics)
  target_link_libraries(gfsphysics PRIVATE ccppphys ccpp)
else()
  add_library(gfsphysics ${IPD_SOURCES})
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  if(CMAKE_Platform MATCHES "jet")
    set_property(SOURCE physics/radiation_aerosols.f  APPEND_STRING PROPERTY COMPILE_FLAGS "-axSSE4.2,AVX,CORE-AVX-I")
  else()
    set_property(SOURCE physics/radiation_aerosols.f  APPEND_STRING PROPERTY COMPILE_FLAGS "-xCORE-AVX-I")
  endif()
endif()
set_property(SOURCE GFS_layer/GFS_diagnostics.F90 APPEND_STRING PROPERTY COMPILE_FLAGS "-O0")

target_compile_definitions(gfsphysics PRIVATE -DNEMS_GSM)
target_link_libraries(gfsphysics PRIVATE NetCDF::NetCDF_Fortran)

set_target_properties(gfsphysics PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(gfsphysics PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
