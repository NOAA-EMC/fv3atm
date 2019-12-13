#!/usr/bin/env python

# CCPP prebuild config for GFDL Finite-Volume Cubed-Sphere Model (FV3)


###############################################################################
# Definitions                                                                 #
###############################################################################

HOST_MODEL_IDENTIFIER = "FV3"

# Add all files with metadata tables on the host model side,
# relative to basedir = top-level directory of host model
VARIABLE_DEFINITION_FILES = [
    'FV3/ccpp/physics/physics/machine.F',
    'FV3/ccpp/physics/physics/radsw_param.f',
    'FV3/ccpp/physics/physics/radlw_param.f',
    'FV3/gfsphysics/CCPP_layer/CCPP_typedefs.F90',
    'FV3/gfsphysics/GFS_layer/GFS_typedefs.F90',
    'FV3/gfsphysics/CCPP_layer/CCPP_data.F90',
    ]

TYPEDEFS_NEW_METADATA = {
    'ccpp_types' : {
        'ccpp_types' : '',
        'ccpp_t' : 'cdata',
        },
    'machine' : {
        'machine' : '',
        },
    'module_radlw_parameters' : {
        'module_radsw_parameters' : '',
        },
    'module_radlw_parameters' : {
        'module_radlw_parameters' : '',
        },
    'CCPP_typedefs' : {
        'CCPP_interstitial_type' : 'CCPP_interstitial',
        'CCPP_typedefs' : '',
        },
    'CCPP_data' : {
        'CCPP_data' : '',
        },
    'GFS_typedefs' : {
        'GFS_control_type'      : 'GFS_Control',
        'GFS_interstitial_type' : 'GFS_Interstitial(cdata%thrd_no)',
        'GFS_data_type'         : 'GFS_Data(cdata%blk_no)',
        'GFS_diag_type'         : 'GFS_Data(cdata%blk_no)%Intdiag',
        'GFS_tbd_type'          : 'GFS_Data(cdata%blk_no)%Tbd',
        'GFS_sfcprop_type'      : 'GFS_Data(cdata%blk_no)%Sfcprop',
        'GFS_coupling_type'     : 'GFS_Data(cdata%blk_no)%Coupling',
        'GFS_statein_type'      : 'GFS_Data(cdata%blk_no)%Statein',
        'GFS_cldprop_type'      : 'GFS_Data(cdata%blk_no)%Cldprop',
        'GFS_radtend_type'      : 'GFS_Data(cdata%blk_no)%Radtend',
        'GFS_grid_type'         : 'GFS_Data(cdata%blk_no)%Grid',
        'GFS_stateout_type'     : 'GFS_Data(cdata%blk_no)%Stateout',
        'GFS_typedefs' : '',
        },
    }

# Add all physics scheme dependencies relative to basedir - note that the CCPP
# rules stipulate that dependencies are not shared between the schemes!
SCHEME_FILES_DEPENDENCIES = [
    'FV3/ccpp/physics/physics/GFDL_parse_tracers.F90',
    'FV3/ccpp/physics/physics/aer_cloud.F',
    'FV3/ccpp/physics/physics/aerclm_def.F',
    'FV3/ccpp/physics/physics/aerinterp.F90',
    'FV3/ccpp/physics/physics/calpreciptype.f90',
    'FV3/ccpp/physics/physics/cldwat2m_micro.F',
    'FV3/ccpp/physics/physics/cldmacro.F',
    'FV3/ccpp/physics/physics/date_def.f',
    'FV3/ccpp/physics/physics/funcphys.f90',
    'FV3/ccpp/physics/physics/gcycle.F90',
    'FV3/ccpp/physics/physics/gfs_phy_tracer_config.F',
    'FV3/ccpp/physics/physics/gocart_tracer_config_stub.f',
    'FV3/ccpp/physics/physics/h2o_def.f',
    'FV3/ccpp/physics/physics/h2ointerp.f90',
    'FV3/ccpp/physics/physics/iccn_def.F',
    'FV3/ccpp/physics/physics/iccninterp.F90',
    'FV3/ccpp/physics/physics/iounitdef.f',
    'FV3/ccpp/physics/physics/machine.F',
    'FV3/ccpp/physics/physics/mersenne_twister.f',
    'FV3/ccpp/physics/physics/mfpbl.f',
    'FV3/ccpp/physics/physics/micro_mg_utils.F90',
    'FV3/ccpp/physics/physics/micro_mg2_0.F90',
    'FV3/ccpp/physics/physics/micro_mg3_0.F90',
    'FV3/ccpp/physics/physics/module_bfmicrophysics.f',
    'FV3/ccpp/physics/physics/multi_gases.F90',
    'FV3/ccpp/physics/physics/module_gfdl_cloud_microphys.F90',
    'FV3/ccpp/physics/physics/module_nst_model.f90',
    'FV3/ccpp/physics/physics/module_nst_parameters.f90',
    'FV3/ccpp/physics/physics/module_nst_water_prop.f90',
    'FV3/ccpp/physics/physics/module_mp_radar.F90',
    'FV3/ccpp/physics/physics/module_mp_thompson.F90',
    'FV3/ccpp/physics/physics/module_mp_thompson_make_number_concentrations.F90',
    'FV3/ccpp/physics/physics/module_bl_mynn.F90',
    'FV3/ccpp/physics/physics/module_sf_mynn.F90',
    'FV3/ccpp/physics/physics/module_SF_JSFC.F90',
    'FV3/ccpp/physics/physics/module_BL_MYJPBL.F90',
    'FV3/ccpp/physics/physics/module_sf_noahmp_glacier.f90',
    'FV3/ccpp/physics/physics/module_sf_noahmplsm.f90',
    'FV3/ccpp/physics/physics/cires_ugwp_module.F90',
    'FV3/ccpp/physics/physics/ugwp_driver_v0.F',
    'FV3/ccpp/physics/physics/cires_ugwp_triggers.F90',
    'FV3/ccpp/physics/physics/cires_ugwp_initialize.F90',
    'FV3/ccpp/physics/physics/cires_ugwp_solvers.F90',
    'FV3/ccpp/physics/physics/cires_ugwp_utils.F90',
    'FV3/ccpp/physics/physics/cires_orowam2017.f',
    'FV3/ccpp/physics/physics/cires_vert_lsatdis.F90',
    'FV3/ccpp/physics/physics/cires_vert_orodis.F90',
    'FV3/ccpp/physics/physics/cires_vert_wmsdis.F90',
    'FV3/ccpp/physics/physics/namelist_soilveg.f',
    'FV3/ccpp/physics/physics/mfpblt.f',
    'FV3/ccpp/physics/physics/mfpbltq.f',
    'FV3/ccpp/physics/physics/mfscu.f',
    'FV3/ccpp/physics/physics/mfscuq.f',
    'FV3/ccpp/physics/physics/noahmp_tables.f90',
    'FV3/ccpp/physics/physics/num_parthds.F',
    'FV3/ccpp/physics/physics/ozne_def.f',
    'FV3/ccpp/physics/physics/ozinterp.f90',
    'FV3/ccpp/physics/physics/physcons.F90',
    'FV3/ccpp/physics/physics/physparam.f',
    'FV3/ccpp/physics/physics/radcons.f90',
    'FV3/ccpp/physics/physics/radiation_aerosols.f',
    'FV3/ccpp/physics/physics/radiation_astronomy.f',
    'FV3/ccpp/physics/physics/radiation_clouds.f',
    'FV3/ccpp/physics/physics/radiation_gases.f',
    'FV3/ccpp/physics/physics/radiation_surface.f',
    'FV3/ccpp/physics/physics/radlw_datatb.f',
    'FV3/ccpp/physics/physics/radlw_param.f',
    'FV3/ccpp/physics/physics/radsw_datatb.f',
    'FV3/ccpp/physics/physics/radsw_param.f',
    'FV3/ccpp/physics/physics/samfaerosols.F',
    'FV3/ccpp/physics/physics/sfcsub.F',
    'FV3/ccpp/physics/physics/sflx.f',
    'FV3/ccpp/physics/physics/set_soilveg.f',
    'FV3/ccpp/physics/physics/surface_perturbation.F90',
    'FV3/ccpp/physics/physics/cu_gf_deep.F90',
    'FV3/ccpp/physics/physics/cu_gf_sh.F90',
    'FV3/ccpp/physics/physics/tridi.f',
    'FV3/ccpp/physics/physics/wv_saturation.F',
    'FV3/ccpp/physics/physics/module_sf_ruclsm.F90',
    'FV3/ccpp/physics/physics/namelist_soilveg_ruc.F90',
    'FV3/ccpp/physics/physics/set_soilveg_ruc.F90',
    'FV3/ccpp/physics/physics/module_soil_pre.F90',
    # derived data type definitions
    'FV3/gfsphysics/GFS_layer/GFS_typedefs.F90',
    'FV3/gfsphysics/CCPP_layer/CCPP_typedefs.F90',
    ]

# Add all physics scheme files relative to basedir
SCHEME_FILES = {
    # Relative path to source (from where ccpp_prebuild.py is called) : [ list of physics sets in which scheme may be called ];
    # current restrictions are that each scheme can only belong to one physics set, and all schemes within one group in the
    # suite definition file have to belong to the same physics set
    'FV3/ccpp/physics/physics/GFS_DCNV_generic.F90'              : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/GFS_GWD_generic.F90'               : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/GFS_MP_generic.F90'                : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/GFS_PBL_generic.F90'               : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/GFS_SCNV_generic.F90'              : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/GFS_debug.F90'                     : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/GFS_phys_time_vary.fv3.F90'        : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/GFS_rad_time_vary.fv3.F90'         : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/GFS_rrtmg_post.F90'                : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/GFS_rrtmg_pre.F90'                 : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/GFS_rrtmg_setup.F90'               : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/GFS_stochastics.F90'               : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/GFS_suite_interstitial.F90'        : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/GFS_surface_generic.F90'           : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/GFS_surface_composites.F90'        : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/GFS_surface_loop_control.F90'      : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/GFS_time_vary_pre.fv3.F90'         : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/cires_ugwp.F90'                    : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/cires_ugwp_post.F90'               : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/cnvc90.f'                          : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/cs_conv.F90'                       : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/cs_conv_aw_adj.F90'                : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/cu_ntiedtke_pre.F90'               : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/cu_ntiedtke.F90'                   : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/cu_ntiedtke_post.F90'              : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/dcyc2.f'                           : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/drag_suite.F90'                    : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/gcm_shoc.F90'                      : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/get_prs_fv3.F90'                   : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/gfdl_cloud_microphys.F90'          : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/gfdl_fv_sat_adj.F90'               : [ 'fast_physics' ],
    'FV3/ccpp/physics/physics/gscond.f'                          : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/gwdc.f'                            : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/gwdps.f'                           : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/h2ophys.f'                         : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/samfdeepcnv.f'                     : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/samfshalcnv.f'                     : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/maximum_hourly_diagnostics.F90'    : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/m_micro.F90'                       : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/m_micro_interstitial.F90'          : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/cu_gf_driver_pre.F90'              : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/cu_gf_driver.F90'                  : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/cu_gf_driver_post.F90'             : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/moninedmf.f'                       : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/moninshoc.f'                       : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/satmedmfvdif.F'                    : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/satmedmfvdifq.F'                   : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/shinhongvdif.F90'                  : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/ysuvdif.F90'                       : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/module_MYNNPBL_wrapper.F90'        : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/module_MYNNSFC_wrapper.F90'        : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/module_MYNNrad_pre.F90'            : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/module_MYNNrad_post.F90'           : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/module_MYJSFC_wrapper.F90'         : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/module_MYJPBL_wrapper.F90'         : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/mp_thompson_pre.F90'               : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/mp_thompson.F90'                   : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/mp_thompson_post.F90'              : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/ozphys.f'                          : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/ozphys_2015.f'                     : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/precpd.f'                          : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/radlw_main.f'                      : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/radsw_main.f'                      : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/rayleigh_damp.f'                   : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/rrtmg_lw_post.F90'                 : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/rrtmg_lw_pre.F90'                  : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/rrtmg_sw_post.F90'                 : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/rrtmg_sw_pre.F90'                  : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/sfc_diag.f'                        : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/sfc_diag_post.F90'                 : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/sfc_drv_ruc.F90'                   : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/lsm_ruc_sfc_sice_interstitial.F90' : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/sfc_cice.f'                        : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/sfc_diff.f'                        : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/sfc_drv.f'                         : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/sfc_noahmp_drv.f'                  : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/sfc_nst.f'                         : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/sfc_ocean.F'                       : [ 'slow_physics' ],
    'FV3/ccpp/physics/physics/sfc_sice.f'                        : [ 'slow_physics' ],
    # for testing the <init> and <finalize> sections
    'FV3/ccpp/physics/physics/GFS_suite_init_finalize_test.F90'  : [ 'slow_physics' ],
    }

# Default build dir, relative to current working directory,
# if not specified as command-line argument
DEFAULT_BUILD_DIR = 'FV3'

# Auto-generated makefile/cmakefile snippets that contain all schemes
SCHEMES_MAKEFILE   = '{build_dir}/ccpp/physics/CCPP_SCHEMES.mk'
SCHEMES_CMAKEFILE  = '{build_dir}/ccpp/physics/CCPP_SCHEMES.cmake'
SCHEMES_SOURCEFILE = '{build_dir}/ccpp/physics/CCPP_SCHEMES.sh'

# CCPP host cap in which to insert the ccpp_field_add statements;
# determines the directory to place ccpp_{modules,fields}.inc
TARGET_FILES = [
    'FV3/atmos_cubed_sphere/driver/fvGFS/atmosphere.F90',
    'FV3/ccpp/driver/CCPP_Driver.F90',
    ]

# Auto-generated makefile/cmakefile snippets that contain all caps
CAPS_MAKEFILE   = '{build_dir}/ccpp/physics/CCPP_CAPS.mk'
CAPS_CMAKEFILE  = '{build_dir}/ccpp/physics/CCPP_CAPS.cmake'
CAPS_SOURCEFILE = '{build_dir}/ccpp/physics/CCPP_CAPS.sh'

# Directory where to put all auto-generated physics caps
CAPS_DIR = '{build_dir}/ccpp/physics'

# Directory where the suite definition files are stored
SUITES_DIR = 'FV3/ccpp/suites'

# Optional arguments - only required for schemes that use
# optional arguments. ccpp_prebuild.py will throw an exception
# if it encounters a scheme subroutine with optional arguments
# if no entry is made here. Possible values are: 'all', 'none',
# or a list of standard_names: [ 'var1', 'var3' ].
OPTIONAL_ARGUMENTS = {
    'rrtmg_sw' : {
        'rrtmg_sw_run' : [
            'tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step',
            'components_of_surface_downward_shortwave_fluxes',
            'cloud_liquid_water_path',
            'mean_effective_radius_for_liquid_cloud',
            'cloud_ice_water_path',
            'mean_effective_radius_for_ice_cloud',
            'cloud_rain_water_path',
            'mean_effective_radius_for_rain_drop',
            'cloud_snow_water_path',
            'mean_effective_radius_for_snow_flake',
            ],
        },
    'rrtmg_lw' : {
        'rrtmg_lw_run' : [
            'tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step',
            'cloud_liquid_water_path',
            'mean_effective_radius_for_liquid_cloud',
            'cloud_ice_water_path',
            'mean_effective_radius_for_ice_cloud',
            'cloud_rain_water_path',
            'mean_effective_radius_for_rain_drop',
            'cloud_snow_water_path',
            'mean_effective_radius_for_snow_flake',
            ],
        },
    'mp_thompson' : {
        'mp_thompson_init' : [
            'water_friendly_aerosol_number_concentration',
            'ice_friendly_aerosol_number_concentration',
            'tendency_of_water_friendly_aerosols_at_surface',
            'tendency_of_ice_friendly_aerosols_at_surface',
            ],
        'mp_thompson_run' : [
            'cloud_droplet_number_concentration_updated_by_physics',
            'water_friendly_aerosol_number_concentration_updated_by_physics',
            'ice_friendly_aerosol_number_concentration_updated_by_physics',
            'tendency_of_water_friendly_aerosols_at_surface',
            'tendency_of_ice_friendly_aerosols_at_surface',
            'mean_effective_radius_for_liquid_cloud',
            'mean_effective_radius_for_ice_cloud',
            'mean_effective_radius_for_snow_flake',
            ],
        },
    'mp_thompson_pre' : {
        'mp_thompson_pre_run' : [
            'cloud_droplet_number_concentration_updated_by_physics',
            'water_friendly_aerosol_number_concentration_updated_by_physics',
            'ice_friendly_aerosol_number_concentration_updated_by_physics',
            'tendency_of_water_friendly_aerosols_at_surface',
            'tendency_of_ice_friendly_aerosols_at_surface',
            ],
        },
    #'subroutine_name_1' : 'all',
    #'subroutine_name_2' : 'none',
    #'subroutine_name_2' : [ 'var1', 'var3'],
    }

# Names of Fortran include files in the host model cap (do not change);
# both files will be written to the directory of each target file, only
# used by the dynamic builds
MODULE_INCLUDE_FILE = 'ccpp_modules_{set}.inc'
FIELDS_INCLUDE_FILE = 'ccpp_fields_{set}.inc'

# Directory where to write static API to
STATIC_API_DIR = '{build_dir}/ccpp/physics'
STATIC_API_SRCFILE = '{build_dir}/ccpp/physics/CCPP_STATIC_API.sh'

# Directory for writing HTML pages generated from metadata files
METADATA_HTML_OUTPUT_DIR = '{build_dir}/ccpp/physics/physics/docs'

# HTML document containing the model-defined CCPP variables
HTML_VARTABLE_FILE = '{build_dir}/ccpp/physics/CCPP_VARIABLES_FV3.html'

# LaTeX document containing the provided vs requested CCPP variables
LATEX_VARTABLE_FILE = '{build_dir}/ccpp/framework/doc/DevelopersGuide/CCPP_VARIABLES_FV3.tex'


###############################################################################
# Template code to generate include files                                     #
###############################################################################

# Name of the CCPP data structure in the host model cap;
# in the case of FV3, this is a 2-dimensional array with
# the number of blocks as the first and the number of
# OpenMP threads as the second dimension; nb is the loop
# index for the current block, nt for the current thread.
# Internally, the model uses an associate construct to
# reference cdata(nb,nt) with cdata (recommended).
CCPP_DATA_STRUCTURE = 'cdata'
