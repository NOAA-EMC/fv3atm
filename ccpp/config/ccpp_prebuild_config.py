#!/usr/bin/env python

# CCPP prebuild config for GFDL Finite-Volume Cubed-Sphere Model (FV3)


###############################################################################
# Definitions                                                                 #
###############################################################################

HOST_MODEL_IDENTIFIER = "FV3"

# Add all files with metadata tables on the host model side and in CCPP,
# relative to basedir = top-level directory of host model. This includes
# kind and type definitions used in CCPP physics. Also add any internal
# dependencies of these files to the list.
VARIABLE_DEFINITION_FILES = [
    # actual variable definition files
    'ccpp/physics/physics/machine.F',
    'ccpp/physics/physics/radsw_param.f',
    'ccpp/physics/physics/h2o_def.f',
    'ccpp/physics/physics/ozne_def.f',
    'ccpp/physics/physics/radlw_param.f',
    'ccpp/data/CCPP_typedefs.F90',
    'ccpp/data/GFS_typedefs.F90',
    'ccpp/data/CCPP_data.F90',
    'ccpp/physics/physics/rte-rrtmgp/rrtmgp/mo_gas_optics_rrtmgp.F90',
    'ccpp/physics/physics/rte-rrtmgp/rrtmgp/mo_gas_concentrations.F90',
    'ccpp/physics/physics/rte-rrtmgp/rte/mo_optical_props.F90',
    'ccpp/physics/physics/rte-rrtmgp/extensions/cloud_optics/mo_cloud_optics.F90',
    'ccpp/physics/physics/rte-rrtmgp/rte/mo_source_functions.F90',
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
    'mo_gas_concentrations' : {
        'ty_gas_concs' : '',
        'mo_gas_concentrations' : '',
        },
    'mo_gas_optics_rrtmgp' : {
        'ty_gas_optics_rrtmgp' : '',
        'mo_gas_optics_rrtmgp' : '',
        },
    'mo_optical_props' : {
        'ty_optical_props_1scl' : '',
        'ty_optical_props_2str' : '',
        'mo_optical_props' : '',
        },
    'mo_cloud_optics' : {
        'ty_cloud_optics' : '',
        'mo_cloud_optics' : '',
        },
    'mo_source_functions' : {
        'ty_source_func_lw' : '',
        'mo_source_functions' : '',
        },
    }

# Add all physics scheme files relative to basedir
SCHEME_FILES = [
    # Relative path to source (from where ccpp_prebuild.py is called) : [ list of physics sets in which scheme may be called ];
    # current restrictions are that each scheme can only belong to one physics set, and all schemes within one group in the
    # suite definition file have to belong to the same physics set
    'ccpp/physics/physics/GFS_DCNV_generic.F90',
    'ccpp/physics/physics/GFS_GWD_generic.F90',
    'ccpp/physics/physics/GFS_MP_generic.F90',
    'ccpp/physics/physics/GFS_PBL_generic.F90',
    'ccpp/physics/physics/GFS_SCNV_generic.F90',
    'ccpp/physics/physics/GFS_debug.F90',
    'ccpp/physics/physics/GFS_phys_time_vary.fv3.F90',
    'ccpp/physics/physics/GFS_rad_time_vary.fv3.F90',
    'ccpp/physics/physics/GFS_rrtmg_post.F90',
    'ccpp/physics/physics/GFS_rrtmg_pre.F90',
    'ccpp/physics/physics/GFS_rrtmg_setup.F90',
    'ccpp/physics/physics/GFS_stochastics.F90',
    'ccpp/physics/physics/GFS_suite_interstitial.F90',
    'ccpp/physics/physics/GFS_surface_generic.F90',
    'ccpp/physics/physics/GFS_surface_composites.F90',
    'ccpp/physics/physics/GFS_surface_loop_control.F90',
    'ccpp/physics/physics/GFS_time_vary_pre.fv3.F90',
    'ccpp/physics/physics/cires_ugwp.F90',
    'ccpp/physics/physics/cires_ugwp_post.F90',
    'ccpp/physics/physics/unified_ugwp.F90',
    'ccpp/physics/physics/unified_ugwp_post.F90',
    'ccpp/physics/physics/ugwpv1_gsldrag.F90',
    'ccpp/physics/physics/ugwpv1_gsldrag_post.F90',        
    'ccpp/physics/physics/cnvc90.f',
    'ccpp/physics/physics/cs_conv.F90',
    'ccpp/physics/physics/cs_conv_aw_adj.F90',
    'ccpp/physics/physics/cu_ntiedtke_pre.F90',
    'ccpp/physics/physics/cu_ntiedtke.F90',
    'ccpp/physics/physics/cu_ntiedtke_post.F90',
    'ccpp/physics/physics/dcyc2.f',
    'ccpp/physics/physics/drag_suite.F90',
    'ccpp/physics/physics/gcm_shoc.F90',
    'ccpp/physics/physics/get_prs_fv3.F90',
    'ccpp/physics/physics/gfdl_cloud_microphys.F90',
    'ccpp/physics/physics/gfdl_fv_sat_adj.F90',
    'ccpp/physics/physics/gfdl_sfc_layer.F90',
    'ccpp/physics/physics/gscond.f',
    'ccpp/physics/physics/gwdc.f',
    'ccpp/physics/physics/gwdps.f',
    'ccpp/physics/physics/h2ophys.f',
    'ccpp/physics/physics/samfdeepcnv.f',
    'ccpp/physics/physics/samfshalcnv.f',
    'ccpp/physics/physics/sascnvn.F',
    'ccpp/physics/physics/shalcnv.F',
    'ccpp/physics/physics/maximum_hourly_diagnostics.F90',
    'ccpp/physics/physics/m_micro.F90',
    'ccpp/physics/physics/m_micro_interstitial.F90',
    'ccpp/physics/physics/cu_gf_driver_pre.F90',
    'ccpp/physics/physics/cu_gf_driver.F90',
    'ccpp/physics/physics/cu_gf_driver_post.F90',
    'ccpp/physics/physics/moninedmf.f',
    'ccpp/physics/physics/moninshoc.f',
    'ccpp/physics/physics/satmedmfvdif.F',
    'ccpp/physics/physics/satmedmfvdifq.F',
    'ccpp/physics/physics/shinhongvdif.F90',
    'ccpp/physics/physics/ysuvdif.F90',
    'ccpp/physics/physics/module_MYNNPBL_wrapper.F90',
    'ccpp/physics/physics/module_MYNNSFC_wrapper.F90',
    'ccpp/physics/physics/module_SGSCloud_RadPre.F90',
    'ccpp/physics/physics/module_SGSCloud_RadPost.F90',
    'ccpp/physics/physics/module_MYJSFC_wrapper.F90',
    'ccpp/physics/physics/module_MYJPBL_wrapper.F90',
    'ccpp/physics/physics/mp_thompson_pre.F90',
    'ccpp/physics/physics/mp_thompson.F90',
    'ccpp/physics/physics/mp_thompson_post.F90',
    'ccpp/physics/physics/ozphys.f',
    'ccpp/physics/physics/ozphys_2015.f',
    'ccpp/physics/physics/precpd.f',
    'ccpp/physics/physics/phys_tend.F90',
    'ccpp/physics/physics/radlw_main.F90',
    'ccpp/physics/physics/radsw_main.F90',
    'ccpp/physics/physics/rascnv.F90',
    'ccpp/physics/physics/rayleigh_damp.f',
    'ccpp/physics/physics/rrtmg_lw_post.F90',
    'ccpp/physics/physics/rrtmg_lw_pre.F90',
    'ccpp/physics/physics/rrtmg_sw_post.F90',
    'ccpp/physics/physics/rrtmg_sw_pre.F90',
    'ccpp/physics/physics/sfc_diag.f',
    'ccpp/physics/physics/sfc_diag_post.F90',
    'ccpp/physics/physics/sfc_drv_ruc.F90',
    'ccpp/physics/physics/sfc_cice.f',
    'ccpp/physics/physics/sfc_diff.f',
    'ccpp/physics/physics/sfc_drv.f',
    'ccpp/physics/physics/sfc_noah_wrfv4_interstitial.F90',
    'ccpp/physics/physics/sfc_noah_wrfv4.F90',
    'ccpp/physics/physics/sfc_noahmp_drv.F90',
    'ccpp/physics/physics/flake_driver.F90',
    'ccpp/physics/physics/sfc_nst.f',
    'ccpp/physics/physics/sfc_ocean.F',
    'ccpp/physics/physics/sfc_sice.f',
    # HAFS FER_HIRES
    'ccpp/physics/physics/mp_fer_hires.F90',
    # RRTMGP
    'ccpp/physics/physics/rrtmgp_lw_gas_optics.F90',
    'ccpp/physics/physics/rrtmgp_lw_cloud_optics.F90',
    'ccpp/physics/physics/rrtmgp_sw_gas_optics.F90',
    'ccpp/physics/physics/rrtmgp_sw_cloud_optics.F90',
    'ccpp/physics/physics/rrtmgp_sw_aerosol_optics.F90',
    'ccpp/physics/physics/rrtmgp_lw_rte.F90',
    'ccpp/physics/physics/rrtmgp_sw_rte.F90',
    'ccpp/physics/physics/rrtmgp_lw_aerosol_optics.F90',
    'ccpp/physics/physics/GFS_rrtmgp_setup.F90',
    'ccpp/physics/physics/GFS_rrtmgp_pre.F90',
    'ccpp/physics/physics/rrtmgp_lw_pre.F90',
    'ccpp/physics/physics/GFS_rrtmgp_sw_pre.F90',
    'ccpp/physics/physics/GFS_rrtmgp_lw_post.F90',
    'ccpp/physics/physics/rrtmgp_lw_cloud_sampling.F90',
    'ccpp/physics/physics/rrtmgp_sw_cloud_sampling.F90',
    'ccpp/physics/physics/GFS_cloud_diagnostics.F90',
    'ccpp/physics/physics/GFS_rrtmgp_thompsonmp_pre.F90',
    'ccpp/physics/physics/GFS_rrtmgp_gfdlmp_pre.F90',
    'ccpp/physics/physics/GFS_rrtmgp_zhaocarr_pre.F90',
    'ccpp/physics/physics/GFS_rrtmgp_cloud_overlap_pre.F90',
    'ccpp/physics/physics/GFS_rrtmgp_sw_post.F90'
    ]

# Default build dir, relative to current working directory,
# if not specified as command-line argument
DEFAULT_BUILD_DIR = 'FV3'

# Auto-generated makefile/cmakefile snippets that contain all type definitions
TYPEDEFS_MAKEFILE   = '{build_dir}/ccpp/physics/CCPP_TYPEDEFS.mk'
TYPEDEFS_CMAKEFILE  = '{build_dir}/ccpp/physics/CCPP_TYPEDEFS.cmake'
TYPEDEFS_SOURCEFILE = '{build_dir}/ccpp/physics/CCPP_TYPEDEFS.sh'

# Auto-generated makefile/cmakefile snippets that contain all schemes
SCHEMES_MAKEFILE   = '{build_dir}/ccpp/physics/CCPP_SCHEMES.mk'
SCHEMES_CMAKEFILE  = '{build_dir}/ccpp/physics/CCPP_SCHEMES.cmake'
SCHEMES_SOURCEFILE = '{build_dir}/ccpp/physics/CCPP_SCHEMES.sh'

# Auto-generated makefile/cmakefile snippets that contain all caps
CAPS_MAKEFILE   = '{build_dir}/ccpp/physics/CCPP_CAPS.mk'
CAPS_CMAKEFILE  = '{build_dir}/ccpp/physics/CCPP_CAPS.cmake'
CAPS_SOURCEFILE = '{build_dir}/ccpp/physics/CCPP_CAPS.sh'

# Directory where to put all auto-generated physics caps
CAPS_DIR = '{build_dir}/ccpp/physics'

# Directory where the suite definition files are stored
SUITES_DIR = 'ccpp/suites'

# Optional arguments - only required for schemes that use
# optional arguments. ccpp_prebuild.py will throw an exception
# if it encounters a scheme subroutine with optional arguments
# if no entry is made here. Possible values are: 'all', 'none',
# or a list of standard_names: [ 'var1', 'var3' ].
OPTIONAL_ARGUMENTS = {
    'rrtmg_sw' : {
        'rrtmg_sw_run' : [
            'tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step_and_radiation_levels',
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
            'tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step_and_radiation_levels',
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
            'cloud_droplet_number_concentration',
            'water_friendly_aerosol_number_concentration',
            'ice_friendly_aerosol_number_concentration',
            'tendency_of_water_friendly_aerosols_at_surface',
            'tendency_of_ice_friendly_aerosols_at_surface',
            # DH* 2020-06-01: turn off calculation of effective radii, now done in GFS_rrtmg_pre
            #'effective_radius_of_stratiform_cloud_liquid_water_particle_in_um',
            #'effective_radius_of_stratiform_cloud_ice_particle_in_um',
            #'effective_radius_of_stratiform_cloud_snow_particle_in_um',
            # *DH 2020-06-01
            ],
        'mp_thompson_run' : [
            'cloud_droplet_number_concentration_updated_by_physics',
            'water_friendly_aerosol_number_concentration_updated_by_physics',
            'ice_friendly_aerosol_number_concentration_updated_by_physics',
            'tendency_of_water_friendly_aerosols_at_surface',
            'tendency_of_ice_friendly_aerosols_at_surface',
            # DH* 2020-06-01: turn off calculation of effective radii, now done in GFS_rrtmg_pre
            #'effective_radius_of_stratiform_cloud_liquid_water_particle_in_um',
            #'effective_radius_of_stratiform_cloud_ice_particle_in_um',
            #'effective_radius_of_stratiform_cloud_snow_particle_in_um',
            # *DH 2020-06-01
            ],
        },
    'mp_fer_hires' : {
        'mp_fer_hires_init' : [
            'fraction_of_ice_water_cloud',
            'fraction_of_rain_water_cloud',
            'rime_factor',
            ],
        },
    'rrtmgp_sw_rte' : {
         'rrtmgp_sw_rte_run' : [
             'components_of_surface_downward_shortwave_fluxes',
             ],
         },
    'GFS_rrtmgp_sw_post' : {
         'GFS_rrtmgp_sw_post_run' : [
             'tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step',
             'components_of_surface_downward_shortwave_fluxes',
             ],
         },
    'GFS_rrtmgp_lw_post' : {
         'GFS_rrtmgp_lw_post_run' : [
             'tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step',
             ],
         },
    #'subroutine_name_1' : 'all',
    #'subroutine_name_2' : 'none',
    #'subroutine_name_2' : [ 'var1', 'var3'],
    }

# Directory where to write static API to
STATIC_API_DIR = '{build_dir}/ccpp/physics'
STATIC_API_SRCFILE = '{build_dir}/ccpp/physics/CCPP_STATIC_API.sh'

# Directory for writing HTML pages generated from metadata files
METADATA_HTML_OUTPUT_DIR = '{build_dir}/ccpp/physics/physics/docs'

# HTML document containing the model-defined CCPP variables
HTML_VARTABLE_FILE = '{build_dir}/ccpp/physics/CCPP_VARIABLES_FV3.html'

# LaTeX document containing the provided vs requested CCPP variables
LATEX_VARTABLE_FILE = '{build_dir}/ccpp/framework/doc/DevelopersGuide/CCPP_VARIABLES_FV3.tex'
