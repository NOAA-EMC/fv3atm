SHELL = /bin/sh

inside_nems := $(wildcard ../../../conf/configure.nems)
ifneq ($(strip $(inside_nems)),)
    include ../../../conf/configure.nems
else
    exist_configure_fv3 := $(wildcard ../conf/configure.fv3)
    ifneq ($(strip $(exist_configure_fv3)),)
        include ../conf/configure.fv3
    else
        $(error "../conf/configure.fv3 file is missing. Run ./configure")
    endif
    $(info )
    $(info Build standalone FV3 gfsphysics ...)
    $(info )
endif

LIBRARY  = libgfsphys.a

FFLAGS   += -I$(FMS_DIR) -I../cpl

CPPDEFS += -DNEW_TAUCTMAX -DSMALL_PE -DNEMS_GSM -DINTERNAL_FILE_NML

# CCPP dynamic and static builds
ifneq (,$(findstring CCPP,$(CPPDEFS)))
# Set flags for 32-bit dynamics build
ifeq ($(DYN32),Y)
CPPDEFS += -DOVERLOAD_R4
endif
# Set CCPP static api for static build
ifneq (,$(findstring STATIC,$(CPPDEFS)))
    CCPP_STATIC_API = ../ccpp/physics/ccpp_static_api.F90
else
    CCPP_STATIC_API = 
endif
# Set physics source files
# DH* 20190916 - can remove namelist_soilveg and set_soilveg
# once calculation of snocvr is moved from FV3GFS_io.F90/GFS_driver.F90
# to CCPP (GFS_phys_time_vary.fv3.F90) *DH
SRCS_f = \
	   ./physics/mersenne_twister.f                                              \
	   ./physics/namelist_soilveg.f                                              \
	   ./physics/physparam.f                                                     \
	   ./physics/radlw_param.f                                                   \
	   ./physics/radsw_param.f                                                   \
	   ./physics/set_soilveg.f
SRCS_f90 = \
	     ./physics/noahmp_tables.f90
SRCS_F = \
	   ./physics/machine.F
SRCS_F90 = \
	   ./physics/GFDL_parse_tracers.F90                                          \
	   ./physics/physcons.F90                                                    \
	   ./CCPP_layer/CCPP_typedefs.F90                                            \
	   ./CCPP_layer/CCPP_data.F90                                                \
	   $(CCPP_STATIC_API)                                                        \
	   ./GFS_layer/GFS_abstraction_layer.F90                                     \
	   ./GFS_layer/GFS_diagnostics.F90                                           \
	   ./GFS_layer/GFS_driver.F90                                                \
	   ./GFS_layer/GFS_restart.F90                                               \
	   ./GFS_layer/GFS_typedefs.F90
# non-CCPP build
else
SRCS_f =  \
	   ./physics/cnvc90.f                                                        \
	   ./physics/co2hc.f                                                         \
	   ./physics/date_def.f                                                      \
	   ./physics/dcyc2.f                                                         \
	   ./physics/dcyc2.pre.rad.f                                                 \
	   ./physics/efield.f                                                        \
	   ./physics/get_prs.f                                                       \
	   ./physics/gocart_tracer_config_stub.f                                     \
	   ./physics/gscond.f                                                        \
	   ./physics/gscondp.f                                                       \
	   ./physics/gwdc.f                                                          \
	   ./physics/gwdps.f                                                         \
	   ./physics/ugwp_driver_v0.f                                                \
	   ./physics/cires_orowam2017.f                                              \
	   ./physics/h2o_def.f                                                       \
	   ./physics/h2oc.f                                                          \
	   ./physics/h2ohdc.f                                                        \
	   ./physics/h2ophys.f                                                       \
	   ./physics/ideaca.f                                                        \
	   ./physics/idea_co2.f                                                      \
	   ./physics/idea_composition.f                                              \
	   ./physics/idea_dissipation.f                                              \
	   ./physics/idea_h2o.f                                                      \
	   ./physics/idea_ion.f                                                      \
	   ./physics/idea_o2_o3.f                                                    \
	   ./physics/idea_phys.f                                                     \
	   ./physics/idea_solar_heating.f                                            \
	   ./physics/idea_tracer.f                                                   \
	   ./physics/iounitdef.f                                                     \
	   ./physics/lrgsclr.f                                                       \
	   ./physics/mersenne_twister.f                                              \
	   ./physics/mfpbl.f                                                         \
	   ./physics/mfpblt.f                                                        \
           ./physics/mfpbltq.f                                                       \
	   ./physics/mfscu.f                                                         \
           ./physics/mfscuq.f                                                        \
	   ./physics/module_bfmicrophysics.f                                         \
	   ./physics/moninedmf.f                                                     \
	   ./physics/moninedmf_hafs.f                                                \
	   ./physics/moninp.f                                                        \
	   ./physics/moninp1.f                                                       \
	   ./physics/moninq.f                                                        \
	   ./physics/moninq1.f                                                       \
	   ./physics/moninshoc.f                                                     \
	   ./physics/mstadb.f                                                        \
	   ./physics/mstadbtn.f                                                      \
	   ./physics/mstadbtn2.f                                                     \
	   ./physics/mstcnv.f                                                        \
	   ./physics/namelist_soilveg.f                                              \
	   ./physics/ozne_def.f                                                      \
	   ./physics/iccn_def.f                                                      \
	   ./physics/aerclm_def.f                                                    \
	   ./physics/ozphys.f                                                        \
	   ./physics/ozphys_2015.f                                                   \
	   ./physics/physparam.f                                                     \
	   ./physics/precpd.f                                                        \
	   ./physics/precpd_shoc.f                                                   \
	   ./physics/precpdp.f                                                       \
	   ./physics/precpd_shoc.f                                                   \
	   ./physics/progt2.f                                                        \
	   ./physics/progtm_module.f                                                 \
	   ./physics/rad_initialize.f                                                \
	   ./physics/radiation_aerosols.f                                            \
	   ./physics/radiation_astronomy.f                                           \
	   ./physics/radiation_clouds.f                                              \
	   ./physics/radiation_gases.f                                               \
	   ./physics/radiation_surface.f                                             \
	   ./physics/radlw_datatb.f                                                  \
	   ./physics/radlw_main.f                                                    \
	   ./physics/radlw_param.f                                                   \
	   ./physics/radsw_datatb.f                                                  \
	   ./physics/radsw_main.f                                                    \
	   ./physics/radsw_param.f                                                   \
	   ./physics/rascnvv2.f                                                      \
	   ./physics/rayleigh_damp.f                                                 \
	   ./physics/rayleigh_damp_mesopause.f                                       \
	   ./physics/samfaerosols.f                                                  \
	   ./physics/samfdeepcnv.f                                                   \
	   ./physics/samfshalcnv.f                                                   \
	   ./physics/sascnv.f                                                        \
	   ./physics/sascnvn.f                                                       \
	   ./physics/satmedmfvdif.f                                                  \
	   ./physics/satmedmfvdifq.f                                                 \
	   ./physics/set_soilveg.f                                                   \
	   ./physics/sfc_cice.f                                                      \
	   ./physics/sfc_diag.f                                                      \
	   ./physics/sfc_diff.f                                                      \
	   ./physics/sfc_drv.f                                                       \
	   ./physics/sfc_noahmp_drv.f                                                \
	   ./physics/sfc_nst.f                                                       \
	   ./physics/sfc_ocean.f                                                     \
	   ./physics/sfc_sice.f                                                      \
	   ./physics/sflx.f                                                          \
	   ./physics/shalcnv.f                                                       \
	   ./physics/shalcv.f                                                        \
	   ./physics/shalcv_1lyr.f                                                   \
	   ./physics/shalcv_fixdp.f                                                  \
	   ./physics/shalcv_opr.f                                                    \
	   ./physics/tracer_const_h.f                                                \
	   ./physics/tridi2t3.f

SRCS_f90 = \
	   ./physics/calpreciptype.f90                                               \
	   ./physics/funcphys.f90                                                    \
	   ./physics/gcm_shoc.f90                                                    \
	   ./physics/get_prs_fv3.f90                                                 \
	   ./physics/h2ointerp.f90                                                   \
	   ./physics/module_nst_model.f90                                            \
	   ./physics/module_nst_parameters.f90                                       \
	   ./physics/module_nst_water_prop.f90                                       \
	   ./physics/ozinterp.f90                                                    \
	   ./physics/module_wrf_utl.f90                                              \
	   ./physics/noahmp_tables.f90                                               \
	   ./physics/module_sf_noahmplsm.f90                                         \
	   ./physics/module_sf_noahmp_glacier.f90                                    \
	   ./physics/iccninterp.f90                                                  \
	   ./physics/aerinterp.f90                                                   \
	   ./physics/wam_f107_kp_mod.f90

SRCS_F = \
	   ./physics/aer_cloud.F                                                     \
	   ./physics/cldmacro.F                                                      \
	   ./physics/cldwat2m_micro.F                                                \
	   ./physics/gfs_phy_tracer_config.F                                         \
	   ./physics/machine.F                                                       \
	   ./physics/num_parthds.F                                                   \
	   ./physics/sfcsub.F                                                        \
	   ./physics/wv_saturation.F

SRCS_F90 = \
	   ./physics/GFDL_parse_tracers.F90                                          \
	   ./physics/gcycle.F90                                                      \
  	   ./physics/cires_ugwp_initialize.F90                                       \
	   ./physics/cires_ugwp_module.F90                                           \
	   ./physics/cires_ugwp_utils.F90                                            \
	   ./physics/cires_ugwp_triggers.F90                                         \
	   ./physics/cires_ugwp_solvers.F90                                          \
	   ./physics/cires_vert_lsatdis.F90                                          \
	   ./physics/cires_vert_orodis.F90                                           \
	   ./physics/cires_vert_wmsdis.F90                                           \
	   ./physics/gfdl_cloud_microphys.F90                                        \
	   ./physics/micro_mg_utils.F90                                              \
	   ./physics/micro_mg2_0.F90                                                 \
	   ./physics/micro_mg3_0.F90                                                 \
	   ./physics/m_micro_driver.F90                                              \
	   ./physics/cs_conv.F90                                                     \
	   ./physics/GFS_debug.F90                                                   \
	   ./physics/module_mp_radar.F90                                             \
	   ./physics/module_mp_thompson_gfs.F90                                      \
	   ./physics/module_mp_wsm6_fv3.F90                                          \
	   ./physics/physcons.F90                                                    \
	   ./physics/surface_perturbation.F90                                        \
	   ./GFS_layer/GFS_abstraction_layer.F90                                     \
	   ./GFS_layer/GFS_diagnostics.F90                                           \
	   ./GFS_layer/GFS_driver.F90                                                \
	   ./GFS_layer/GFS_physics_driver.F90                                        \
	   ./GFS_layer/GFS_radiation_driver.F90                                      \
	   ./GFS_layer/GFS_restart.F90                                               \
	   ./GFS_layer/GFS_typedefs.F90
endif

SRCS_c   =

DEPEND_FILES = $(SRCS_f) $(SRCS_f90) $(SRCS_F) $(SRCS_F90)

OBJS_f   = $(SRCS_f:.f=.o)
OBJS_f90 = $(SRCS_f90:.f90=.o)
OBJS_F   = $(SRCS_F:.F=.o)
OBJS_F90 = $(SRCS_F90:.F90=.o)
OBJS_c   = $(SRCS_c:.c=.o)

OBJS = $(OBJS_f) $(OBJS_f90) $(OBJS_F) $(OBJS_F90) $(OBJS_c)

all default: depend $(LIBRARY)

$(LIBRARY): $(OBJS)
	$(AR) $(ARFLAGS) $@ $?

# this is the place to override default (implicit) compilation rules
# and create specific (explicit) rules

# Reduce optimization (substitute (x)CORE-AVX2, (x)CORE-AVX512 or xHOST with (x)CORE-AVX-I) for radiation_aerosols.f
FFLAGS_LOPT1=$(subst CORE-AVX512,CORE-AVX-I,\
               $(subst CORE-AVX2,CORE-AVX-I,\
                 $(subst xHOST,xCORE-AVX-I,$(FFLAGS))))
./physics/radiation_aerosols.o : ./physics/radiation_aerosols.f
	$(FC) $(CPPDEFS) $(FFLAGS_LOPT1) $(OTHER_FFLAGS) -c $< -o $@

# Reduce optimization (add -O0) for GFS_diagnsostics.F90
./GFS_layer/GFS_diagnostics.o : ./GFS_layer/GFS_diagnostics.F90
	$(FC) $(CPPDEFS) $(FFLAGS) $(OTHER_FFLAGS) -O0 -c $< -o $@

# Do preprocessing of the GFS drivers in two steps to be able
# to look at the actual .f90 file that is compiled
./GFS_layer/GFS_driver.o: ./GFS_layer/GFS_driver.F90
	$(CPP) $(CPPDEFS) $(CPPFLAGS) $< > $*.tmp.f90
	$(FC) $(FFLAGS) $(OTHER_FFLAGS) -c $*.tmp.f90 -o $@
./GFS_layer/GFS_physics_driver.o: ./GFS_layer/GFS_physics_driver.F90
	$(CPP) $(CPPDEFS) $(CPPFLAGS) $< > $*.tmp.f90
	$(FC) $(FFLAGS) $(OTHER_FFLAGS) -c $*.tmp.f90 -o $@
./GFS_layer/GFS_radiation_driver.o: ./GFS_layer/GFS_radiation_driver.F90
	$(CPP) $(CPPDEFS) $(CPPFLAGS) $< > $*.tmp.f90
	$(FC) $(FFLAGS) $(OTHER_FFLAGS) -c $*.tmp.f90 -o $@

.PHONY: clean
clean:
	@echo "Cleaning gfsphysics  ... "
	@echo
	$(RM) -f $(LIBRARY) *__genmod.f90 *.o */*.o *.mod *.i90 *.lst *.i depend */*.tmp.f90

MKDEPENDS = ../mkDepends.pl
include ../conf/make.rules

include ./depend

# do not include 'depend' file if the target contains string 'clean'
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
   -include depend
endif

