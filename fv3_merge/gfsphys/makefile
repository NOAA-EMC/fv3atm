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
    $(info Build standalone FV3 gfsphys ...)
    $(info )
endif

LIBRARY  = libgfsphys.a

FFLAGS   += -I../fms -I../fms/include

SRCS_f   =  \
           ./cnvc90.f                                                        \
           ./co2hc.f                                                         \
           ./date_def.f                                                      \
           ./dcyc2.f                                                         \
           ./dcyc2.pre.rad.f                                                 \
           ./efield.f                                                        \
           ./get_prs.f                                                       \
           ./gfs_phy_tracer_config.f                                         \
           ./gocart_tracer_config_stub.f                                     \
           ./grrad.f                                                         \
           ./gscond.f                                                        \
           ./gscondp.f                                                       \
           ./gsmddrive.f                                                     \
           ./gwdps.f                                                         \
           ./h2oc.f                                                          \
           ./h2ohdc.f                                                        \
           ./ideaca.f                                                        \
           ./idea_co2.f                                                      \
           ./idea_composition.f                                              \
           ./idea_dissipation.f                                              \
           ./idea_h2o.f                                                      \
           ./idea_ion.f                                                      \
           ./idea_o2_o3.f                                                    \
           ./idea_phys.f                                                     \
           ./idea_solar_heating.f                                            \
           ./idea_tracer.f                                                   \
           ./iounitdef.f                                                     \
           ./lrgsclr.f                                                       \
           ./mersenne_twister.f                                              \
           ./mfpbl.f                                                         \
           ./module_bfmicrophysics.f                                         \
           ./moninedmf.f                                                     \
           ./moninp1.f                                                       \
           ./moninp.f                                                        \
           ./moninq1.f                                                       \
           ./moninq.f                                                        \
           ./moninshoc.f                                                     \
           ./mstadb.f                                                        \
           ./mstadbtn.f                                                      \
           ./mstcnv.f                                                        \
           ./namelist_soilveg.f                                              \
           ./ozne_def.f                                                      \
           ./ozphys.f                                                        \
           ./physparam.f                                                     \
           ./precpd.f                                                        \
           ./precpdp.f                                                       \
           ./precpd_shoc.f                                                   \
           ./progt2.f                                                        \
           ./progtm_module.f                                                 \
           ./radiation_aerosols.f                                            \
           ./radiation_astronomy.f                                           \
           ./radiation_clouds.f                                              \
           ./radiation_gases.f                                               \
           ./radiation_surface.f                                             \
           ./rad_initialize.f                                                \
           ./radlw_datatb.f                                                  \
           ./radlw_main.f                                                    \
           ./radlw_param.f                                                   \
           ./radsw_datatb.f                                                  \
           ./radsw_main.f                                                    \
           ./radsw_param.f                                                   \
           ./rascnvv2.f                                                      \
           ./rayleigh_damp.f                                                 \
           ./rayleigh_damp_mesopause.f                                       \
           ./sascnv.f                                                        \
           ./sascnvn.f                                                       \
           ./set_soilveg.f                                                   \
           ./sfc_diag.f                                                      \
           ./sfc_diff.f                                                      \
           ./sfc_drv.f                                                       \
           ./sfc_land.f                                                      \
           ./sfc_nst.f                                                       \
           ./sfc_ocean.f                                                     \
           ./sfc_sice.f                                                      \
           ./sfcsub.f                                                        \
           ./sflx.f                                                          \
           ./shalcnv.f                                                       \
           ./shalcv.f                                                        \
           ./shalcv_opr.f                                                    \
           ./tracer_const_h.f                                                \
           ./tridi2t3.f

SRCS_f90 = \
           ./calpreciptype.f90                                               \
           ./cs_conv.f90                                                     \
           ./funcphys.f90                                                    \
           ./gcm_shoc.f90                                                    \
           ./gcycle.f90                                                      \
           ./get_prs_fv3.f90                                                 \
           ./module_nst_model.f90                                            \
           ./module_nst_parameters.f90                                       \
           ./module_nst_water_prop.f90                                       \
           ./ozinterp.f90                                                    \
           ./physcons.f90

SRCS_F   = ./gwdc.F                                                          \
           ./gbphys.F                                                        \
           ./machine.F                                                       \
           ./num_parthds.F

SRCS_F90 = \
           ./nuopc_physics.F90

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

./radiation_aerosols.o : ./radiation_aerosols.f
	$(FC) $(FFLAGS) $(OTHER_FFLAGS) -msse3 -c $< -o $@

./radiation_astronomy.o : ./radiation_astronomy.f
	$(FC) $(FFLAGS) $(OTHER_FFLAGS) -xCORE-AVX-I -c $< -o $@

.PHONY: clean
clean:
	@echo "Cleaning gfsphys ... "
	@echo
	$(RM) -f $(LIBRARY) *__genmod.f90 *.o */*.o *.mod *.lst *.i depend

MKDEPENDS = ../mkDepends.pl
include ../conf/make.rules

# do not include 'depend' file if the target contains string 'clean'
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
    -include depend
endif

