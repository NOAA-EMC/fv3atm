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
    $(info Build standalone FV3 namphysics ...)
    $(info )
endif

LIBRARY  = libnamphys.a

FFLAGS   += -I$(FMS_DIR) -I../cpl $(ESMF_INC)

CPPDEFS = -DNEW_TAUCTMAX -DSMALL_PE -DNEMS_GSM -DINTERNAL_FILE_NML

SRCS_f   =  \
           ./physics/grrad_nmmb.f                                                    \
           ./physics/gscond.f                                                        \
           ./physics/iounitdef.f                                                     \
           ./physics/mersenne_twister.f                                              \
           ./physics/mfdeepcnv.f                                                     \
           ./physics/mfpbl.f                                                         \
           ./physics/mfshalcnv.f                                                     \
           ./physics/moninedmf.f                                                     \
           ./physics/moninq.f                                                        \
           ./physics/namelist_soilveg.f                                              \
           ./physics/ozne_def.f                                                      \
           ./physics/physparam.f                                                     \
           ./physics/precpd_nmmb.f                                                   \
           ./physics/rad_initialize_nmmb.f                                           \
           ./physics/radiation_aerosols_nmmb.f                                       \
           ./physics/radiation_astronomy_nmmb.f                                      \
           ./physics/radiation_clouds_nmmb.f                                         \
           ./physics/radiation_gases_nmmb.f                                          \
           ./physics/radiation_surface_nmmb.f                                        \
           ./physics/radlw_datatb.f                                                  \
           ./physics/radlw_main_nmmb.f                                               \
           ./physics/radlw_param.f                                                   \
           ./physics/radsw_datatb.f                                                  \
           ./physics/radsw_main_nmmb.f                                               \
           ./physics/radsw_param.f                                                   \
           ./physics/sascnvn.f                                                       \
           ./physics/h2o_def.f                                                       \
           ./physics/shalcnv.f

SRCS_F   = ./physics/machine.F

SRCS_f90 = \
           ./physics/funcphys.f90                                                    \
           ./physics/physcons.f90

SRCS_F90 = \
           ./physics/module_BL_GFSPBL.F90                                            \
           ./physics/module_BL_GFSPBLEDMFHUR.F90                                     \
           ./physics/module_BL_GFSPBLHUR.F90                                         \
           ./physics/module_BL_MYJPBL.F90                                            \
           ./physics/module_CONSTANTS.F90                                            \
           ./physics/module_CONTROL.F90                                              \
           ./physics/module_CONVECTION.F90                                           \
           ./physics/module_CU_BMJ.F90                                               \
           ./physics/module_CU_SAS.F90                                               \
           ./physics/module_CU_SASHUR.F90                                            \
           ./physics/module_CU_SCALE.F90                                             \
           ./physics/module_DIAGNOSE.F90                                             \
           ./physics/module_GWD.F90                                                  \
           ./physics/module_KINDS.F90                                                \
           ./physics/module_LS_LISS.F90                                              \
           ./physics/module_LS_NOAHLSM.F90                                           \
           ./physics/module_MICROPHYSICS.F90                                         \
           ./physics/module_MP_ETANEW.F90                                            \
           ./physics/module_MP_FER_HIRES.F90                                         \
           ./physics/module_MP_GFS.F90                                               \
           ./physics/module_MP_WSM6.F90                                              \
           ./physics/module_RADIATION.F90                                            \
           ./physics/module_RA_GFDL.F90                                              \
           ./physics/module_RA_RRTM.F90                                              \
           ./physics/module_SF_GFDL.F90                                              \
           ./physics/module_SF_JSFC.F90                                              \
           ./physics/module_SF_URBAN.F90                                             \
           ./physics/module_TURBULENCE.F90                                           \
           ./physics/module_mp_radar.F90                                             \
           ./physics/module_mp_thompson.F90                                          \
           ./physics/GFDL_parse_tracers.F90                                          \
           ./physics/gfdl_cloud_microphys.F90                                        \
           ./NAM_layer/module_SOLVER_GRID_COMP.F90                                   \
	   ./NAM_layer/NAM_abstraction_layer.F90 				     \
	   ./NAM_layer/NAM_diagnostics.F90       				     \
	   ./NAM_layer/NAM_driver.F90            				     \
           ./NAM_layer/NAM_radiation_driver.F90                                      \
	   ./NAM_layer/NAM_restart.F90           				     \
	   ./NAM_layer/NAM_typedefs.F90

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

./physics/grrad_nmmb.o : ./physics/grrad_nmmb.f
	$(FC) $(FFLAGS) -free -r8 -c $< -o $@
./physics/radiation_astronomy_nmmb.o : ./physics/radiation_astronomy_nmmb.f
	$(FC) $(FFLAGS) -free -r8 -c $< -o $@
./physics/radiation_clouds_nmmb.o : ./physics/radiation_clouds_nmmb.f
	$(FC) $(FFLAGS) -free -r8 -c $< -o $@
./physics/radsw_main_nmmb.o : ./physics/radsw_main_nmmb.f
	$(FC) $(FFLAGS) -free -r8 -c $< -o $@

./NAM_layer/NAM_diagnostics.o : ./NAM_layer/NAM_diagnostics.F90
	$(FC) $(FFLAGS) $(OTHER_FFLAGS) -O0 -c $< -o $@

.PHONY: clean
clean:
	@echo "Cleaning namphysics  ... "
	@echo
	$(RM) -f $(LIBRARY) *__genmod.f90 *.o */*.o *.mod *.i90 *.lst *.i depend

MKDEPENDS = ../mkDepends.pl
include ../conf/make.rules

include ./depend

# do not include 'depend' file if the target contains string 'clean'
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
   -include depend
endif

