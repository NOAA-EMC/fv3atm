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
    $(info Build standalone FV3 fv3core ...)
    $(info )
endif

LIBRARY  = libfv3core.a

FFLAGS   += -I../fms -I../fms/include -I../gfsphys

SRCS_f   =

SRCS_f90 =

SRCS_F   =

SRCS_F90 = \
		   ./atmos_cubed_sphere/model/a2b_edge.F90                           \
		   ./atmos_cubed_sphere/model/boundary.F90                           \
		   ./atmos_cubed_sphere/model/dyn_core.F90                           \
		   ./atmos_cubed_sphere/model/fv_arrays.F90                          \
		   ./atmos_cubed_sphere/model/fv_control.F90                         \
		   ./atmos_cubed_sphere/model/fv_current_grid.F90                    \
		   ./atmos_cubed_sphere/model/fv_dynamics.F90                        \
		   ./atmos_cubed_sphere/model/fv_fill.F90                            \
		   ./atmos_cubed_sphere/model/fv_grid_utils.F90                      \
		   ./atmos_cubed_sphere/model/fv_mapz.F90                            \
		   ./atmos_cubed_sphere/model/fv_nesting.F90                         \
		   ./atmos_cubed_sphere/model/fv_sg.F90                              \
		   ./atmos_cubed_sphere/model/fv_tracer2d.F90                        \
		   ./atmos_cubed_sphere/model/fv_update_phys.F90                     \
		   ./atmos_cubed_sphere/model/sw_core.F90                            \
		   ./atmos_cubed_sphere/model/tp_core.F90                            \
		   ./atmos_cubed_sphere/model_nh/nh_core.F90                         \
		   ./atmos_cubed_sphere/model_nh/nh_utils.F90                        \
		   ./atmos_cubed_sphere/tools/external_ic.F90                        \
		   ./atmos_cubed_sphere/tools/external_sst.F90                       \
		   ./atmos_cubed_sphere/tools/fv_diagnostics.F90                     \
		   ./atmos_cubed_sphere/tools/fv_eta.F90                             \
		   ./atmos_cubed_sphere/tools/fv_grid_tools.F90                      \
		   ./atmos_cubed_sphere/tools/fv_io.F90                              \
		   ./atmos_cubed_sphere/tools/fv_mp_mod.F90                          \
		   ./atmos_cubed_sphere/tools/fv_nggps_diag.F90                      \
		   ./atmos_cubed_sphere/tools/fv_nudge.F90                           \
		   ./atmos_cubed_sphere/tools/fv_restart.F90                         \
		   ./atmos_cubed_sphere/tools/fv_surf_map.F90                        \
		   ./atmos_cubed_sphere/tools/fv_timing.F90                          \
		   ./atmos_cubed_sphere/tools/init_hydro.F90                         \
		   ./atmos_cubed_sphere/tools/sim_nc_mod.F90                         \
		   ./atmos_cubed_sphere/tools/sorted_index.F90                       \
		   ./atmos_cubed_sphere/tools/test_cases.F90                         \
		   ./atmos_cubed_sphere/driver/coupled/atmosphere.F90                \
		   ./atmos_cubed_sphere/lin_cloud_microphys/lin_cloud_microphys.F90  \
		   ./atmos_drivers/coupled/atmos_model.F90                           \
           ./gfsphys_driver/gfs_physics_driver.F90

SRCS_c   = \
		   ./avec_timer/avec_timer_lib.c

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

./atmos_cubed_sphere/model_nh/nh_utils.o : ./atmos_cubed_sphere/model_nh/nh_utils.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHER_FFLAGS) $(FAST) -c $< -o $@

./atmos_cubed_sphere/model/fv_mapz.o : ./atmos_cubed_sphere/model/fv_mapz.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHER_FFLAGS) $(FAST) -c $< -o $@

.PHONY: clean
clean:
	@echo "Cleaning fv3core ... "
	@echo
	$(RM) -f $(LIBRARY) *__genmod.f90 *.o */*.o */*/*.o *.mod *.lst *.i depend

MKDEPENDS = ../mkDepends.pl
include ../conf/make.rules

# do not include 'depend' file if the target contains string 'clean'
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
    -include depend
endif

