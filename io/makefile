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
    $(info Build standalone FV3 io ...)
    $(info )
endif
$(info $$ESMF_INC is [${ESMF_INC}])

LIBRARY  = libfv3io.a

FFLAGS   += -I$(FMS_DIR) -I../gfsphysics -I../ipd -I../namphysics

ifneq (,$(findstring NO_INLINE_POST,$(CPPDEFS)))
POST_SRC = \
		   ./post_gfs_stub.F90
else
POST_SRC = \
		   ./post_gfs.F90                      \
		   ./post_nems_routines.F90
endif

SRCS_f   =

SRCS_f90 = 

SRCS_F   =

SRCS_F90 = \
		   ./ffsync.F90                        \
		   ./FV3GFS_io.F90                     \
		   $(POST_SRC)                         \
		   ./module_write_nemsio.F90           \
		   ./module_write_netcdf.F90           \
		   ./module_fv3_io_def.F90             \
		   ./module_write_internal_state.F90   \
		   ./module_wrt_grid_comp.F90

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

FV3GFS_io.o: FV3GFS_io.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) $(ESMF_INC) -c FV3GFS_io.F90
post_nems_routines.o: post_nems_routines.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -I$(POST_INC) -c post_nems_routines.F90
module_write_nemsio.o: module_write_nemsio.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) $(ESMF_INC) $(NEMSIOINC) -c module_write_nemsio.F90
module_write_netcdf.o: module_write_netcdf.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) $(ESMF_INC) $(NEMSIOINC) -c module_write_netcdf.F90
module_write_internal_state.o: module_write_internal_state.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) $(ESMF_INC) -c module_write_internal_state.F90
module_wrt_grid_comp.o: module_wrt_grid_comp.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) $(ESMF_INC) -c module_wrt_grid_comp.F90
post_gfs.o: post_gfs.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) $(ESMF_INC) -I$(POST_INC) -c post_gfs.F90
post_gfs_stub.o: post_gfs_stub.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) $(ESMF_INC) -c post_gfs_stub.F90


.PHONY: clean
clean:
	@echo "Cleaning io ... "
	@echo
	$(RM) -f $(LIBRARY) *__genmod.f90 *.o *.mod *.i90 *.lst *.i depend

MKDEPENDS = ../mkDepends.pl
include ../conf/make.rules

# do not include 'depend' file if the target contains string 'clean'
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
    -include depend
endif

