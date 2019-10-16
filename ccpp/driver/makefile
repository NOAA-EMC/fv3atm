SHELL = /bin/sh

inside_nems := $(wildcard ../../../conf/configure.nems)
ifneq ($(strip $(inside_nems)),)
    include ../../../conf/configure.nems
else
    exist_configure_fv3 := $(wildcard ../../conf/configure.fv3)
    ifneq ($(strip $(exist_configure_fv3)),)
        include ../../conf/configure.fv3
    else
        $(error "../../conf/configure.fv3 file is missing. Run ./configure")
    endif
    $(info )
    $(info Build CCPP layer ...)
    $(info )
endif

CCPP_DRIVER = CCPP_driver.F90

LIBRARY  = libccppdriver.a

# Needed for ccpp_data.mod, fv_arrays_mod.mod, ...
FFLAGS   += -I$(FMS_DIR) -I../../gfsphysics -I../../atmos_cubed_sphere

#CPPDEFS += -DNEW_TAUCTMAX -DSMALL_PE -DNEMS_GSM -DINTERNAL_FILE_NML

# Set flags for 32-bit dynamics build
ifeq ($(DYN32),Y)
CPPDEFS += -DOVERLOAD_R4
endif

SRCS_F90 = \
	   CCPP_driver.F90

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

# Do preprocessing of the CCPP driver in two steps to be
# able to look at the actual .f90 file that gets compiled
./CCPP_driver.o: ./CCPP_driver.F90
	$(CPP) $(CPPDEFS) $(CPPFLAGS) $< > $*.tmp.f90
	$(FC) $(FFLAGS) $(OTHER_FFLAGS) -c $*.tmp.f90 -o $@

.PHONY: clean
clean:
	@echo "Cleaning CCPP_layer  ... "
	@echo
	$(RM) -f $(LIBRARY) *__genmod.f90 *.tmp.f90 *.o */*.o *.mod *.i90 *.lst *.i depend

MKDEPENDS = ../../mkDepends.pl
include ../../conf/make.rules

include ./depend

# do not include 'depend' file if the target contains string 'clean'
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
   -include depend
endif

