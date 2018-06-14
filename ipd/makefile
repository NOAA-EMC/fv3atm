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

LIBRARY  = libipd.a

FFLAGS   += -I../fms -I../fms/include -I../gfsphysics -I../namphysics

CPPDEFS = -DNEW_TAUCTMAX -DSMALL_PE -DNEMS_GSM -DINTERNAL_FILE_NML

SRCS_F90 = \
	   ./IPD_driver.F90						     \
	   ./IPD_typedefs.F90

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

ifneq (,$(findstring axSSE4.2,AVX,CORE-AVX2,$(FFLAGS)))
./physics/radiation_aerosols.o : ./physics/radiation_aerosols.f
	$(FC) $(FFLAGS) $(OTHER_FFLAGS) -axSSE4.2,AVX,CORE-AVX-I -c $< -o $@
else
./physics/radiation_aerosols.o : ./physics/radiation_aerosols.f
	$(FC) $(FFLAGS) $(OTHER_FFLAGS) -xCORE-AVX-I -c $< -o $@
endif

./GFS_layer/GFS_diagnostics.o : ./GFS_layer/GFS_diagnostics.F90
	$(FC) $(FFLAGS) $(OTHER_FFLAGS) -O0 -c $< -o $@

.PHONY: clean
clean:
	@echo "Cleaning ipd  ... "
	@echo
	$(RM) -f $(LIBRARY) *__genmod.f90 *.o */*.o *.mod *.lst *.i depend

MKDEPENDS = ../mkDepends.pl
include ../conf/make.rules

include ./depend

# do not include 'depend' file if the target contains string 'clean'
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
   -include depend
endif

