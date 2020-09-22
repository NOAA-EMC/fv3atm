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

LIBRARY  = libstochastic_physics_wrapper.a

FFLAGS   += -I$(FMS_DIR) -I ../../stochastic_physics -I../ccpp/physics -I../ccpp/build/physics -I../atmos_cubed_sphere

SRCS_f   =

SRCS_f90 = 

SRCS_F   =

SRCS_F90 = stochastic_physics_wrapper.F90

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

.PHONY: clean
clean:
	@echo "Cleaning io ... "
	@echo
	$(RM) -f $(LIBRARY) *__genmod.f90 *.o *.mod *.i90 *.lst *.i

include ../conf/make.rules
