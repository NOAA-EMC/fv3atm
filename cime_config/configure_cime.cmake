message("")
message("Setting configuration for $ENV{CMAKE_Platform}")
message("")
include( $ENV{CASEROOT}/Macros.cmake )

get_filename_component (C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)
get_filename_component (CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message("C       compiler: ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION} (${C_COMPILER_NAME})")
message("CXX     compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} (${CXX_COMPILER_NAME})")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION} (${Fortran_COMPILER_NAME})")
message("")

option(DEBUG   "Enable DEBUG mode" OFF)
option(REPRO   "Enable REPRO mode" OFF)
option(VERBOSE "Enable VERBOSE mode" OFF)
option(32BIT   "Enable 32BIT (single precision arithmetic in dycore)" OFF)
option(OPENMP  "Enable OpenMP threading" ON)
option(AVX2    "Enable AVX2 instruction set" OFF)

option(INLINE_POST "Enable inline post" OFF)

message("1: CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS}")
message("2: CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS}")

if(32BIT)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FC_32BIT}")
    add_definitions(-DOVERLOAD_R4)
    add_definitions(-DOVERLOAD_R8)
else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FC_64BIT}")
    if(NOT REPRO)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-prec-div -no-prec-sqrt")
    endif()
endif()
message("3: CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS}")

if(REPRO)
    add_definitions(-DREPRO)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FC_REPRO}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CC_REPRO}")
endif()
message("4: CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS}")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__IFC ${CFLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${LDFLAGS} ${SLIBS}")
# print build options

if(DEBUG)
    message("DEBUG  is      ENABLED")
else()
    message("DEBUG  is      disabled")
endif()

if(REPRO)
    message("REPRO  is      ENABLED")
else()
    message("REPRO  is      disabled")
endif()

if(32BIT)
    message("32BIT  is      ENABLED")
else()
    message("32BIT  is      disabled")
endif()

if(OPENMP)
    message("OPENMP is      ENABLED")
else()
    message("OPENMP is      disabled")
endif()

if(AVX2)
    message("AVX2 is        ENABLED")
else()
    message("AVX2 is        disabled")
endif()

if(INLINE_POST)
    message("INLINE_POST is ENABLED")
else()
    message("INLINE_POST is disabled")
endif()


set(NEMSIO_INC $ENV{NEMSIO_INC})
set(POST_INC $ENV{POST_INC})
set(NCEP_LIBS $ENV{POST_LIB} $ENV{NEMSIO_LIB} $ENV{G2_LIB4} $ENV{G2TMPL_LIB} $ENV{BACIO_LIB4} $ENV{SP_LIBd} $ENV{W3EMC_LIBd} $ENV{W3NCO_LIBd} $ENV{CRTM_LIB} $ENV{PNG_LIB} $ENV{JASPER_LIB} $ENV{Z_LIB})

set(ESMF_MOD ${ESMF_F90COMPILEPATHS})
set(ESMF_LIBS "${ESMF_F90ESMFLINKRPATHS} ${ESMF_F90ESMFLINKPATHS} ${ESMF_F90ESMFLINKLIBS}")

set(NETCDF_INC_DIR $ENV{NETCDF}/include)
set(NETCDF_LIBDIR $ENV{NETCDF}/lib)
set(NETCDF_LIBS -L$ENV{NETCDF}/lib -lnetcdff -lnetcdf)

message("")
