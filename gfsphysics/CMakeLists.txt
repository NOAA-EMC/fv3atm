if(32BIT)
  message ("Force 64 bits in gfsphysics")
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    if(REPRO)
      string (REPLACE "-i4 -real-size 32" "-i4 -real-size 64" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
    else()
      string (REPLACE "-i4 -real-size 32" "-i4 -real-size 64 -no-prec-div -no-prec-sqrt" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
    endif()
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
  endif()
endif()

set(CCPP_SOURCES
    physics/noahmp_tables.f90
)

list(APPEND _gfsphysics_defs_private NEMS_GSM
                                     MOIST_CAPPA
                                     USE_COND
                                     INTERNAL_FILE_NML)

if(MULTI_GASES)
  list(APPEND _gfsphysics_defs_private MULTI_GASES)
endif()

list(APPEND _gfsphysics_srcs ${CCPP_SOURCES})
list(APPEND _gfsphysics_defs_private CCPP)
if(DYN32)
  list(APPEND _gfsphysics_defs_private OVERLOAD_R4)
endif()

add_library(gfsphysics ${_gfsphysics_srcs})

target_include_directories(gfsphysics PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/framework/src
                                              ${CMAKE_BINARY_DIR}/FV3/ccpp/physics)
target_link_libraries(gfsphysics PRIVATE ccppphys ccpp)

set_property(SOURCE GFS_layer/GFS_diagnostics.F90 APPEND_STRING PROPERTY COMPILE_FLAGS "-O0")

set_target_properties(gfsphysics PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_compile_definitions(gfsphysics PRIVATE "${_gfsphysics_defs_private}")
target_include_directories(gfsphysics PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
target_link_libraries(gfsphysics PRIVATE NetCDF::NetCDF_Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(gfsphysics PUBLIC OpenMP::OpenMP_Fortran)
endif()
