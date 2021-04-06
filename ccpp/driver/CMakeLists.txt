
if(NOT DYN32)
  remove_definitions(-DOVERLOAD_R8)
  remove_definitions(-DOVERLOAD_R4)
endif()

message ("Force 64 bits in ccpp/driver")
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  if(REPRO)
    string (REPLACE "-i4 -real-size 32" "-i4 -real-size 64" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  else()
    string (REPLACE "-i4 -real-size 32" "-i4 -real-size 64 -no-prec-div -no-prec-sqrt" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
endif()

list(APPEND _ccppdriver_defs_private NEMS_GSM
                                    MOIST_CAPPA
                                    USE_COND
                                    INTERNAL_FILE_NML)

if(MULTI_GASES)
  list(APPEND _ccppdriver_defs_private MULTI_GASES)
endif()

add_library(
    ccppdriver

    GFS_diagnostics.F90
    GFS_restart.F90
    GFS_init.F90

    ${CMAKE_BINARY_DIR}/FV3/ccpp/physics/ccpp_static_api.F90
    CCPP_driver.F90
)

# Compile GFS_diagnostics.F90 without optimization, this leads to out of memory errors on wcoss_dell_p3
set_property(SOURCE GFS_diagnostics.F90 APPEND_STRING PROPERTY COMPILE_FLAGS "-O0")

target_link_libraries(ccppdriver PUBLIC ccpp)
target_link_libraries(ccppdriver PUBLIC ccppphys)
target_link_libraries(ccppdriver PUBLIC ccppdata)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(ccppdriver PUBLIC OpenMP::OpenMP_Fortran)
endif()

target_include_directories(ccppdriver PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/framework/src
                                              ${CMAKE_BINARY_DIR}/FV3/ccpp/physics)

set_target_properties(ccppdriver PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_compile_definitions(ccppdriver PRIVATE "${_ccppdriver_defs_private}")
target_include_directories(ccppdriver PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
