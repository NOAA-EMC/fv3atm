
if(NOT DYN32)
  remove_definitions(-DOVERLOAD_R8)
  remove_definitions(-DOVERLOAD_R4)
endif()

message ("Force 64 bits in ccpp/data")
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  if(REPRO)
    string (REPLACE "-i4 -real-size 32" "-i4 -real-size 64" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  else()
    string (REPLACE "-i4 -real-size 32" "-i4 -real-size 64 -no-prec-div -no-prec-sqrt" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
endif()

list(APPEND _ccppdata_defs_private NEMS_GSM
                                    MOIST_CAPPA
                                    USE_COND
                                    INTERNAL_FILE_NML)

if(MULTI_GASES)
  list(APPEND _ccppdata_defs_private MULTI_GASES)
endif()

add_library(
    ccppdata

    CCPP_typedefs.F90
    GFS_typedefs.F90
    CCPP_data.F90
)

target_link_libraries(ccppdata PUBLIC ccpp)
target_link_libraries(ccppdata PUBLIC ccppphys)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(ccppdata PUBLIC OpenMP::OpenMP_Fortran)
endif()

target_include_directories(ccppdata PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/framework/src
                                            ${CMAKE_BINARY_DIR}/FV3/ccpp/physics)

set_target_properties(ccppdata PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_compile_definitions(ccppdata PRIVATE "${_ccppdata_defs_private}")
target_include_directories(ccppdata PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
