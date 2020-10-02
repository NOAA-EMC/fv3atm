if(32BIT)
  message ("Force 64 bits in ipd")
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

if(CCPP)
  list(APPEND _ipd_defs_private CCPP)
endif()

add_library(
    ipd

    IPD_driver.F90
    IPD_typedefs.F90
)
set_target_properties(ipd PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_compile_definitions(ipd PRIVATE "${_ipd_defs_private}")
target_include_directories(ipd PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)

if(CCPP)
  target_include_directories(ipd PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/framework/src
                                         ${CMAKE_BINARY_DIR}/FV3/ccpp/physics)
endif()

target_link_libraries(ipd PRIVATE gfsphysics)
