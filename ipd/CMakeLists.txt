
if(32BIT)
remove_definitions(-DOVERLOAD_R4)
remove_definitions(-DOVERLOAD_R8)
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

add_library(
    ipd

    IPD_driver.F90
    IPD_typedefs.F90
)
target_link_libraries(ipd PRIVATE gfsphysics)

if(CCPP)
  target_include_directories(ipd PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/framework/src)
  target_include_directories(ipd PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/physics)
endif()

set_target_properties(ipd PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(ipd PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
