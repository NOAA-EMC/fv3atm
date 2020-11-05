
if(NOT DYN32)
remove_definitions(-DOVERLOAD_R8)
remove_definitions(-DOVERLOAD_R4)
endif()

message ("Force 64 bits in CCPP_layer")
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    string (REPLACE "-i4 -real-size 32" "-i4 -real-size 64 -no-prec-div -no-prec-sqrt" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
endif()

add_library(
    ccppdriver

    CCPP_driver.F90
)

target_link_libraries(ccppdriver gfsphysics)
target_link_libraries(ccppdriver ccpp)
target_link_libraries(ccppdriver ccppphys)

target_include_directories(ccppdriver PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/framework/src)

set_target_properties(ccppdriver PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(ccppdriver PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
