
if(NOT DYN32)
remove_definitions(-DOVERLOAD_R8)
remove_definitions(-DOVERLOAD_R4)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    message ("Force 64 bits in CCPP_layer")
    string (REPLACE "-i4 -real-size 32" "-i4 -real-size 64 -no-prec-div -no-prec-sqrt" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
endif()

add_library(
    ccppdriver

    CCPP_driver.F90
)

target_link_libraries(ccppdriver ${PHYS}physics)
target_link_libraries(ccppdriver ccpp)
target_link_libraries(ccppdriver ccppphys)

target_include_directories(ccppdriver PRIVATE ${CMAKE_BINARY_DIR}/ccpp/framework/src)

set_target_properties(ccppdriver PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(ccppdriver PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/mod)
