
add_library(
    fv3cpl

    module_cplfields.F90
    module_cap_cpl.F90
)

target_include_directories(fv3cpl PRIVATE ${ESMF_MOD})

set_target_properties(fv3cpl PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(fv3cpl PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/mod)
