add_library(stochastic_physics_wrapper stochastic_physics_wrapper.F90)

target_include_directories(stochastic_physics_wrapper PRIVATE ${CMAKE_BINARY_DIR}/stochastic_physics)
target_include_directories(stochastic_physics_wrapper PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/framework/src
                                                              ${CMAKE_BINARY_DIR}/FV3/ccpp/physics)
target_include_directories(stochastic_physics_wrapper PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
target_link_libraries(stochastic_physics_wrapper PUBLIC fms
                                                        stochastic_physics
                                                        fv3dycore)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(stochastic_physics_wrapper PUBLIC OpenMP::OpenMP_Fortran)
endif()
