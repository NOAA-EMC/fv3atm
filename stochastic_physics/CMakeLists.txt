add_library(stochastic_physics_wrapper stochastic_physics_wrapper.F90)

target_include_directories(stochastic_physics_wrapper PRIVATE ${CMAKE_BINARY_DIR}/stochastic_physics)
target_include_directories(stochastic_physics_wrapper PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/framework/src
                                                              ${CMAKE_BINARY_DIR}/FV3/ccpp/physics)
target_include_directories(stochastic_physics_wrapper PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
if(32BIT)
   target_link_libraries(stochastic_physics_wrapper PUBLIC fms_r4
                                                           stochastic_physics
                                                           gfsphysics
                                                           fv3dycore)
else()
   target_link_libraries(stochastic_physics_wrapper PUBLIC fms_r8
                                                           stochastic_physics
                                                           gfsphysics
                                                           fv3dycore)
endif()

if(OpenMP_Fortran_FOUND)
  target_link_libraries(stochastic_physics_wrapper PUBLIC OpenMP::OpenMP_Fortran)
endif()
