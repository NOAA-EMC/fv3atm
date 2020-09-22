add_library(
    stochastic_physics_wrapper

    ./stochastic_physics_wrapper.F90
)

target_include_directories(stochastic_physics_wrapper PRIVATE ${CMAKE_BINARY_DIR}/stochastic_physics)
if(CCPP)
  target_include_directories(stochastic_physics_wrapper PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/framework/src)
  target_include_directories(stochastic_physics_wrapper PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/physics)
endif()
# DH* is there a better way to do this?
target_include_directories(stochastic_physics_wrapper PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../mod>)
# *DH
target_link_libraries(stochastic_physics_wrapper fms)
target_link_libraries(stochastic_physics_wrapper stochastic_physics)
target_link_libraries(stochastic_physics_wrapper gfsphysics)
target_link_libraries(stochastic_physics_wrapper fv3dycore)

