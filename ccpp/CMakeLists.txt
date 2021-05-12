cmake_minimum_required(VERSION 3.0)

project(CCPP-FV3
        LANGUAGES C CXX Fortran)
set(PROJECT "CCPP-FV3")

# Attempt to add link library "NetCDF::NetCDF_Fortran" to target "ccppphys"
# which is not built in this directory.
cmake_policy(SET CMP0079 NEW)

#------------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "Coverage")
endif()

#------------------------------------------------------------------------------
# CMake Modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/framework/cmake")

#------------------------------------------------------------------------------
# Call to CCPP code generator
if(DEBUG)
  set(_ccpp_debug_arg "--debug")
endif()
if(DEFINED CCPP_SUITES)
  set(_ccpp_suites_arg "--suites=${CCPP_SUITES}")
  message("Calling CCPP code generator (ccpp_prebuild.py) for suites ${_ccpp_suites_arg} ...")
else()
  message("Calling CCPP code generator (ccpp_prebuild.py) for all available suites ...")
endif()
execute_process(COMMAND ${Python_EXECUTABLE}
                        "framework/scripts/ccpp_prebuild.py"
                        "--config=config/ccpp_prebuild_config.py"
                        "--builddir=${CMAKE_CURRENT_BINARY_DIR}" ${_ccpp_suites_arg} ${_ccpp_debug_arg}
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/ccpp_prebuild.out
                ERROR_FILE  ${CMAKE_CURRENT_BINARY_DIR}/ccpp_prebuild.err
                RESULT_VARIABLE RC)
# Check return code from ccpp_prebuild.py
if(NOT RC EQUAL 0)
  message(FATAL_ERROR "An error occured while running ccpp_prebuild.py, check ${CMAKE_CURRENT_BINARY_DIR}/ccpp_prebuild.{out,err}")
endif()

#------------------------------------------------------------------------------
# Add host-model specific preprocessor flag (needed for some physics schemes)
add_definitions(-DFV3)

#------------------------------------------------------------------------------
# Set MPI flags for C/C++/Fortran preprocessor
if(MPI)
  add_definitions(-DMPI)
endif()

#------------------------------------------------------------------------------
# Set flag for 32bit dynamics build
if(32BIT)
  set(DYN32 ON CACHE BOOL "Enable support for 32bit fast physics in CCPP")
  message(STATUS "Dynamics compiled with 32-bit option, adjust fv_sat_adj types")
  add_definitions(-DOVERLOAD_R4)
  message ("Force 64 bits in ccpp")
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -real-size 64")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
  endif()
else()
  set(DYN32 OFF CACHE BOOL "Disable support for 32bit fast physics in CCPP")
  remove_definitions(-DOVERLOAD_R8)
  remove_definitions(-DOVERLOAD_R4)
endif()

#------------------------------------------------------------------------------
# Add model-specific flags for C/C++/Fortran preprocessor
add_definitions(-DMOIST_CAPPA -DUSE_COND -DNEMS_GSM)
add_definitions(-DINTERNAL_FILE_NML)

if(MULTI_GASES)
  add_definitions(-DMULTI_GASES)
endif()

#------------------------------------------------------------------------------
# Build CCPP framework and physics

add_subdirectory(framework)

add_subdirectory(physics)
add_dependencies(ccppphys ccpp)
target_link_libraries(ccppphys PUBLIC NetCDF::NetCDF_Fortran)
# This should not be necessary once framework and physics targets define BUILD_INTERFACE
target_include_directories(ccppphys PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/framework/src)

#------------------------------------------------------------------------------
# Build fv3ccpp

# Can we move data/*.F90 and driver/*.F90 to this directory ???
add_library(
    fv3ccpp

    data/CCPP_typedefs.F90
    data/GFS_typedefs.F90
    data/CCPP_data.F90

    driver/GFS_diagnostics.F90
    driver/GFS_restart.F90
    driver/GFS_init.F90
    driver/CCPP_driver.F90

    ${CMAKE_CURRENT_BINARY_DIR}/physics/ccpp_static_api.F90
)

# Compile GFS_diagnostics.F90 without optimization, this leads to out of memory errors on wcoss_dell_p3
set_property(SOURCE driver/GFS_diagnostics.F90 APPEND_STRING PROPERTY COMPILE_FLAGS "-O0")


target_link_libraries(fv3ccpp PUBLIC ccpp)
target_link_libraries(fv3ccpp PUBLIC ccppphys)

if(OPENMP)
  target_link_libraries(fv3ccpp PUBLIC OpenMP::OpenMP_Fortran)
endif()

# This should not be necessary once framework and physics targets define BUILD_INTERFACE
target_include_directories(fv3ccpp PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/framework/src
                                           ${CMAKE_CURRENT_BINARY_DIR}/physics)

set_target_properties(fv3ccpp PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(fv3ccpp PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
