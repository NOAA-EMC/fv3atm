cmake_minimum_required(VERSION 3.18)

project(MOVING_NEST
  VERSION 1.0.0
  DESCRIPTION  "HRD Moving Nest"
  HOMEPAGE_URL "https://www.gfdl.noaa.gov/gfdl_atmos_cubed_sphere"
  LANGUAGES Fortran)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(GNUInstallDirs)

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|Repro|MinSizeRel|RelWithDebInfo)$")
  message(STATUS "No build type specified.")
  set(CMAKE_BUILD_TYPE "Release"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                     "Debug" "Release" "Repro" "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "Setting build type to '${CMAKE_BUILD_TYPE}'.")

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|GNU)$")
  message(WARNING "Compiler not officially supported: ${CMAKE_Fortran_COMPILER_ID}")
endif()

option(OPENMP      "Enable OpenMP threading"                               OFF)
option(32BIT       "Enable single precision (r4) arithmetic in FV3 dycore" ON)
option(DEBUG       "Enable compiler definition -DDEBUG"                    OFF)
option(MOVING_NEST "Enable compiler definition -DMOVING_NEST"              OFF)
option(MULTI_GASES "Enable compiler definition -DMULTI_GASES"              OFF)
option(USE_GFSL63  "Enable compiler definition -DUSE_GFSL63"               ON)
option(GFS_PHYS    "Enable compiler definition -DGFS_PHYS"                 ON)
option(GFS_TYPES   "Enable compiler definition -DGFS_TYPES"                ON)
option(use_WRTCOMP "Enable compiler definition -Duse_WRTCOMP"              OFF)
option(INTERNAL_FILE_NML     "Enable compiler definition -DINTERNAL_FILE_NML"      ON)
option(ENABLE_QUAD_PRECISION "Enable compiler definition -DENABLE_QUAD_PRECISION"  ON)

find_package(MPI REQUIRED)
if(OPENMP)
  find_package(OpenMP REQUIRED)
endif()
find_package(NetCDF REQUIRED C Fortran)
if(use_WRTCOMP)
  if(NOT ESMF_FOUND)
  find_package(ESMF REQUIRED)
  endif()
endif()

if(32BIT)
  set(kind "R4")
else()
  set(kind "R8")
endif()

if(NOT FMS_FOUND)
  find_package(FMS REQUIRED COMPONENTS ${kind})
  string(TOLOWER ${kind} kind_lower)
  add_library(fms ALIAS FMS::fms_${kind_lower})
endif()


list(APPEND moving_srcs
  bounding_box.F90
  fv_tracker.F90
  fv_moving_nest.F90
  fv_moving_nest_main.F90
  fv_moving_nest_physics.F90
  fv_moving_nest_types.F90
  fv_moving_nest_utils.F90)


list(APPEND mn_defs SPMD
                     MOIST_CAPPA
                     USE_COND)

# Additional (optional) compiler definitions
if(DEBUG)
  list(APPEND mn_defs DEBUG)
endif()

if(USE_GFSL63)
  list(APPEND mn_defs USE_GFSL63)
endif()

if(GFS_PHYS)
  list(APPEND mn_defs GFS_PHYS)
endif()

if(GFS_TYPES)
  list(APPEND mn_defs GFS_TYPES)
endif()

if(use_WRTCOMP)
  list(APPEND mn_defs use_WRTCOMP)
endif()

if(MOVING_NEST)
  list(APPEND mn_defs MOVING_NEST)
endif()

if(MULTI_GASES)
  list(APPEND mn_defs MULTI_GASES)
endif()

if(32BIT)
  list(APPEND mn_defs OVERLOAD_R4
                       OVERLOAD_R8)
endif()

if(INTERNAL_FILE_NML)
  list(APPEND mn_defs INTERNAL_FILE_NML)
endif()

if(ENABLE_QUAD_PRECISION)
  list(APPEND mn_defs ENABLE_QUAD_PRECISION)
endif()

if(OPENMP)
  list(APPEND mn_defs OPENMP)
endif()

# Obtain compiler-specific flags
include(moving_nest_compiler_flags)

# Enable -fPIC compiler flag for all targets
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(moving_nest ${moving_srcs})
add_library(MOVING_NEST::moving_nest ALIAS moving_nest)

add_dependencies(moving_nest fv3)

set_target_properties(moving_nest PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/moving_nest)
target_include_directories(moving_nest PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include/moving_nest>
                                         $<INSTALL_INTERFACE:include/moving_nest>)

target_include_directories(moving_nest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                       ${CMAKE_CURRENT_BINARY_DIR}/../atmos_cubed_sphere/include/fv3
                                       ${CMAKE_CURRENT_SOURCE_DIR}/../atmos_cubed_sphere/model
                                       ${CMAKE_CURRENT_SOURCE_DIR}/../atmos_cubed_sphere/tools)


target_compile_definitions(moving_nest PRIVATE "${mn_defs}")

target_link_libraries(moving_nest PUBLIC fms)

if(GFS_PHYS)
  target_link_libraries(moving_nest PUBLIC fv3ccpp)
endif()

if(use_WRTCOMP)
  target_link_libraries(moving_nest PUBLIC esmf)
endif()

if(OPENMP)
  target_link_libraries(moving_nest PUBLIC OpenMP::OpenMP_Fortran)
endif()

# Install compiled Fortran module files
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})

install(
  TARGETS moving_nest
  EXPORT MOVING_NESTExports
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

### Package config
include(CMakePackageConfigHelpers)
set(CONFIG_INSTALL_DESTINATION lib/cmake/moving_nest)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/MOVING_NESTConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/moving_nest-config.cmake
  INSTALL_DESTINATION ${CONFIG_INSTALL_DESTINATION})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/moving_nest-config.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/moving_nest-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/moving_nest-config-version.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

install(EXPORT MOVING_NESTExports
  NAMESPACE MOVING_NEST::
  FILE moving_nest-targets.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

