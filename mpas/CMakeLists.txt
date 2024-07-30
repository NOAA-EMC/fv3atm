cmake_minimum_required(VERSION 3.19)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i4 -r8")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

# Build product output locations
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Set default build type to RelWithDebInfo
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting default build type to Release.  Specify CMAKE_BUILD_TYPE to override.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "CMake Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Detect MPAS git version
if(NOT MPAS_GIT_VERSION)
    find_package(Git QUIET)
    if(GIT_FOUND)
        execute_process(COMMAND ${GIT_EXECUTABLE} describe --dirty
                        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
                        OUTPUT_VARIABLE _mpas_git_version
                        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    else()
        set(_mpas_git_version "Unknown")
    endif()
    set(MPAS_GIT_VERSION ${_mpas_git_version} CACHE STRING "MPAS-Model git version")
endif()
 
# Dependencies
find_package(OpenMP COMPONENTS Fortran)
find_package(MPI REQUIRED COMPONENTS Fortran)
find_package(NetCDF REQUIRED COMPONENTS Fortran C)
find_package(PnetCDF REQUIRED COMPONENTS Fortran)
find_package(PIO REQUIRED COMPONENTS Fortran C)

# Find C pre-processor
if(CMAKE_C_COMPILER_ID MATCHES GNU)
  find_program(CPP_EXECUTABLE NAMES cpp REQUIRED)
  set(CPP_EXTRA_FLAGS -traditional)
elseif(CMAKE_C_COMPILER_ID MATCHES "(Apple)?Clang" )
  find_program(CPP_EXECUTABLE NAMES cpp REQUIRED)
else()
  message(STATUS "Unknown compiler: ${CMAKE_C_COMPILER_ID}")
  set(CPP_EXECUTABLE ${CMAKE_C_COMPILER})
endif()

# Fortran module output directory for build interface
set(MPAS_MODULE_DIR ${PROJECT_NAME}/module/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION})
# Install Fortran module directory
install(DIRECTORY ${CMAKE_BINARY_DIR}/${MPAS_MODULE_DIR}/ DESTINATION ${CMAKE_INSTALL_LIBDIR}/${MPAS_MODULE_DIR}/)

# Registry file.
set(REGISTRY_FILE MPAS-Model/src/core_atmosphere/Registry.xml)
set(REGISTRY_FILE ${CMAKE_CURRENT_SOURCE_DIR}/MPAS-Model/src/core_atmosphere/Registry.xml)

#
# Build libraries...
#

# MPAS Utilities (Externals)
add_subdirectory(MPAS-Model/src/external/ezxml)

# ESMF libraries.
if(NOT ESMF_FOUND)
  find_package(ESMF REQUIRED)
endif()
add_definitions(-DMPAS_EXTERNAL_ESMF_LIB -DMPAS_NO_ESMF_INIT)
add_library(${PROJECT_NAME}::external::esmf ALIAS esmf)

# MPAS Namelist
add_subdirectory(MPAS-Model/src/tools/input_gen) # Targets: namelist_gen, streams_gen

# MPAS Registry
add_subdirectory(MPAS-Model/src/tools/registry) # Targets: mpas_parse_<core_name>

# MPAS framework
add_subdirectory(MPAS-Model/src/framework) # Target: MPAS::framework
add_subdirectory(MPAS-Model/src/operators) # Target: MPAS::operators

# MPAS-Atmosphere

# MPAS-Model/src/core_atmosphere/diagnostics/
set(MPAS_DYCORE_DIAGNOSTIC_SOURCES
  mpas_atm_diagnostic_template.F
  mpas_atm_diagnostics_manager.F
  mpas_atm_diagnostics_utils.F
  mpas_cloud_diagnostics.F
  mpas_convective_diagnostics.F
  mpas_isobaric_diagnostics.F
  mpas_pv_diagnostics.F
  mpas_soundings.F
  )
list(TRANSFORM MPAS_DYCORE_DIAGNOSTIC_SOURCES PREPEND MPAS-Model/src/core_atmosphere/diagnostics/)

# MPAS-Model/src/core_atmosphere/dynamics/
set(MPAS_DYCORE_DYNAMICS_SOURCES
  mpas_atm_boundaries.F
  mpas_atm_iau.F
  mpas_atm_time_integration.F
  )
list(TRANSFORM MPAS_DYCORE_DYNAMICS_SOURCES PREPEND MPAS-Model/src/core_atmosphere/dynamics/)

# MPAS-Model/src/core_atmosphere/core_atosphere
set(MPAS_DYCORE_SOURCES
  mpas_atm_dimensions.F
  mpas_atm_threading.F
  mpas_atm_core.F
  mpas_atm_core_interface.F
  mpas_atm_halos.F
  )
list(TRANSFORM MPAS_DYCORE_SOURCES PREPEND MPAS-Model/src/core_atmosphere/)

# Generated includes
set(MPAS_DYCORE_INCLUDES
  block_dimension_routines.inc
  core_variables.inc
  define_packages.inc
  domain_variables.inc
  namelist_call.inc
  namelist_defines.inc
  setup_immutable_streams.inc
  structs_and_variables.inc
  )
      
# Create library
add_library(mpas_dycore ${MPAS_DYCORE_SOURCES}
                        ${MPAS_DYCORE_DIAGNOSTIC_SOURCES}
			${MPAS_DYCORE_DYNAMICS_SOURCES})
		      
# Compilation definitions
set(MPAS_DYCORE_COMPILE_DEFINITIONS
  mpas=1
  MPAS_NATIVE_TIMERS
  )

# Link library
target_compile_definitions(mpas_dycore PRIVATE ${MPAS_DYCORE_COMPILE_DEFINITIONS})
set_MPAS_DEBUG_flag(mpas_dycore)
set(CORE_MODULE_DIR ${MPAS_MODULE_DIR}/mpas_dycore)
set_target_properties(mpas_dycore PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${CORE_MODULE_DIR})
target_include_directories(mpas_dycore INTERFACE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${CORE_MODULE_DIR}>
  $<INSTALL_INTERFACE:${CORE_MODULE_DIR}>)
#
target_link_libraries(mpas_dycore PUBLIC framework operators)


#Per-core generated output and tables directory location
set(CORE_DATADIR ${CMAKE_BINARY_DIR}/${PROJECT_NAME}/mpas_dycore)
file(MAKE_DIRECTORY ${CORE_DATADIR})
option(DO_PHYSICS "Use built-in physics schemes." FALSE)

# Process registry.
add_custom_command(OUTPUT Registry_processed.xml
  COMMAND ${CPP_EXECUTABLE} -E -P ${CPP_EXTRA_FLAGS} ${REGISTRY_FILE} > Registry_processed.xml
  COMMENT "CORE atmosphere: Pre-Process Registry"
  DEPENDS ${REGISTRY_FILE})
add_custom_command(OUTPUT ${MPAS_DYCORE_INCLUDES}
  COMMAND mpas_parse_atmosphere Registry_processed.xml
  COMMENT "CORE atmosphere: Parse Registry"
  DEPENDS mpas_parse_atmosphere Registry_processed.xml)
add_custom_command(OUTPUT namelist.atmosphere
  WORKING_DIRECTORY ${CORE_DATADIR}
  COMMAND mpas_namelist_gen ${CMAKE_CURRENT_BINARY_DIR}/Registry_processed.xml namelist.atmosphere in_defaults=true
  COMMENT "CORE atmosphere: Generate Namelist"
  DEPENDS mpas_namelist_gen Registry_processed.xml)
add_custom_command(OUTPUT streams.atmosphere
  WORKING_DIRECTORY ${CORE_DATADIR}
  COMMAND mpas_streams_gen ${CMAKE_CURRENT_BINARY_DIR}/Registry_processed.xml streams.atmosphere stream_list.atmosphere. listed
  COMMENT "CORE atmosphere: Generate Streams"
  DEPENDS mpas_streams_gen Registry_processed.xml)
add_custom_target(gen_atmosphere DEPENDS ${MPAS_DYCORE_INCLUDES} namelist.atmosphere streams.atmosphere)
add_dependencies(mpas_dycore gen_atmosphere)
