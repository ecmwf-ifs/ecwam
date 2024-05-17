# (C) Copyright 2022- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

cmake_minimum_required( VERSION 3.24 FATAL_ERROR )
find_package( ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )

# If the directory name matches ecwam_dp or ecwam_sp, that will be chosen as PROJECT_NAME instead
get_filename_component( ecwam-source_dir ${CMAKE_CURRENT_SOURCE_DIR} NAME )
if( ecwam-source_dir MATCHES "ecwam_(dp|sp)" )
  set( ECWAM_PROJECT_NAME ${ecwam-source_dir} )
endif()

if( NOT ECWAM_PROJECT_NAME )
  set( ECWAM_PROJECT_NAME ecwam )
endif()

string( TOUPPER ${ECWAM_PROJECT_NAME} PNAME )
set( ${PNAME}_ECBUILD_COMPILE_FLAGS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/ecwam_init_flags.cmake CACHE FILEPATH "")
project( ${ECWAM_PROJECT_NAME} LANGUAGES Fortran C CXX )

include( ecwam_macros )
ecbuild_enable_fortran( REQUIRED NO_MODULE_DIRECTORY )

### Dependencies

ecbuild_find_package( fiat      REQUIRED )
ecbuild_find_package( eccodes   REQUIRED )
ecbuild_find_package( field_api REQUIRED )

### Options

ecbuild_add_option( FEATURE OMP
                    DESCRIPTION "support for OpenMP shared memory parallelism"
                    REQUIRED_PACKAGES "OpenMP COMPONENTS Fortran" )

ecbuild_add_option( FEATURE MPI
                    DESCRIPTION "Support for MPI distributed parallelism"
                    REQUIRED_PACKAGES "MPI COMPONENTS Fortran" )

ecbuild_add_option( FEATURE MULTIO
                    DESCRIPTION "Support for IO to the FDB"
                    REQUIRED_PACKAGES "multio" )

ecbuild_add_option( FEATURE SINGLE_PRECISION
                    DEFAULT OFF
                    DESCRIPTION "Support for Single Precision" )

ecbuild_add_option( FEATURE UNWAM
                    DEFAULT OFF
                    DESCRIPTION "Support for UNWAM" )

### Determine ocean model precision
if( NOT DEFINED OCEAN_PREC )
  if( HAVE_SINGLE_PRECISION )
     set(OCEAN_PREC SP)
  else()
     set(OCEAN_PREC DP)
  endif()
endif()
string( TOLOWER "${OCEAN_PREC}" ocean_prec )

ecbuild_add_option( FEATURE OCEAN_COUPLING 
                    DEFAULT ON
                    DESCRIPTION "Support ocean coupling"
                    REQUIRED_PACKAGES "nemo_${ocean_prec}" )

ecbuild_add_option( FEATURE ECWAM_LAUNCH
                    DEFAULT ON
                    DESCRIPTION "Use share/ecwam/scripts/ecwam-launch script to run ecwam tests" )

ecbuild_add_option( FEATURE ECFLOW
                    DESCRIPTION "ecflow meter updates"
                    REQUIRED_PACKAGES "ecflow_light")

if( HAVE_ECWAM_LAUNCH )
  # Overwrite meaning of mpiexec
  set( ECWAM_LAUNCH ${CMAKE_CURRENT_SOURCE_DIR}/share/ecwam/scripts/ecwam-launch )
  set( MPIEXEC_EXECUTABLE ${ECWAM_LAUNCH} )
  set( MPIEXEC_NUMPROC_FLAG   -n )
  set( MPIEXEC_NUMTHREAD_FLAG -c )
  set( CMAKE_CROSSCOMPILING_EMULATOR ${ECWAM_LAUNCH} ) # For serial runs
endif()

### Find python dependencies

ecwam_find_python_mods()

### Set third party libraries

if( HAVE_SINGLE_PRECISION )
  set( prec sp )
else()
  set( prec dp )
endif()

unset( OpenMP_Fortran_LIBRARIES )
if( HAVE_OMP )
  set( OpenMP_Fortran_LIBRARIES OpenMP::OpenMP_Fortran )
endif()

### Apply Loki source-to-source transformations to the ecWam physics
ecbuild_add_option( FEATURE LOKI
                    DEFAULT OFF
                    DESCRIPTION "Apply Loki transformations to ecWam physics"
                    REQUIRED_PACKAGES "loki" )

# Set default Loki transformation mode
set( LOKI_MODE "scc-hoist" CACHE STRING "Transformation mode for Loki source transformations" )

### OpenACC
if( ${CMAKE_VERSION} VERSION_LESS "3.25" AND (NOT DEFINED ENABLE_ACC OR ENABLE_ACC ) )
    # Incredibly inconvenient: FindOpenACC does _not_ set OpenACC_FOUND, only
    # the language-specific components OpenACC_Fortran_FOUND and OpenACC_C_FOUND.
    # This means, even internally CMake considers OpenACC as not found.
    # (See eg get_property(... GLOBAL PROPERTY PACKAGES_NOT_FOUND))
    # Therefore, we search for OpenACC, set OpenACC_FOUND ourselves according to
    # the result, and then, trigger a second find_package via ecbuild_add_option.
    # This then conveniently takes the previously set OpenACC_FOUND into account
    # and rectifies CMake's internal bookkeeping in the process.
    # This has been fixed in CMake 3.25
    find_package( OpenACC )
    if( OpenACC_Fortran_FOUND AND OpenACC_C_FOUND )
        set( OpenACC_FOUND ON )
    endif()
endif()
ecbuild_add_option( FEATURE ACC
                    DESCRIPTION "OpenACC" DEFAULT OFF
                    REQUIRED_PACKAGES "OpenACC" 
                    CONDITION HAVE_LOKI AND NOT LOKI_MODE MATCHES "idem|idem-stack" )

### CUDA-aware MPI
ecbuild_add_option( FEATURE GPU_AWARE_MPI
                    DEFAULT OFF
                    DESCRIPTION "Enable GPU-aware MPI" 
                    CONDITION HAVE_ACC AND MPI_Fortran_HAVE_F08_MODULE )

### Sources

include(ecwam_compile_flags)
add_subdirectory(src)
add_subdirectory(share)

### Tests

add_subdirectory(tests)

get_property(ecwam_test_configs GLOBAL PROPERTY ecwam_test_configs)
ecwam_cache_bathymetry( ${ecwam_test_configs} )

### Documentation

add_subdirectory(doc)

### Export

ecbuild_install_project( NAME ${PROJECT_NAME} )

ecbuild_print_summary()

