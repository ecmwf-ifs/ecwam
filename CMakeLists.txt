# (C) Copyright 2022- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )
find_package( ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )

if( NOT ECWAM_PROJECT_NAME )
  set( ECWAM_PROJECT_NAME ecwam )
endif()

project( ${ECWAM_PROJECT_NAME} LANGUAGES Fortran C CXX )

include( ecwam_macros )
ecbuild_enable_fortran( REQUIRED NO_MODULE_DIRECTORY )

### Dependencies

ecbuild_find_package( fiat    REQUIRED )
ecbuild_find_package( eccodes REQUIRED )


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

ecbuild_add_option( FEATURE OCEAN_COUPLING 
                    DEFAULT ON
                    DESCRIPTION "Support ocean coupling"
                    CONDITION ${PNAME}_OCEANMODEL_LIBRARIES )

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

### Sources

include(ecwam_compile_flags)
add_subdirectory(src)

### Tests

add_subdirectory(tests)

### Documentation

add_subdirectory(doc)

### Export

ecbuild_install_project( NAME ${PROJECT_NAME} )

ecbuild_print_summary()
