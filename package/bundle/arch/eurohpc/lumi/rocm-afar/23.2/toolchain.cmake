# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

####################################################################
# OpenMP flags
####################################################################

set( OpenMP_Fortran_FLAGS     "-fopenmp --offload-arch=gfx90a" CACHE STRING "" )

if(NOT DEFINED CMAKE_HIP_ARCHITECTURES)
  set(CMAKE_HIP_ARCHITECTURES gfx90a)
endif()

set( ECBUILD_Fortran_FLAGS_BIT "-g -O2" )
