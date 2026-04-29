# (C) Copyright 2022- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

####################################################################
# OpenMP FLAGS
####################################################################

set( OpenMP_Fortran_FLAGS       "-mp -mp=bind,allcores,numa" CACHE STRING "" FORCE)

####################################################################
# OpenAcc FLAGS
####################################################################

set( OpenACC_Fortran_FLAGS "-acc=gpu -gpu=cc90,lineinfo,fastmath,rdc" CACHE STRING "" FORCE)

####################################################################
# CUDA FLAGS
####################################################################

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 90)
endif()

#####################################################################
## COMMON FLAGS
#####################################################################
#
#set( ECBUILD_Fortran_FLAGS "-Mframe" )
#set( ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mbyteswapio" )
#set( ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mstack_arrays" )
#set( ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mrecursive" )
#set( ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Kieee" )
#set( ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mdaz" )
#
####################################################################
# DEFINE COMPILERS
####################################################################

set(CMAKE_C_COMPILER nvc)
set(CMAKE_Fortran_COMPILER nvfortran)
set(CMAKE_CXX_COMPILER nvc++)

