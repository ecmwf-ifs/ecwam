# (C) Copyright 2022- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

if(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  set(autopromote_flags   "-sreal64")
  set(convertendian_flags "-hbyteswapio")
  set(checkbounds_flags   "-Rb")
  set(fpe_flags           "-Ktrap=fp")
  set(initsnan_flags      "-ei")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(autopromote_flags   "-fdefault-real-8 -fdefault-double-8")
  set(convertendian_flags "-fconvert=big-endian" )
  set(checkbounds_flags   "-fcheck=bounds")
  set(fpe_flags           "-ffpe-trap=invalid,zero,overflow")
  set(initsnan_flags      "-finit-real=snan")
  if( NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10 )
    ecbuild_add_fortran_flags( "-fallow-argument-mismatch" NAME argument_mismatch )
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(autopromote_flags   "-real-size 64")
  set(convertendian_flags "-convert big_endian")
  set(checkbounds_flags   "-check bounds")
  set(initsnan_flags      "-init=snan")
  set(fpe_flags           "-fpe0")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC")
  set(autopromote_flags   "-r8")
  set(convertendian_flags "-Mbyteswapio")
  set(checkbounds_flags   "-Mbounds")
  set(fpe_flags           "-Ktrap=fp")
endif()

if( NOT HAVE_SINGLE_PRECISION )
  ecbuild_add_fortran_flags( "${autopromote_flags}"   NAME autopromote )
endif()
#ecbuild_add_fortran_flags( "${convertendian_flags}" NAME convert_bigendian )

if( CMAKE_BUILD_TYPE MATCHES "Debug" )
  if( fpe_flags )
    ecbuild_add_fortran_flags( "${fpe_flags}"         NAME fpe)
  endif()
  if( initsnan_flags )
    ecbuild_add_fortran_flags( "${initsnan_flags}"    NAME initsnan)
  endif()
  if( checkbounds_flags )
    ecbuild_add_fortran_flags( "${checkbounds_flags}" NAME checkbounds)
  endif()
endif()
