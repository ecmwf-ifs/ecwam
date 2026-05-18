# (C) Copyright 2022- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

### The file mubuf.F90, which is only used for "preproc" is sensitive to optimisations
#   possibly leading to different wam_grid_<1,2,3> files.
#   This in turn leads to non-neglibible differences
#   of average 'swh' when running "chief".

if( CMAKE_Fortran_COMPILER_ID MATCHES Intel )
  set_source_files_properties( mubuf.F90 PROPERTIES COMPILE_OPTIONS "-fp-model;strict" )
  set_source_files_properties( propconnect.F90 PROPERTIES COMPILE_OPTIONS "-fp-model;strict" )
  set_source_files_properties( stresso.F90 PROPERTIES COMPILE_OPTIONS "-fp-model;strict" )
  set_source_files_properties( create_wam_bathymetry_ETOPO1.F90 PROPERTIES COMPILE_FLAGS "-check nobounds" )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES GNU )
  # TODO: This is only needed due to a known bug in the vectorization characteristics of GNU 15, and should
  # fixed in GNU 16. See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=118955 for more information.
  set_source_files_properties( sdice3.F90 PROPERTIES COMPILE_OPTIONS "-fdisable-ipa-simdclone" )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC" AND CMAKE_BUILD_TYPE MATCHES "Bit")
  set_source_files_properties(
      w_maxh.F90 sbottom.F90 PROPERTIES COMPILE_FLAGS " -g -O1 -Mflushz -Mno-signed-zeros "
  )
  set_source_files_properties( mubuf.F90 PROPERTIES COMPILE_OPTIONS "-Mnofma" )
  set_source_files_properties( initgc.F90 PROPERTIES COMPILE_FLAGS " -g -O1 -Mflushz -Mno-signed-zeros " )
  set_source_files_properties( iniwcst.F90 PROPERTIES COMPILE_FLAGS " -g -O1 -Mflushz -Mno-signed-zeros " )
  set_source_files_properties( depthprpt.F90 PROPERTIES COMPILE_FLAGS " -g -O1 -Mflushz -Mno-signed-zeros " )
  if( HAVE_SINGLE_PRECISION )
     set_source_files_properties( aki.F90 PROPERTIES COMPILE_FLAGS " -g -O1 -Mflushz -Mno-signed-zeros " )
     set_source_files_properties( kurtosis.F90 PROPERTIES COMPILE_FLAGS " -g -O1 -Mflushz -Mno-signed-zeros " )
     set_source_files_properties( stat_nl.F90 PROPERTIES COMPILE_FLAGS " -g -O1 -Mflushz -Mno-signed-zeros " )
     set_source_files_properties( transf_bfi.F90 PROPERTIES COMPILE_FLAGS " -g -O1 -Mflushz -Mno-signed-zeros " )
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC" AND CMAKE_BUILD_TYPE MATCHES "Debug")
  if( DEFINED ${PNAME}_Fortran_FLAGS_DEBUG )
    string(REPLACE "-Ktrap=fp" "" ${PNAME}_Fortran_FLAGS_DEBUG ${${PNAME}_Fortran_FLAGS_DEBUG})
  endif()
  set_source_files_properties( outbeta.F90 PROPERTIES COMPILE_OPTIONS "${${PNAME}_Fortran_FLAGS_DEBUG} -Ktrap=divz")
  set_source_files_properties( secondhh.F90 PROPERTIES COMPILE_OPTIONS "${${PNAME}_Fortran_FLAGS_DEBUG} -Ktrap=inv,ovf")
endif()

### The file grib2wgrid.F90 is sensitive to optimizations in single precision builds.
#   This leads to non-neglibible differences
#   of average 'swh' when running "chief".

if( CMAKE_Fortran_COMPILER_ID MATCHES Intel AND HAVE_SINGLE_PRECISION )
  set_source_files_properties( grib2wgrid.F90 PROPERTIES COMPILE_OPTIONS "-fp-model;strict" )
endif()
