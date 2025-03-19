# (C) Copyright 2022- ECMWF.
# (C) Copyright 2022- Meteo-France.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

if(OpenMP_FOUND)

   try_compile(
       _HAVE_OMP_OFFLOAD
       ${CMAKE_CURRENT_BINARY_DIR}
       ${CMAKE_CURRENT_SOURCE_DIR}/cmake/features/OMP/test_openmp_target.F90
       LINK_LIBRARIES OpenMP::OpenMP_Fortran
       LINK_OPTIONS SHELL:${OpenMP_Fortran_FLAGS}
       OUTPUT_VARIABLE _HAVE_OMP_OFFLOAD_OUTPUT
   )

endif()