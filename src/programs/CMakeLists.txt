# (C) Copyright 2022- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

foreach( program
    chief
    bouint
    preproc
    preset
    create_wam_bathymetry
    create_wam_bathymetry_ETOPO1 )

    ecbuild_add_executable(
        TARGET          ${PROJECT_NAME}-${program}
        SOURCES         ${program}.F90
        LIBS            ${PROJECT_NAME} ${OpenMP_Fortran_LIBRARIES}
    )
    ecwam_target_compile_definitions_FILENAME(${PROJECT_NAME}-${program})

endforeach()
