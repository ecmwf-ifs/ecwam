# (C) Copyright 2022- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

foreach( prec sp dp )
  if( HAVE_${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_aqua_era5_O48 CONFIG aqua_era5_O48.yml PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_aqua_era5_O48 CONFIG aqua_era5_O48.yml MPI 4 OMP 1 PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_aqua_era5_O48 CONFIG aqua_era5_O48.yml MPI 1 OMP 4 PREC ${prec} )

    ecwam_add_test( ${PROJECT_NAME}_test_aqua_oper_an_fc_O48 CONFIG aqua_oper_an_fc_O48.yml PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_aqua_oper_an_fc_O48 CONFIG aqua_oper_an_fc_O48.yml MPI 4 OMP 1 PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_aqua_oper_an_fc_O48 CONFIG aqua_oper_an_fc_O48.yml MPI 1 OMP 4 PREC ${prec} )

    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_era5_O48 CONFIG etopo1_era5_O48.yml PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_era5_O48 CONFIG etopo1_era5_O48.yml MPI 4 OMP 1 PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_era5_O48 CONFIG etopo1_era5_O48.yml MPI 1 OMP 4 PREC ${prec} )

    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_O48 CONFIG etopo1_oper_an_fc_O48.yml PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_O48 CONFIG etopo1_oper_an_fc_O48.yml MPI 4 OMP 1 PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_O48 CONFIG etopo1_oper_an_fc_O48.yml MPI 1 OMP 4 PREC ${prec} )

    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_O48_iphys_0 CONFIG etopo1_oper_an_fc_O48_iphys_0.yml PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_O48_iphys_0 CONFIG etopo1_oper_an_fc_O48_iphys_0.yml MPI 4 OMP 1 PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_O48_iphys_0 CONFIG etopo1_oper_an_fc_O48_iphys_0.yml MPI 1 OMP 4 PREC ${prec} )

    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_O48_cy49r1 CONFIG etopo1_oper_an_fc_O48_cy49r1.yml PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_O48_cy49r1 CONFIG etopo1_oper_an_fc_O48_cy49r1.yml MPI 4 OMP 1 PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_O48_cy49r1 CONFIG etopo1_oper_an_fc_O48_cy49r1.yml MPI 1 OMP 4 PREC ${prec} )

    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_O48_cy50r1 CONFIG etopo1_oper_an_fc_O48_cy50r1.yml PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_O48_cy50r1 CONFIG etopo1_oper_an_fc_O48_cy50r1.yml MPI 4 OMP 1 PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_O48_cy50r1 CONFIG etopo1_oper_an_fc_O48_cy50r1.yml MPI 1 OMP 4 PREC ${prec} )

  endif()
endforeach()
