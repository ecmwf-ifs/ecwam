! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

       MODULE WAM_INIT_GPU_MOD
       CONTAINS
          SUBROUTINE WAM_INIT_GPU(IRANK)
#ifdef _OPENACC
             USE OPENACC
#endif
             USE PARKIND_WAVE, ONLY : JWIM
             IMPLICIT NONE
 
             INTEGER(KIND=JWIM), INTENT(IN) :: IRANK
             INTEGER :: DEVTYPE, DEVNUM, DEV
        
        
#ifdef _OPENACC
             DEVTYPE = ACC_GET_DEVICE_TYPE()
             DEVNUM = ACC_GET_NUM_DEVICES(DEVTYPE)
             DEV = MOD(IRANK-1, DEVNUM)
             CALL ACC_SET_DEVICE_NUM(DEV, DEVTYPE)
#endif
          END SUBROUTINE WAM_INIT_GPU
       END MODULE WAM_INIT_GPU_MOD
