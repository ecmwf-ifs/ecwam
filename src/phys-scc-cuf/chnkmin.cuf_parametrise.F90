! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE CHNKMIN_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) FUNCTION CHNKMIN_CUF_PARAMETRISE (U10)
    
    ! ----------------------------------------------------------------------
    
    !**** *CHNKMIN* - FUNCTION TO COMPUTE THE MINMUM CHARNOCK
    
    !*    PURPOSE.
    !     -------
    
    
    !**   INTERFACE.
    !     ----------
    
    !       *FUNCTION* *CHNKMIN (U10)*
    
    !     METHOD.
    !     -------
    
    !     CHNKMIN = ALPHAMIN + (ALPHA-ALPHAMIN)*0.5_JWRB*(1.0_JWRB-TANH(U10-A))
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    !       NONE.
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWPHYS, ONLY: ALPHA_D, ALPHAMIN_D, CHNKMIN_U_D
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
!$loki routine seq
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 24
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    REAL(KIND=JWRB) :: CHNKMIN_CUF_PARAMETRISE
    REAL(KIND=JWRB), INTENT(IN) :: U10
    
    ! ----------------------------------------------------------------------
    
    
    CHNKMIN_CUF_PARAMETRISE = ALPHAMIN_D + (ALPHA_D - ALPHAMIN_D)*0.5_JWRB*(1.0_JWRB - TANH(U10 - CHNKMIN_U_D))
    
    
  END FUNCTION CHNKMIN_CUF_PARAMETRISE
END MODULE CHNKMIN_CUF_PARAMETRISE_MOD
