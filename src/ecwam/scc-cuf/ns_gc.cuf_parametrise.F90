! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE NS_GC_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) FUNCTION NS_GC_CUF_PARAMETRISE (USTAR)
    
    ! ----------------------------------------------------------------------
    
    !**** *NS_GC* - FUNCTION TO DETERMINE THE CUT-OFF ANGULAR FREQUENCY INDEX
    !               FOR THE GRAVITY-CAPILLARY MODEL
    !               !!!! rounded to the closest index of XK_GC  !!!!!
    
    !**   INTERFACE.
    !     ----------
    
    !       *FUNCTION* *NS_GC (USTAR)*
    
    !       *USTAR*  - FRICTION VELOCITY.
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWFRED, ONLY: XLOGKRATIOM1_GC, NWAV_GC_D, XKM_GC_D
    USE YOWPCONS, ONLY: SQRTGOSURFT_D
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER :: NS_GC_CUF_PARAMETRISE
    REAL(KIND=JWRB), INTENT(IN) :: USTAR
    
!$loki routine seq
    REAL(KIND=JWRB) :: Y, XKS
    
    ! ----------------------------------------------------------------------
    
    
    !!!Y = 1.0_JWRB/(1.48_JWRB+2.05_JWRB*UST)
    !!!Y = (1.0_JWRB + UST**2)/(1.0_JWRB+10.0_JWRB*UST**2)
    
    XKS = SQRTGOSURFT_D / (1.48_JWRB + 2.05_JWRB*USTAR)
    
    NS_GC_CUF_PARAMETRISE = MIN(INT(LOG(XKS*XKM_GC_D(1))*XLOGKRATIOM1_GC) + 1, NWAV_GC_D - 1)
    
    
  END FUNCTION NS_GC_CUF_PARAMETRISE
END MODULE NS_GC_CUF_PARAMETRISE_MOD
