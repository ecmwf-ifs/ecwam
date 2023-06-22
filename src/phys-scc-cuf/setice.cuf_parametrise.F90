! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!-----------------------------------------------------------------------
MODULE SETICE_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE SETICE_CUF_PARAMETRISE (KIJS, KIJL, FL1, CICOVER, COSWDIF, IJ)
    
    !-----------------------------------------------------------------------
    
    !**** *SETICE* ROUTINE TO SET SPECTRA ON ICE TO NOISE LEVEL.
    
    !     R.PORTZ      MPI         OKT.1992
    !     J. BIDLOT    ECMWF       JUNE 1996  MESSAGE PASSING
    
    !     PURPOSE.
    !     -------
    
    !          *SETICE* SET ICE SPECTRA (FL1) TO NOISE LEVEL
    
    !**   INTERFACE.
    !     ----------
    
    !         *CALL* *SETICE(KIJS, KIJL, FL1, CICOVER, WSWAVE, COSWDIF)*
    !          *KIJS*    - LOCAL INDEX OF FIRST GRIDPOINT
    !          *KIJL*    - LOCAL  INDEX OF LAST GRIDPOINT
    !          *FL1*     - SPECTRA
    !          *CICOVER* - SEA ICE COVER
    !          *WSWAVE*  - WIND SPEED.
    !          *COSWDIF* - COS(TH(K)-WDWAVE(IJ))
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWICE, ONLY: FLMIN, CITHRSH_D
    USE YOWPARAM, ONLY: NANG_D, NFRE_D
    USE YOWPCONS, ONLY: EPSMIN
    
    
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(INOUT), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: CICOVER
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(NANG_LOKI_PARAM) :: COSWDIF
    
    
    INTEGER(KIND=JWIM) :: M, K
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    REAL(KIND=JWRB) :: CIREDUC, TEMP, ICEFREE
    ! ----------------------------------------------------------------------
    
    
    !*    1. SET SPECTRA TO NOISE LEVEL OVER ICE POINTS.
    !     ----------------------------------------------
    
    IF (CICOVER(IJ) > CITHRSH_D) THEN
      CIREDUC = MAX(EPSMIN, (1.0_JWRB - CICOVER(IJ)))
      ICEFREE = 0.0_JWRB
    ELSE
      CIREDUC = 0.0_JWRB
      ICEFREE = 1.0_JWRB
    END IF
    
    TEMP = CIREDUC*FLMIN
    DO M=1,NFRE_D
      DO K=1,NANG_D
        FL1(IJ, K, M) = FL1(IJ, K, M)*ICEFREE + TEMP*MAX(0.0_JWRB, COSWDIF(K))**2
      END DO
    END DO
    
    
  END SUBROUTINE SETICE_CUF_PARAMETRISE
END MODULE SETICE_CUF_PARAMETRISE_MOD
