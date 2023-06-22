! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE CIWABR_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE CIWABR_CUF_PARAMETRISE (KIJS, KIJL, CICOVER, FL1, WAVNUM, CGROUP, CIWAB, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *CIWABR* - COMPUTE SEA ICE WAVE ATTENUATION FACTORS DUE TO ICE FLOES
    !                BOTTOM FRICTION.
    
    !*    PURPOSE.
    !     --------
    
    !       CIWABR COMPUTES SEA ICE WAVE ATTENUATION FACTORS DUE TO ICE FLOES
    !              BOTTOM FRICTION.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *CIWABR (KIJS,KIJL,CICOVER,FL1,CIWAB)
    
    !          *KIJS*     - INDEX OF FIRST POINT.
    !          *KIJL*     - INDEX OF LAST POINT.
    !          *CICOVER*  -SEA ICE COVER.
    !          *FL1*      -ENERGY SPECTRUM.
    !          *CIWAB*    -SEA ICE WAVE ATTENUATION FACTOR DUE TO ICE FLOE BOTTOM FRICTION
    
    !     METHOD.
    !     -------
    
    !     EXTERNALS.
    !     ----------
    
    !     REFERENCES.
    !     -----------
    
    !     KOHOUT A., M. MEYLAN, D PLEW, 2011: ANNALS OF GLACIOLOGY, 2011.
    
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWFRED, ONLY: FR_D, DFIM_D, DELTH_D
    USE YOWICE, ONLY: LICERUN_D, LMASKICE_D, CDICWA_D
    USE YOWPARAM, ONLY: NANG_D, NFRE_D
    USE YOWPCONS, ONLY: G_D, ZPI_D, ZPI4GM2_D, EPSMIN
    USE YOWSTAT, ONLY: IDELT_D
    
    
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: CICOVER
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NFRE_LOKI_PARAM) :: WAVNUM, CGROUP
    REAL(KIND=JWRB), INTENT(OUT), DEVICE, DIMENSION(NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: CIWAB
    
    
    INTEGER(KIND=JWIM) :: K, M
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    REAL(KIND=JWRB) :: EWH
    REAL(KIND=JWRB) :: X, ALP
    REAL(KIND=JWRB) :: XK2
    
    ! ----------------------------------------------------------------------
    
    
    IF (.not.LICERUN_D .or. LMASKICE_D) THEN
      
      DO M=1,NFRE_D
        DO K=1,NANG_D
          CIWAB(K, M) = 1.0_JWRB
        END DO
      END DO
      
    ELSE
      
      DO M=1,NFRE_D
        DO K=1,NANG_D
          EWH = 4.0_JWRB*SQRT(MAX(EPSMIN, FL1(IJ, K, M)*DFIM_D(M)))
          XK2 = WAVNUM(IJ, M)**2
          ALP = CDICWA_D*XK2*EWH
          X = ALP*CGROUP(IJ, M)*IDELT_D
          CIWAB(K, M) = 1.0_JWRB - CICOVER(IJ)*(1.0_JWRB - EXP(-MIN(X, 50.0_JWRB)))
        END DO
      END DO
      
    END IF
    
    
  END SUBROUTINE CIWABR_CUF_PARAMETRISE
END MODULE CIWABR_CUF_PARAMETRISE_MOD
