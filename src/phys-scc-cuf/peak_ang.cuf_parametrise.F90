! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE PEAK_ANG_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE PEAK_ANG_CUF_PARAMETRISE (KIJS, KIJL, FL1, XNU, SIG_TH, IJ)
    
    !***  *PEAK_ANG*   DETERMINES ANGULAR WIDTH NEAR PEAK OF SPECTRUM
    
    !     PETER JANSSEN
    
    !     PURPOSE.
    !     --------
    
    !              DETERMINATION OF PEAK PARAMETERS
    
    !     INTERFACE.
    !     ----------
    !              *CALL*  *PEAK_ANG(KIJS,KIJL,FL1,XNU,SIG_TH)*
    
    !               INPUT:
    !                  *KIJS*   - FIRST GRIDPOINT
    !                  *KIJL*   - LAST GRIDPOINT
    !                  *FL1*    - SPECTRUM
    !               OUTPUT:
    !                  *XNU*    - RELATIVE SPECTRAL WIDTH
    !                  *SIG_TH* - RELATIVE WIDTH IN DIRECTION
    
    !     METHOD.
    !     -------
    !              NONE
    
    !     EXTERNALS.
    !     ----------
    !              NONE
    
    !-----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWFRED, ONLY: FR_D, DFIM_D, DFIMFR_D, DFIMOFR_D, DFIMFR2_D, DELTH_D, TH_D, SINTH_D, COSTH_D, WETAIL, WP1TAIL, WP2TAIL,  &
    & FRATIO
    USE YOWPARAM, ONLY: NANG_D, NFRE_D
    
    
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    REAL(KIND=JWRB), INTENT(OUT), DEVICE :: XNU, SIG_TH
    
    
    INTEGER(KIND=JWIM) :: NSH
    INTEGER(KIND=JWIM) :: M, K
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: MMAX
    INTEGER(KIND=JWIM) :: MMSTART, MMSTOP
    REAL(KIND=JWRB), PARAMETER :: CONST_SIG = 1.0_JWRB
    REAL(KIND=JWRB) :: R1
    REAL(KIND=JWRB) :: DELT25, COEF_FR, COEF_FR2
    REAL(KIND=JWRB) :: ZEPSILON
    REAL(KIND=JWRB) :: SUM0, SUM1, SUM2, XMAX, TEMP
    REAL(KIND=JWRB) :: THMEAN, SUM_S, SUM_C
    
    ! ----------------------------------------------------------------------
    
    !***  1. DETERMINE L-H SPECTRAL WIDTH OF THE 2-D SPECTRUM.
    !     ---------------------------------------------------
    
    ZEPSILON = 10._JWRB*EPSILON(ZEPSILON)
    NSH = 1 + INT(LOG(1.5_JWRB) / LOG(FRATIO))
    
    SUM0 = ZEPSILON
    SUM1 = 0._JWRB
    SUM2 = 0._JWRB
    
    DO M=1,NFRE_D
      K = 1
      TEMP = FL1(IJ, K, M)
      DO K=2,NANG_D
        TEMP = TEMP + FL1(IJ, K, M)
      END DO
      SUM0 = SUM0 + TEMP*DFIM_D(M)
      SUM1 = SUM1 + TEMP*DFIMFR_D(M)
      SUM2 = SUM2 + TEMP*DFIMFR2_D(M)
    END DO
    
    !     ADD TAIL CORRECTIONS
    DELT25 = WETAIL*FR_D(NFRE_D)*DELTH_D
    COEF_FR = WP1TAIL*DELTH_D*FR_D(NFRE_D)**2
    COEF_FR2 = WP2TAIL*DELTH_D*FR_D(NFRE_D)**3
    SUM0 = SUM0 + DELT25*TEMP
    SUM1 = SUM1 + COEF_FR*TEMP
    SUM2 = SUM2 + COEF_FR2*TEMP
    
    IF (SUM0 > ZEPSILON) THEN
      XNU = SQRT(MAX(ZEPSILON, SUM2*SUM0 / SUM1**2 - 1._JWRB))
    ELSE
      XNU = ZEPSILON
    END IF
    
    !***  2. DETERMINE ANGULAR WIDTH OF THE 2-D SPECTRUM.
    !     ----------------------------------------------
    
    !     MAX OF 2-D SPECTRUM
    XMAX = 0._JWRB
    MMAX = 2
    
    DO M=2,NFRE_D - 1
      DO K=1,NANG_D
        IF (FL1(IJ, K, M) > XMAX) THEN
          MMAX = M
          XMAX = FL1(IJ, K, M)
        END IF
      END DO
    END DO
    
    SUM1 = ZEPSILON
    SUM2 = 0._JWRB
    
    MMSTART = MAX(1, MMAX - NSH)
    MMSTOP = MIN(NFRE_D, MMAX + NSH)
    DO M=MMSTART,MMSTOP
      SUM_S = 0._JWRB
      SUM_C = ZEPSILON
      DO K=1,NANG_D
        SUM_S = SUM_S + SINTH_D(K)*FL1(IJ, K, M)
        SUM_C = SUM_C + COSTH_D(K)*FL1(IJ, K, M)
      END DO
      THMEAN = ATAN2(SUM_S, SUM_C)
      DO K=1,NANG_D
        SUM1 = SUM1 + FL1(IJ, K, M)*DFIM_D(M)
        SUM2 = SUM2 + COS(TH_D(K) - THMEAN)*FL1(IJ, K, M)*DFIM_D(M)
      END DO
    END DO
    
    IF (SUM1 > ZEPSILON) THEN
      R1 = SUM2 / SUM1
      SIG_TH = CONST_SIG*SQRT(2._JWRB*(1._JWRB - R1))
    ELSE
      SIG_TH = 0._JWRB
    END IF
    
    
  END SUBROUTINE PEAK_ANG_CUF_PARAMETRISE
END MODULE PEAK_ANG_CUF_PARAMETRISE_MOD
