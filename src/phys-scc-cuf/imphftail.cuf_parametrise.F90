! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE IMPHFTAIL_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE IMPHFTAIL_CUF_PARAMETRISE (KIJS, KIJL, MIJ, FLM, WAVNUM, XK2CG, FL1, IJ)
    ! ----------------------------------------------------------------------
    
    !**** *IMPHFTAIL* - IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM
    
    
    !*    PURPOSE.
    !     --------
    
    !     IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM ABOVE FREQUENCY INDEX MIJ
    
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *IMPHFTAIL (KIJS, KIJL, MIJ, FLM, WAVNUM, XK2CG, FL1)
    !          *KIJS*    - INDEX OF FIRST GRIDPOINT
    !          *KIJL*    - INDEX OF LAST GRIDPOINT
    !          *MIJ*     - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
    !          *FLM*     - SPECTAL DENSITY MINIMUM VALUE
    !          *WAVNUM*  - WAVENUMBER
    !          *XK2CG*   - (WAVNUM)**2 * GROUP SPEED
    !          *FL1*     - SPECTRUM (INPUT AND OUTPUT).
    
    !     METHOD.
    !     -------
    
    !     EXTERNALS.
    !     ---------
    
    !     REFERENCE.
    !     ----------
    
    ! ----------------------------------------------------------------------
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWPARAM, ONLY: NANG_D, NFRE_D
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 24
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    INTEGER(KIND=JWIM), INTENT(IN), DEVICE, DIMENSION(KIJL) :: MIJ
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(NANG_LOKI_PARAM) :: FLM
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NFRE_LOKI_PARAM) :: WAVNUM, XK2CG
    REAL(KIND=JWRB), INTENT(INOUT), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    
    
    INTEGER(KIND=JWIM) :: K, M
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    REAL(KIND=JWRB) :: AKM1, TFAC
    REAL(KIND=JWRB) :: TEMP1, TEMP2
    
    ! ----------------------------------------------------------------------
    
    
    !*    DIAGNOSTIC TAIL.
    !     ----------------
    
    TEMP1 = 1.0_JWRB / XK2CG(IJ, MIJ(IJ)) / WAVNUM(IJ, MIJ(IJ))
    
    DO M=MIJ(IJ) + 1,NFRE_D
      TEMP2 = 1.0_JWRB / XK2CG(IJ, M) / WAVNUM(IJ, M)
      TEMP2 = TEMP2 / TEMP1
      
      !*    MERGE TAIL INTO SPECTRA.
      !     ------------------------
      DO K=1,NANG_D
        TFAC = FL1(IJ, K, MIJ(IJ))
        FL1(IJ, K, M) = MAX(TEMP2*TFAC, FLM(K))
      END DO
    END DO
    
    ! ----------------------------------------------------------------------
    
    
  END SUBROUTINE IMPHFTAIL_CUF_PARAMETRISE
END MODULE IMPHFTAIL_CUF_PARAMETRISE_MOD
