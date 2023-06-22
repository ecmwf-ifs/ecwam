! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDIWBK_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE SDIWBK_CUF_PARAMETRISE (KIJS, KIJL, FL1, FLD, SL, DEPTH, EMAXDPT, EMEAN, F1MEAN, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *SDIWBK* - COMPUTATION OF BOTTOM-INDUCED WAVE BREAKING DISSIPATION
    
    
    !*    PURPOSE.
    !     --------
    !       COMPUTE BOTTOM-INDUCED DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
    !       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
    !       OF DISSIPATION SOURCE FUNCTION.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *SDIWBK (KIJS, KIJL, FL1, FLD, SL, DEPTH, EMAXDPT, EMEAN, F1MEAN)*
    !          *KIJS*    - INDEX OF FIRST GRIDPOINT
    !          *KIJL*    - INDEX OF LAST GRIDPOINT
    !          *FL1*     - SPECTRUM.
    !          *FLD*     - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*      - TOTAL SOURCE FUNCTION ARRAY
    !          *DEPTH*   - WATER DEPTH
    !          *EMAXDPT* - MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
    !          *EMEAN*   - MEAN ENERGY DENSITY
    !          *F1MEAN*  - MEAN FREQUENCY BASED ON 1st MOMENT.
    
    !     METHOD.
    !     -------
    
    !       SEE REFERENCES.
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWPARAM, ONLY: NANG_D, NFRE_D, NFRE_RED_D
    USE YOWSTAT, ONLY: LBIWBK_D
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 24
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    REAL(KIND=JWRB), INTENT(INOUT), DEVICE, DIMENSION(NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FLD, SL
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: DEPTH, EMAXDPT
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: EMEAN, F1MEAN
    
    
    INTEGER(KIND=JWIM) :: K, M, IC
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    REAL(KIND=JWRB) :: ALPH, ARG, Q, Q_OLD, REL_ERR, EXPQ
    REAL(KIND=JWRB) :: SDS
    
    REAL, PARAMETER :: ALPH_B_J = 1.0_JWRB
    REAL, PARAMETER :: COEF_B_J = 2*ALPH_B_J
    REAL, PARAMETER :: DEPTHTRS = 50.0_JWRB
    
    ! ----------------------------------------------------------------------
    
    
    !*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
    !*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
    !        --------------------------------------------------------------
    
    IF (LBIWBK_D) THEN
      !       (FOLLOWING BATTJES-JANSSEN AND BEJI)
      IF (DEPTH(IJ) < DEPTHTRS) THEN
        ALPH = 2.0_JWRB*EMAXDPT(IJ) / EMEAN
        ARG = MIN(ALPH, 50.0_JWRB)
        Q_OLD = EXP(-ARG)
        !            USE NEWTON-RAPHSON METHOD
        DO IC=1,15
          EXPQ = EXP(-ARG*(1.0_JWRB - Q_OLD))
          Q = Q_OLD - (EXPQ - Q_OLD) / (ARG*EXPQ - 1.0_JWRB)
          REL_ERR = ABS(Q - Q_OLD) / Q_OLD
          IF (REL_ERR < 0.00001_JWRB) EXIT
          Q_OLD = Q
        END DO
        Q = MIN(Q, 1.0_JWRB)
        SDS = COEF_B_J*ALPH*Q*F1MEAN
      END IF
      
      DO M=1,NFRE_RED_D
        DO K=1,NANG_D
          IF (DEPTH(IJ) < DEPTHTRS) THEN
            SL(K, M) = SL(K, M) - SDS*FL1(IJ, K, M)
            FLD(K, M) = FLD(K, M) - SDS
          END IF
        END DO
      END DO
      
    END IF
    
    
  END SUBROUTINE SDIWBK_CUF_PARAMETRISE
END MODULE SDIWBK_CUF_PARAMETRISE_MOD
