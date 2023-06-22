! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SBOTTOM_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE SBOTTOM_CUF_PARAMETRISE (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, IJ)
    
    !SHALLOW
    ! ----------------------------------------------------------------------
    
    !**** *SBOTTOM* - COMPUTATION OF BOTTOM FRICTION.
    
    !     G.J.KOMEN AND Q.D.GAO
    !     OPTIMIZED BY L.F. ZAMBRESKY
    !     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
    
    !*    PURPOSE.
    !     --------
    
    !       COMPUTATION OF BOTTOM FRICTION DISSIPATION
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *SBOTTOM (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH)
    !          *KIJS*    - INDEX OF FIRST GRIDPOINT
    !          *KIJL*    - INDEX OF LAST GRIDPOINT
    !          *FL1*     - SPECTRUM.
    !          *FLD*     - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*      - TOTAL SOURCE FUNCTION ARRAY
    !          *WAVNUM*  - WAVE NUMBER
    !          *DEPTH*   - WATER DEPTH
    
    !     METHOD.
    !     -------
    
    !       SEE REFERENCES.
    
    !     REFERENCES.
    !     -----------
    
    !       HASSELMANN ET AL, D. HYDR. Z SUPPL A12(1973) (JONSWAP)
    !       BOUWS AND KOMEN, JPO 13(1983)1653-1658
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWPARAM, ONLY: NANG_D, NFRE_D, NFRE_RED_D
    USE YOWPCONS, ONLY: GM1_D
    USE YOWSHAL, ONLY: BATHYMAX
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    REAL(KIND=JWRB), INTENT(INOUT), DEVICE, DIMENSION(NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FLD, SL
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NFRE_LOKI_PARAM) :: WAVNUM
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: DEPTH
    
    
    INTEGER(KIND=JWIM) :: K, M
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    REAL(KIND=JWRB) :: CONST, ARG
    REAL(KIND=JWRB) :: SBO
    
    ! ----------------------------------------------------------------------
    
    
    CONST = -2.0_JWRB*0.038_JWRB*GM1_D
    DO M=1,NFRE_RED_D
      IF (DEPTH(IJ) < BATHYMAX) THEN
        ARG = 2.0_JWRB*DEPTH(IJ)*WAVNUM(IJ, M)
        ARG = MIN(ARG, 50.0_JWRB)
        SBO = CONST*WAVNUM(IJ, M) / SINH(ARG)
      ELSE
        SBO = 0.0_JWRB
      END IF
      
      DO K=1,NANG_D
        SL(K, M) = SL(K, M) + SBO*FL1(IJ, K, M)
        FLD(K, M) = FLD(K, M) + SBO
      END DO
    END DO
    
    
  END SUBROUTINE SBOTTOM_CUF_PARAMETRISE
END MODULE SBOTTOM_CUF_PARAMETRISE_MOD
