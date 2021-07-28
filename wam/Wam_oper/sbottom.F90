      SUBROUTINE SBOTTOM (IJS, IJL, KIJS, KIJL, GFL, FLD, SL, DEPTH)

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

!       *CALL* *SBOTTOM (IJS, IJL, KIJS, KIJL, GFL, FLD, SL, DEPTH)
!          *IJS:IJL* - 1st DIMENSION OF GFL
!          *KIJS*    - INDEX OF FIRST GRIDPOINT
!          *KIJL*    - INDEX OF LAST GRIDPOINT
!          *GFL*     - SPECTRUM.
!          *FLD*     - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*      - TOTAL SOURCE FUNCTION ARRAY
!          *DEPTH*   - WATER DEPTH

!     METHOD.
!     -------

!       SEE REFERENCES.

!     REFERENCES.
!     -----------

!       HASSELMANN ET AL, D. HYDR. Z SUPPL A12(1973) (JONSWAP)
!       BOUWS AND KOMEN, JPO 13(1983)1653-1658

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : GM1
      USE YOWSHAL  , ONLY : TFAK     ,INDEP   ,BATHYMAX
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: GFL

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH 

      INTEGER(KIND=JWIM):: IJ, K, M
      REAL(KIND=JWRB) :: CONST, ARG
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: SBO

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SBOTTOM',0,ZHOOK_HANDLE)

      CONST = -2.0_JWRB*0.038_JWRB*GM1
      DO M=1,NFRE
        DO IJ=KIJS,KIJL
          IF(DEPTH(IJ).LT.BATHYMAX) THEN
            ARG = 2.0_JWRB* DEPTH(IJ)*TFAK(INDEP(IJ),M)
            ARG = MIN(ARG,50.0_JWRB)
            SBO(IJ,M) = CONST*TFAK(INDEP(IJ),M)/SINH(ARG)
          ELSE
            SBO(IJ,M) = 0.0_JWRB
          ENDIF
        ENDDO
      ENDDO

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            SL(IJ,K,M) = SL(IJ,K,M)+SBO(IJ,M)*GFL(IJ,K,M)
            FLD(IJ,K,M) = FLD(IJ,K,M)+SBO(IJ,M)
          ENDDO
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SBOTTOM',1,ZHOOK_HANDLE)

      END SUBROUTINE SBOTTOM
