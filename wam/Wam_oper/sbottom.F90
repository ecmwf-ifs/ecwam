      SUBROUTINE SBOTTOM (IJS, IJL, DPTH, F, FL, SL)

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

!       *CALL* *SBOTTOM (IJS, IJL, DPTH, F, FL, SL)
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *DEPTH* - WATER DEPTH
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY

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

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: DPTH 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL, SL

      INTEGER(KIND=JWIM):: IJ, K, M
      REAL(KIND=JWRB) :: CONST, ARG
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE) :: SBO

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SBOTTOM',0,ZHOOK_HANDLE)

      CONST = -2.0_JWRB*0.038_JWRB*GM1
      DO M=1,NFRE
        DO IJ=IJS,IJL
          IF(DPTH(IJ).LT.BATHYMAX) THEN
            ARG = 2.0_JWRB* DPTH(IJ)*TFAK(INDEP(IJ),M)
            ARG = MIN(ARG,50.0_JWRB)
            SBO(IJ,M) = CONST*TFAK(INDEP(IJ),M)/SINH(ARG)
          ELSE
            SBO(IJ,M) = 0.0_JWRB
          ENDIF
        ENDDO
      ENDDO

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            SL(IJ,K,M) = SL(IJ,K,M)+SBO(IJ,M)*F(IJ,K,M)
            FL(IJ,K,M) = FL(IJ,K,M)+SBO(IJ,M)
          ENDDO
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SBOTTOM',1,ZHOOK_HANDLE)

      END SUBROUTINE SBOTTOM
