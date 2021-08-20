      SUBROUTINE PEAK (IJS, IJL, FETCH, FPMAX, U10)

! ----------------------------------------------------------------------

!**** *PEAK* - COMPUTE JONSWAP PARAMETERS FROM WINDSPEED.

!     S. HASSELMANN  - JULY 1990
!     H. GUNTHER     - DECEMBER 1990   MODIFIED FOR CYCLE_4.
!     J. BIDLOT      - FEBRUARY 1995   MESSAGE PASSING

!*    PURPOSE.
!     --------

!       COMPUTES FOR EACH GRID POINT OF A BLOCK THE PEAK FREQUENCY
!       FROM A FETCH LAW AND THE JONSWAP ALPHAJONS FROM THE ALPHAJONS NY
!       RELATION.

!**   INTERFACE.
!     ----------

!       *CALL* *PEAK (IJL, IJS, FETCH, FPMAX,U10)*
!          *IJS*     INTEGER  FIRST POINT IN BLOCK.
!          *IJL*     INTEGER  LAST  POINT IN BLOCK.
!          *FETCH*   REAL     FETCH TO BE USED (METRES).
!          *FPMAX*   REAL     MAXIMUM PEAK FREQUENCY (HERTZ).
!          *U10*  INTERMEDIATE STORAGE OF MODULUS OF WIND
!                    VELOCITY.

!     METHOD.
!     -------

!       FP = A * (G*FETCH/U_10**2)**D    A = 2.84
!       FP = MAX [FP, 0.13]              D = -3./10.
!       FP = MIN [FP, FRMAX*U_10/G]      b = 0.033
!       ALPHJ = B * FP-**2/3             B = 0.033
!       ALPHJ = MAX [ALPHJ, 0.0081]
!       FP = G/U_10*FP

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCES.
!     -----------

!       K.HASSELMAN,D.B.ROOS,P.MUELLER AND W.SWELL
!          A PARAMETRIC WAVE PREDICTION MODEL
!          JOURNAL OF PHSICAL OCEANOGRAPHY, VOL. 6, NO. 2, MARCH 1976.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWJONS  , ONLY : AJONS    ,BJONS    ,DJONS    ,EJONS    ,    &
     &            FP       ,ALPHJ
      USE YOWPCONS , ONLY : G

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL

      REAL(KIND=JWRB), INTENT(IN) :: FETCH, FPMAX
      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(IN) :: U10


      INTEGER(KIND=JWIM) :: IJ
      REAL(KIND=JWRB) :: GX, GXU, UG

! ----------------------------------------------------------------------

!*    1. COMPUTE VALUES FROM FETCH LAWS.
!        -------------------------------

      GX = G * FETCH
      DO IJ = IJS, IJL
        IF (U10(IJ) > 0.1E-08_JWRB) THEN
          GXU = GX/(U10(IJ)*U10(IJ))
          UG = G/U10(IJ)
          FP(IJ) = AJONS * GXU ** DJONS
          FP(IJ) = MAX(0.13_JWRB, FP(IJ))
          FP(IJ) = MIN(FP(IJ), FPMAX/UG)
          ALPHJ(IJ) = BJONS * FP(IJ)** EJONS
          ALPHJ(IJ) = MAX(ALPHJ(IJ), 0.0081_JWRB)
          FP(IJ) = FP(IJ)*UG
        ELSE
          FP(IJ) = 0.0_JWRB
        ENDIF
      ENDDO

      END SUBROUTINE PEAK
