      SUBROUTINE PEAK (IJS, IJL, IG, FETCH, FPMAX, U10OLD)

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

!       *CALL* *PEAK (IJL, IJS, IG, FETCH, FPMAX,U10OLD)*
!          *IJS*     INTEGER  FIRST POINT IN BLOCK.
!          *IJL*     INTEGER  LAST  POINT IN BLOCK.
!          *IG*      INTEGER  BLOCK NUMBER.
!          *FETCH*   REAL     FETCH TO BE USED (METRES).
!          *FPMAX*   REAL     MAXIMUM PEAK FREQUENCY (HERTZ).
!          *U10OLD*  INTERMEDIATE STORAGE OF MODULUS OF WIND
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

      USE YOWJONS  , ONLY : AJONS    ,BJONS    ,DJONS    ,EJONS    ,
     &            FP       ,ALPHJ
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NBLO
      USE YOWPCONS , ONLY : G

! ----------------------------------------------------------------------

!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD

! ----------------------------------------------------------------------

!*    1. COMPUTE VALUES FROM FETCH LAWS.
!        -------------------------------

      GX = G * FETCH
      DO IJ = IJS, IJL
        IF (U10OLD(IJ,IG) .GT. 0.1E-08) THEN
          U10 = U10OLD(IJ,IG)
          GXU = GX/(U10*U10)
          UG = G/U10
          FP(IJ) = AJONS * GXU ** DJONS
          FP(IJ) = MAX(0.13, FP(IJ))
          FP(IJ) = MIN(FP(IJ), FPMAX/UG)
          ALPHJ(IJ) = BJONS * FP(IJ)** EJONS
          ALPHJ(IJ) = MAX(ALPHJ(IJ), 0.0081)
          FP(IJ) = FP(IJ)*UG
        ELSE
          FP(IJ) = 0.
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE PEAK
