      SUBROUTINE PEAK (IJS, IJL, IG, FETCH, FPMAX,U10OLD)

C ----------------------------------------------------------------------
C
C**** *PEAK* - COMPUTE JONSWAP PARAMETERS FROM WINDSPEED.
C
C     S. HASSELMANN  - JULY 1990
C     H. GUNTHER     - DECEMBER 1990   MODIFIED FOR CYCLE_4.
!     J. BIDLOT      - FEBRUARY 1995   MESSAGE PASSING
C
C*    PURPOSE.
C     --------
C
C       COMPUTES FOR EACH GRID POINT OF A BLOCK THE PEAK FREQUENCY
C       FROM A FETCH LAW AND THE JONSWAP ALPHAJONS FROM THE ALPHAJONS NY
C       RELATION.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *PEAK (IJL, IJS, IG, FETCH, FPMAX,U10OLD)*
C          *IJS*     INTEGER  FIRST POINT IN BLOCK.
C          *IJL*     INTEGER  LAST  POINT IN BLOCK.
C          *IG*      INTEGER  BLOCK NUMBER.
C          *FETCH*   REAL     FETCH TO BE USED (METRES).
C          *FPMAX*   REAL     MAXIMUM PEAK FREQUENCY (HERTZ).
C          *U10OLD*  INTERMEDIATE STORAGE OF MODULUS OF WIND
C                    VELOCITY.
C
C     METHOD.
C     -------
C
C       FP = A * (G*FETCH/U_10**2)**D    A = 2.84
C       FP = MAX [FP, 0.13]              D = -3./10.
C       FP = MIN [FP, FRMAX*U_10/G]      b = 0.033
C       ALPHJ = B * FP-**2/3             B = 0.033
C       ALPHJ = MAX [ALPHJ, 0.0081]
C       FP = G/U_10*FP
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCES.
C     -----------
C
C       K.HASSELMAN,D.B.ROOS,P.MUELLER AND W.SWELL
C          A PARAMETRIC WAVE PREDICTION MODEL
C          JOURNAL OF PHSICAL OCEANOGRAPHY, VOL. 6, NO. 2, MARCH 1976.
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comwind.h"
C
#include "comjons.h"
C
#include "commpp.h"
C
C ----------------------------------------------------------------------

!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD

C ----------------------------------------------------------------------
C
      PARAMETER  (G = 9.806)
      PARAMETER  (A = 2.84, B = 0.033, D = -(3./10.), E = 2./3.)
C
C           FETCH LAW CONSTANTS.
C
C ----------------------------------------------------------------------
C
C*    1. COMPUTE VALUES FROM FETCH LAWS.
C        -------------------------------
C
 1000 CONTINUE
      GX = G * FETCH
      DO 1001 IJ = IJS, IJL
         IF (U10OLD(IJ,IG) .GT. 0.1E-08) THEN
            U10 = U10OLD(IJ,IG)
            GXU = GX/(U10*U10)
            UG = G/U10
            FP(IJ) = A * GXU ** D
            FP(IJ) = MAX(0.13, FP(IJ))
            FP(IJ) = MIN(FP(IJ), FPMAX/UG)
            ALPHJ(IJ) = B * FP(IJ)** E
            ALPHJ(IJ) = MAX(ALPHJ(IJ), 0.0081)
            FP(IJ) = FP(IJ)*UG
         ENDIF
 1001 CONTINUE

      RETURN
      END
