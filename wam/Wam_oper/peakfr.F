      SUBROUTINE PEAKFR (F, IJS, IJL, PEAKF)

C ----------------------------------------------------------------------
C
C**** *PEAKFR* - COMPUTATION OF PEAK FREQUENCY AT EACH GRID POINT.
C
C     H. GUNTHER      ECMWF            DECEMBER 1989
C     (CODE REMOVED FROM SUB. FEMEAN)
C
C*    PURPOSE.
C     --------
C
C       COMPUTE PEAK FREQUENCY AT EACH GRID POINT.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *PEAKFR (F, IJS, IJL, PEAKF)*
C          *F*     - SPECTRUM.
C          *IJS*   - INDEX OF FIRST GRIDPOINT
C          *IJL*   - INDEX OF LAST GRIDPOINT
C          *PEAKF* - PEAKFREQUENCY
C
C     METHOD.
C     -------
C
C       NONE.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comfred.h"
C
#include "commpp.h"
C
      DIMENSION TEMP1(NINF:NSUP), TEMP2(NINF:NSUP)
      DIMENSION F(NINF-1:NSUP,NANG,NFRE), PEAKF(NIBLO)
C
C ----------------------------------------------------------------------
C
C*    1. INITIALIZE ARRAYS
C        -----------------
C
      DO 1000 IJ = IJS,IJL
         TEMP1(IJ) = 0.
         PEAKF(IJ) = 0.
 1000 CONTINUE
C
C*    2. LOOP OVER FREQUENCIES
C        ---------------------
C
      DO 2000 M = 1,NFRE
C
C*    2.1 COMPUTE 1-D SPECTRUM (WITHOUT DELTA THETA)
C        -------------------------------------------
C
         DO 2100 IJ = IJS,IJL
            TEMP2(IJ) = 0.
 2100    CONTINUE
         DO 2110 K = 1,NANG
            DO 2120 IJ = IJS,IJL
              TEMP2(IJ) = TEMP2(IJ) + F(IJ,K,M)
 2120       CONTINUE
 2110    CONTINUE
C
C*    2.2 DEFINE PEAK FREQUENCY
C         ---------------------
C
         DO 2200 IJ = IJS,IJL
            IF (TEMP1(IJ).LT.TEMP2(IJ)) THEN
               TEMP1(IJ) = TEMP2(IJ)
               PEAKF(IJ) = FR(M)
            ENDIF
 2200    CONTINUE

 2000 CONTINUE

      RETURN
      END
