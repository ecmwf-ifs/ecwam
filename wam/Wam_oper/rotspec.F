      SUBROUTINE ROTSPEC (NFRE, NANG, ML, KL, FL1, FL3, RTHET)

C----------------------------------------------------------------------
C
C**** *ROTSPEC* - ROUTINE TO ROTATE THE SPECTRUM.
C
C     EVA BAUER      MPI  HAMBURG    MAY 1990.
C     H. GUNTHER          GKSS/ECMWF JAN. 1991   MODIFIED FOR CYCLE_4
C
C*    PURPOSE.
C     -------
C
C       TO ROTATE THE SPECTRUM.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *ROTSPEC (NFRE, NANG, ML, KL, FL1, FL3, RTHET)*
C          *NFRE*  - FREQUENCY DIMENSION OF SPECTRA.
C          *NANG*  - DIRECTION DIMENSION OF SPECTRA.
C          *ML*    - NUMBER OF FREQUENCIES.
C          *KL*    - NUMBER OF DIRECTIONS.
C          *FL1*   - SPECTRUM TO BE ROTATED.
C          *FL3*   - ROTATED SPECTRUM.
C          *RTHET* - TURNING ANGLE IN RADIANS, CLOCKWISE.
C
C     METHOD.
C     -------
C
C       NONE.
C
C     EXTERNALS.
C     ---------
C
C       NONE.
C
C     REFERENCES.
C     -----------
C
C       NONE.
C
C -------------------------------------------------------------------
C
#include "parcons.h"
C
      DIMENSION FL1(NANG,NFRE), FL3(NANG,NFRE)
      DIMENSION AR1(KL,ML), AR2(ML)
C
C --------------------------------------------------------------------
C
      FTH = MOD(RTHET+ZPI,ZPI)
      FTH = FTH * REAL(KL) / ZPI
      INC = IFIX(FTH)
      ADIF = FTH - INC
      BDIF = 1. - ADIF
C
      DO 1000 K=1,KL
         KC  = K  - INC
         IF (KC .LT. 1) KC  = KC + KL
         KC1 = KC - 1
         IF (KC1.LT. 1) KC1 = KC1+ KL
C
         DO 2000 M=1,ML
            FL3(K,M) = BDIF * FL1(KC,M) + ADIF * FL1(KC1,M)
 2000    CONTINUE
 1000 CONTINUE
C
      RETURN
      END
