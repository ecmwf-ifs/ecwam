      REAL FUNCTION AKI (OM, BETA)

C ----------------------------------------------------------------------
C
C**** *AKI* - FUNCTION TO COMPUTE WAVE NUMBER.
C
C     G. KOMEN, P. JANSSEN   KNMI        01/06/1986
C
C*    PURPOSE.
C     -------
C
C       *AKI* COMPUTES THE WAVE NUMBER AS FUNCTION OF
C             CIRCULAR FREQUENCY AND WATER DEPTH.
C
C**   INTERFACE.
C     ----------
C
C       *FUNCTION* *AKI (OM, BETA)*
C          *OM*      - CIRCULAR FREQUENCY.
C          *BETA*    - WATER DEPTH.
C
C     METHOD.
C     -------
C
C       NEWTONS METHOD TO SOLVE THE DISPERSION RELATION IN SHALLOW
C       WATER.
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
#include "parcons.h"
C
C*    *PARAMETER*  RELATIVE ERROR LIMIT OF NEWTON'S METHOD.
C
      PARAMETER (EBS = 0.0001)
C
C ----------------------------------------------------------------------
C
C*    1. START VALUE:  MAXIMUM FROM DEEP  AND EXTREM SHALLOW WATER
C                      WAVE NUMBER.
C        ---------------------------------------------------------
C
 1000 CONTINUE
      AKM1=OM**2/(4.*G)
      AKM2=OM/(2.*SQRT(G*BETA))
      AO=AMAX1(AKM1,AKM2)
C
C ----------------------------------------------------------------------
C
C*    2. ITERATION LOOP.
C        ---------------
C
 2000 CONTINUE
      AKP = AO
      BO = BETA*AO
      IF (BO.GT.40.) THEN
        AKI = OM**2/G
      ELSE
        TH = G*AO*TANH(BO)
        STH = SQRT(TH)
        AO = AO+(OM-STH)*STH*2./(TH/AO+G*BO/COSH(BO)**2)
        IF (ABS(AKP-AO).GT.EBS*AO) GO TO 2000
        AKI = AO
      ENDIF

      RETURN
      END
