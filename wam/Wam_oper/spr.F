      SUBROUTINE SPR (KL, NANG, THETAQ, THETA, ST)

C ----------------------------------------------------------------------
C
C**** *SPR* - ROUTINE TO COMPUTE SPREADING FACTOR.
C
C     SUSANNE HASSELMANN  JULY 1986.
C
C*    PURPOSE.
C     --------
C
C       COMPUTATION OF COS**2 SPREADING FUNCTION.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *SPR (KL, NANG, THETAQ, THETA, ST)*
C          *KL*      INTEGER  NUMBER OF ANGULAR INTERVALS
C          *NANG*    INTEGER  DIMENSION FOR ANGULAR INTERVALS
C          *THETAQ*  REAL     MEAN WAVE DIRECTION
C          *THETA*   REAL     ANGULAR DIRECTIONS
C          *ST*      REAL     SPREADING FUNCTION
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
C
C     REFERENCES.
C     -----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
      DIMENSION THETA(NANG), ST(NANG)
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
      PARAMETER (PIH = PI/2., PIH3=3.*PI/2., ZDP=2./PI)
C
C     SPREAD FCT. WITH COS**2.
C
      DO 250 K=1,KL
         THE = COS(THETA(K)-THETAQ)
         IF (THE.GT.0.) THEN
            ST(K) = ZDP*THE**2
            IF(ST(K).LT.0.1E-08) ST(K)=0.
         ELSE
            ST(K) = 0.
         ENDIF
  250 CONTINUE

      RETURN
      END
