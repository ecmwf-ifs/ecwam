      SUBROUTINE TAUHF (FRMAX)

C ----------------------------------------------------------------------
C
C**** *TAUHF* - COMPUTATION OF HIGH-FREQUENCY STRESS.
C
C     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90
C
C*    PURPOSE.
C     ---------
C
C       COMPUTE HIGH-FREQUENCY WAVE STRESS
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *TAUHF (FRMAX)*
C          *FRMAX - LAST MODEL FREQUENCY FR(ML).
C
C     METHOD.
C     -------
C
C       SEE REFERENCE FOR WAVE STRESS CALCULATION.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comcoup.h"
C
#include "comtabl.h"
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
      PARAMETER (JTOT = 50)
C
C ----------------------------------------------------------------------
C
C*    1. PRELIMINARY CALCULATIONS.
C        -------------------------
C
       USTARM = 5.
       ALPHAM = 10.*ALPHA
       DELUST = USTARM/FLOAT(IUSTAR)
       DELALP = ALPHAM/FLOAT(IALPHA)
C
       CONST1 = BETAMAX/XKAPPA**2
       OMEGAC = ZPI*FRMAX
C
       DO 1100 L=0,IALPHA
          DO 1200 K=0,IUSTAR
             TAUHFT(K,L) = 0.
 1200     CONTINUE
 1100 CONTINUE
C
C*    2. CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
C        ------------------------------------------------
C
      X0 = 0.05
      DO 2100 L=0,IALPHA
         DO 2200 K=0,IUSTAR
            UST      = MAX(FLOAT(K)*DELUST,0.000001)
            Z0       = UST**2*(ALPHA+FLOAT(L)*DELALP)/G
            OMEGACC  = MAX(OMEGAC,X0*G/UST)
            YC       = OMEGACC*SQRT(Z0/G)
            DELY     = MAX((1.-YC)/FLOAT(JTOT),0.)
            DO 2300 J=1,JTOT
               Y        = YC+FLOAT(J-1)*DELY
               OMEGA    = Y*SQRT(G/Z0)
               CM       = G/OMEGA
               ZX       = UST/CM +ZALP
               ZARG     = MIN(XKAPPA/ZX,20.)
               ZMU      = MIN(G*Z0/CM**2*EXP(ZARG),1.)
C
               ZLOG         = MIN(ALOG(ZMU),0.)
               ZBETA        = CONST1*ZMU*ZLOG**4
               TAUHFT(K,L)  = TAUHFT(K,L)+ZBETA/Y*DELY
 2300       CONTINUE
 2200    CONTINUE
 2100 CONTINUE

      RETURN
      END
