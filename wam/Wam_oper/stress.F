      SUBROUTINE STRESS(IU06)

C ----------------------------------------------------------------------
C
C**** *STRESS* - COMPUTATION OF TOTAL STRESS.
C
C     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
C     J. BIDLOT           ECMWF     SEPTEMBER 1996 : REMOVE Z0 DUE TO 
C                                   VISCOSITY AND ADD ZERO STRESS FOR
C                                   ZERO WIND.     
C     BJORN HANSEN        ECMWF     MAY 1997
C                                   STRESS FOR MORE THAN ONE LEVEL.
C
C*    PURPOSE.
C     ---------
C
C       TO GENERATE STRESS TABLE TAU(TAUW,U10).
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *STRESS(IU06)*
C          *IU06*  -  LOGICAL UNIT FOR PRINTER OUTPUT UNIT.
C
C     METHOD.
C     -------
C
C       A STEADY STATE WIND PROFILE IS ASSUMED.
C       THE WIND STRESS IS COMPUTED USING THE ROUGHNESSLENGTH
C
C                  Z1=Z0/SQRT(1-TAUW/TAU)
C
C       WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-
C       INDUCED STRESS AND TAU IS THE TOTAL STRESS.
C       WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH TAUW/TAU < 1.
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
#include "comcoup.h"
C
#include "comtabl.h"
C
C ----------------------------------------------------------------------
C
      PARAMETER (XM=0.50, XNU=0.00001, G=9.806, NITER=10, EPS1=0.00001)
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      *XM*        REAL      POWER OF TAUW/TAU IN ROUGHNESS LENGTH.
C      *XNU*       REAL      KINEMATIC VISCOSITY OF AIR.
C      *G*         REAL      ACCELERATION OF GRAVITY.
C      *NITER*     INTEGER   NUMBER OF ITERATIONS TO OBTAIN TOTAL STRESS
C      *EPS1*      REAL      SMALL NUMBER TO MAKE SURE THAT A SOLUTION
C                            IS OBTAINED IN ITERATION WITH TAU>TAUW.
C
C ----------------------------------------------------------------------
C
C*    1.DETERMINE TOTAL STRESS.
C       -----------------------
C
C*    1.1 INITIALISE CONSTANTS.
C         ---------------------
C
      UMAX    = 50.
      TAUWMAX = 5.
      DELU    = UMAX/FLOAT(JUMAX)
      DELTAUW = TAUWMAX/FLOAT(ITAUMAX)
C
C*    1.2 DETERMINE STRESS.
C         -----------------
C
      DO 1300 JL=1,MIN(JPLEVT,JPLEVC)
         XL=XNLEV(JL)
         WRITE(IU06,*)' STRESS FOR LEVEL ',XL
         DO 1200 I=0,ITAUMAX
            DO 1100 J=0,JUMAX
               ZTAUW   = FLOAT(I)*DELTAUW
               UTOP    = FLOAT(J)*DELU
               CDRAG   = 0.0012875
               WCD     = SQRT(CDRAG)
               USTOLD  = UTOP*WCD
               TAUOLD  = MAX(USTOLD**2, ZTAUW+EPS1)
               DO 1000 ITER=1,NITER
                  X      = ZTAUW/TAUOLD
                  UST    = SQRT(TAUOLD)
                  Z0     = ALPHA*UST**2/(G)/(1.-X)**XM
ccc               ZNU    = 0.1*XNU/UST
ccc               Z0     = MAX(ZNU,Z0)
                  F      = UST-XKAPPA*UTOP/(ALOG(XL/Z0))
                  DELF   = 1.-XKAPPA*UTOP/(ALOG(XL/Z0))**2*2./UST*
     *                     (1.-(XM+1)*X)/(1.-X)
                  UST    = UST-F/DELF
                  TAUOLD =  MAX(UST**2., ZTAUW+EPS1)
 1000          CONTINUE
               TAUT(I,J,JL)  = TAUOLD
 1100       CONTINUE
C
C*    END DO LOOP OVER INDICES OF TAU-TABLE
C
 1200    CONTINUE
C
C*    END DO LOOP OVER LEVELS.
C
 1300 CONTINUE
C
C*    FORCE ZERO WIND TO HAVE ZERO STRESS
C
      DO 2200 JL=1,JPLEVT
        DO 2000 I=0,ITAUMAX
           TAUT(I,0,JL)=0.0
2000    CONTINUE
2200  CONTINUE
C
      RETURN
      END
