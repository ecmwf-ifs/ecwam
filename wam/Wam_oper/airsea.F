      SUBROUTINE AIRSEA (U10, TAUW, US, Z0, IJS, IJL)

C ----------------------------------------------------------------------
C
C**** *AIRSEA* - DETERMINE TOTAL STRESS IN SURFACE LAYER.
C
C     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
C
C*    PURPOSE.
C     --------
C
C       COMPUTE TOTAL STRESS.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *AIRSEA (U10, TAUW, US, Z0, IJS, IJL)*
C          *U10*  - INPUT BLOCK OF WINDSPEEDS U10.
C          *TAUW* - INPUT BLOCK OF WAVE STRESSES.
C          *US*   - OUTPUT BLOCK OF SURFACE STRESSES.
C          *ZO*   - OUTPUT BLOCK OF ROUGHNESS LENGTH.
C          *IJS*  - INDEX OF FIRST GRIDPOINT.
C          *IJL*  - INDEX OF LAST GRIDPOINT.
C
C     METHOD.
C     -------
C
C       USE TABLE TAUT(TAUW,U) AND LINEAR INTERPOLATION.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCE.
C     ---------
C
C       NONE.
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
!
      REAL,DIMENSION(IJS:IJL) ::  U10,TAUW,US,Z0
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
C ----------------------------------------------------------------------
C
C*    1. SELECT TABLE ACCORDING TO WIND LEVEL.
C        -------------------------------------
C
      IF (LWCOU) THEN
        ILEV=2
      ELSE
        ILEV=1
      ENDIF
C
C
C*    2. DETERMINE TOTAL STRESS FROM TABLE.
C        ----------------------------------
C
 2000 CONTINUE
      DO 2001 IJ=IJS,IJL
         XI      = TAUW(IJ)/DELTAUW
         I       = MIN ( ITAUMAX-1, INT(XI) )
         DELI1   = XI - FLOAT(I)
         DELI2   = 1. - DELI1
         XJ      = U10(IJ)/DELU
         J       = MIN ( JUMAX-1, INT(XJ) )
         DELJ1   = XJ - FLOAT(J)
         DELJ2   = 1. - DELJ1
         US(IJ)  = (SQRT(TAUT(I,J,ILEV))*DELI2
     1              + SQRT(TAUT(I+1,J,ILEV))*DELI1)*DELJ2
     2              +(SQRT(TAUT(I,J+1,ILEV))*DELI2
     3              + SQRT(TAUT(I+1,J+1,ILEV))*DELI1)*DELJ1
 2001 CONTINUE
C
C*    3. DETERMINE ROUGHNESS LENGTH.
C        ---------------------------
C
 3000 CONTINUE
      DO 3001 IJ=IJS,IJL
         SQRTCD  = MIN(U10(IJ)/US(IJ),100.0)
         Z0(IJ)  = XNLEV(ILEV)*EXP(-XKAPPA*SQRTCD)
 3001 CONTINUE
C
      RETURN
      END
