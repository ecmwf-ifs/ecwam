      SUBROUTINE HEADBC (NBOUNC, IDELPRO, TH0, FR1, IU19, IU06)

C ----------------------------------------------------------------------
C
C**** *HEADBC* - OUTPUT OF THE COARSE GRID BOUNDARY VALUE FILE HEADER.
C
C     R. PORTZ     MPI          JANUARY 1991
C     J. BIDLOT    ECMWF        JANUARY 1997 : TH0 IS PASSED AS AN ARGUMENT 
C
C*    PURPOSE.
C     --------
C
C       WRITE A HEADER TO THE BOUNDARY VALUE OUTPUT FILE.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *HEADBC (NBOUNC, IDELPRO, TH0, FR1, IU19, IU06)*
C          *NBOUNC*  - NUMBER OF OUTPUT POINTS.
C          *IDELPRO* - PROPAGATION = OUTPUT TIMESTEP (SECONDS).
C          *TH0*     - FIRST DIRECTION (IN RADIAN)
C          *FR1*     - FIRST FREQUENCY (HERTZ).
C          *IU19*    - OUTPUT UNIT OF BOUNDARY VALUES.
C          *IU06*    - PRINTER OUTPUT UNIT.
C
C     METHOD.
C     -------
C
C       SEQUENTIAL UNFORMATED WRITE TO UNIT.
C
C     EXTERNALS.
C     ----------
C
C       *ABORT1*     - TERMINATES PROCESSING.
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
C ----------------------------------------------------------------------
C
      PARAMETER (CO = 1.1)
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C        *CO*      REAL      FREQUENCY RATIO.
C
C ----------------------------------------------------------------------
C
C*    1. FORMAT AND WRITE HEADER.
C        ------------------------
C
 1000 CONTINUE
      XANG = REAL(NANG)
      XFRE = REAL(NFRE)
      XBOU = REAL(NBOUNC)
      XDEL = REAL(IDELPRO)
      WRITE(IU19,ERR=2000) XANG, XFRE, TH0, FR1, CO, XBOU, XDEL

      RETURN
C
C ----------------------------------------------------------------------
C
C*    2. ERROR MESSAGE.
C        --------------
C
 2000 CONTINUE
      WRITE(IU06,*) '****************************************'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*    FATAL ERROR IN SUB. HEADBC        *'
      WRITE(IU06,*) '*    ===========================       *'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*    WRITE ERROR FROM UNIT : IU19 = ', IU19
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*        PROGRAM ABORTS                *'
      WRITE(IU06,*) '****************************************'
      CALL ABORT1

      END
