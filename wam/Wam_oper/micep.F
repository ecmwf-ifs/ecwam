C-------------------------------------------------------------------
C
      SUBROUTINE MICEP (IU06, NSTART, NEND, TEMPG)
C
C-------------------------------------------------------------------
C
C**** *MICEP* - DEFINE BLOCK AND INDEX OF EACH ICE POINT.
C
C     R. PORTZ     MPI HAMBURG   OCTOBER 1992
C     J. BIDLOT    ECMWF         JUNE 1996    MESSAGE PASSING
C     P. JANSSEN   ECMWF         OCTOBER 1997 SIMPLIFY BY INTRODUCING
C                                             ICEMASK
C     J. BIDLOT    ECMWF           JANURAY 1998 : CORRECT THRESHOLD FOR ICE
C                                                 TO BE LT 271.50
C
C
C     PURPOSE.
C     --------
C            SELECT POINTS (TEMP < 271.5 KELVIN) IN THE TEMPERATURE
C            GRID AND DEFINE AS SEA ICE POINTS.
C
C
C
C**   INTERFACE
C     ---------
C             *CALL MICEP* *( IU06, NICE, NSTART, NEND)*
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      *IU06*      INTEGER   PRINTER OUTPUT UNIT.
C      *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C      *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C      *TEMPG*     REAL      GLOBAL TEMPRERATURES (KELVIN).
C
C
C     METHOD.
C     -------
C             NONE
C
C
C    EXTERNAL.
C    ---------
C             *NONE* 
C
C---------------------------------------------------------------------
C
#include "param.h"
C
#include "comgrid.h"
C
#include "commap.h"
C
#include "comtemp.h"
C
#include "comice.h"
C
#include "commesspass.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comstat.h"
C
C ----------------------------------------------------------------------
C     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(NXT,NYT) :: TEMPG
C
C ----------------------------------------------------------------------
C
C     1. INITIALIZATION
C     ------------------
C
 1000 CONTINUE
C
C-----------------------------------------------------------------
C
C*    2. WE DEFINE ICE FOR THOSE POINTS WHERE SST < 271.5 KELVIN). 
C        AT THOSE POINTS ICEMASK = 0
C        -------------------------------------------------------
C
 2000 CONTINUE
C
      NICE = 0
      ICEMASK = 1
C
      DO 2100 IG=1,IGL
         IF (LMESSPASS) THEN
C           ICEMASK FOR ALL SEAPOINTS IS NEEDED IN GRDATA
            IF(IASSI.EQ.1) THEN
              IJSTART=NSTART(1)
              IJEND=NEND(NPROC)
            ELSE
              IJSTART=NSTART(IRANK)
              IJEND=NEND(IRANK)
            ENDIF
            DO 2200 IJ = IJSTART, IJEND 
               IX = IXLG(IJ,IG)
               KX = KXLT(IJ,IG)
               XI = AMOWEP + REAL(IX-1)*ZDELLO(KX)-TMOWEP
               XI = MOD(XI+720.,360.)
               XI = XI/DELLO+1.00001
               XJ = TMONOP - AMOSOP - REAL(KX-1)*XDELLA
               XJ = XJ/DELLA+1.00001
C         
               I  = NINT(XI)
               J  = NINT(XJ)
               IF (TEMPG(I,J) .LT. 271.5) THEN
                  NICE = NICE + 1
                  ICEMASK(IJ,IG) = 0
               ENDIF
 2200       CONTINUE
         ELSE
            DO 2300 IJ = IJS(IG),IJL(IG)
               IX = IXLG(IJ,IG)
               KX = KXLT(IJ,IG)
               XI = AMOWEP + REAL(IX-1)*ZDELLO(KX)-TMOWEP
               XI = MOD(XI+720.,360.)
               XI = XI/DELLO+1.00001
               XJ = TMONOP - AMOSOP - REAL(KX-1)*XDELLA
               XJ = XJ/DELLA+1.00001
C         
               I  = NINT(XI)
               J  = NINT(XJ)
               IF (TEMPG(I,J) .LT. 271.5) THEN
                  NICE = NICE + 1
                  ICEMASK(IJ,IG) = 0
               ENDIF
 2300       CONTINUE
         ENDIF
 2100 CONTINUE         
C
      IF(IASSI.EQ.1) THEN
        WRITE(IU06,*)'   POINTS WITH TEMP < 271.5 K:',NICE,
     &             ' ON ALL ',NPROC, 'PEs'
      ELSE
        WRITE(IU06,*)'   POINTS WITH TEMP < 271.5 K:',NICE,
     &             ' ON PE ',IRANK, ' OUT OF A TOTAL OF ',NPROC, 'PEs'
      ENDIF
C
      RETURN
      END
