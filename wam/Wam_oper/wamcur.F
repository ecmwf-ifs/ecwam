      SUBROUTINE WAMCUR (IUCUR, NCC, NRC)

C ----------------------------------------------------------------------
C
C**** *WAMCUR* - ROUTINE PROCESS A CURRENT FIELD.
C
C     SUSANNE HASSELMANN        JUNE     1990.
C     H. GUNTHER    ECMWF/GKSS  DECEMBER 1990  MODIFIED FOR CYCLE_4.
C
C*    PURPOSE.
C     --------
C
C       INPUT AND INTERPOLATION TO BLOCKED CURRENT ARRAYS.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *WAMCUR (IUCUR, NCC, NRC)*
C          *IUCUR*  - INPUT UNIT FOR CURRENTS.
C          *NCC*    - NUMBER OF COLUMNES IN CURRENT GRID (DIMENSION)
C          *NRC*    - NUMBER OF ROWS     IN CURRENT GRID (DIMENSION)
C
C     METHOD.
C     -------
C
C       LINEAR INTERPOLATION.
C
C     EXTERNALS.
C     ----------
C
C       *LOCINT*    - INTERPOLATE TO BLOCKS.
C       *READCUR*   - READ A CURRENT FIELD.
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
#include "comcurr.h"
C
#include "comcurg.h"
C
#include "comgrid.h"
C
#include "commap.h"
C
#include "comtest.h"
C
C ----------------------------------------------------------------------
C
      DIMENSION UCUR(NCC,NRC), VCUR(NCC,NRC)
      CHARACTER*12 CDTCUR
C
C ----------------------------------------------------------------------
C
C*    1. READ CURRENT FILE.
C        ------------------
C
 1000 CONTINUE
      CALL READCUR (IUCUR, CDTCR, UCUR, VCUR, NCC, NRC)
C
C ----------------------------------------------------------------------
C
C*    2. INTERPOLATION OF INPUT CURRENT ON MODEL GRID.
C        ---------------------------------------------
C
 2000 CONTINUE
      IJA =1
      DO 2001 IG=1,IGL
         CALL LOCINT (IG, IJA, IJLT(IG), NCC, NRC, KCCI, KRCI, IPCUR,
     1          DUCLO, DUCLA, UCLW, UCLO, UCUR, VCUR, U(IJS(IG),IG),
     2          V(IJS(IG),IG))
 2001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    3. TEST OUTPUT OF WAVE MODEL BLOCKS
C        ---------------------------------
C
 3000 CONTINUE
      IF (ITEST.GE.3) THEN
          WRITE (IU06,*) ' '
          WRITE (IU06,*) '      SUB. WAMCUR..',
     1                   ' CURRENTS CONVERTED TO BLOCKS'
          WRITE (IU06,*) ' '
          DO 3001 IG = 1,ITESTB
             IJA = IJS(IG)
             IJE = IJS(IG)+4
             WRITE (IU06,*) ' BLOCK= ',IG
             WRITE (IU06,*) ' '
             WRITE (IU06,*) ' U(IJS - IJS+4): ',(U(I,IG),I=IJA,IJE)
             WRITE (IU06,*) ' V(IJS - IJS+4): ',(V(I,IG),I=IJA,IJE)
 3001    CONTINUE
      ENDIF

      RETURN
      END
