      SUBROUTINE MAKEGRID (BLOCK, GRID, IJS, IJL, IG)

C ----------------------------------------------------------------------
C
C**** *MAKEGRID* - MAKE GRIDDED WAM MODEL FIELDS FROM BLOCKED FIELDS.
C
C     H. GUNTHER       ECMWF    NOVEMBER 1989
C
C*    PURPOSE.
C     --------
C
C       GRIDDED WAVE FIELDS ARE CREATED FROM BLOCKED WAVE FIELDS.
C
C**   INTERFACE.
C     ----------
C
C        *CALL MAKEGRID*
C           *BLOCK*   REAL   DATA IN BLOCKED FORMAT
C           *GRID*    REAL   DATA IN GRID FORMAT
C           *IJS*     INTEGER INDEX OF FIRST GRIDPOINT IN BLOCK
C           *IJL*     INTEGER INDEX OF LAST  GRIDPOINT IN BLOCK
C           *IG*      INTEGER BLOCK NUMBER
C
C     METHOD.
C     -------
C
C       THE PARAMETER, WHICH IS GIVEN IN BLOCKED FORMAT, IS
C       DISTRIBUTED IN GRID FORMAT. ONE BLOCK IS DONE IN ONE
C       CALL. BEFORE THE FIRST BLOCK IS TRANSFORMED THE GRID ARRAY IS
C       INITIALISED WITH ZMISS.
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
C
#include "param.h"
C
#include "parcons.h"
C
#include "commap.h"
C
C ----------------------------------------------------------------------
C
      DIMENSION BLOCK(NIBLO), GRID(NGX,NGY)
C
C ----------------------------------------------------------------------
C
C*    1. INITIALIZATION
C     -----------------
C
      IF (IG.EQ.1) THEN
         DO J = 1,NGY
           DO I = 1,NGX
             GRID(I,J) = ZMISS
           ENDDO
         ENDDO
      ENDIF
C
C*    2. MAKE GRIDDED FIELD
C     ---------------------
C
      DO 2000 IJ = IJS,IJL
         IX = IXLG(IJ,IG)
         IY = NGY- KXLT(IJ,IG) +1
         GRID(IX,IY) = BLOCK(IJ)
 2000 CONTINUE

      RETURN
      END
