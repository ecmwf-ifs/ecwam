      SUBROUTINE LOCINT (IG, IJS, IJL, NC, NR, KCOL, KROW, IGPER,
     1                  DLAM, DPHI, RLONL, RLATS, U, V, US, DS)

C ----------------------------------------------------------------------
C
C**** *LOCINT* - LOCATE WAM POINT IN INPUT GRID AND INTERPOLATE.
C
C     H. GUNTHER    ECMWF/GKSS  DECEMBER 1990  MODIFIED FOR CYCLE_4.
C
C*    PURPOSE.
C     --------
C
C        *LOCINT* - LOCATE AND INTERPOLATE IN INPUT GRID.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *LOCINT (IG, IJS, IJL, NC, NR, KCOL, KROW, IGPER,
C                       DLAM, DPHI, RLONL, RLATS, U, V, US, DS)*
C         *IG*      - BLOCK NUMBER.
C         *IJS*     - FIRST POINT IN BLOCK.
C         *IJL*     - LAST POINT IN BLOCK.
C         *NC*      - DIMENSION OF INPUT ARRAY, NUMBER OF COLUMNES.
C         *NR*      - DIMENSION OF INPUT ARRAY, NUMBER OF ROWS.
C         *KCOL*    - NUMBER OF COLUMNES IN INPUT (USED).
C         *KROW*    - NUMBER OF ROWS     IN INPUT (USED).
C         *IGPER*   - INDICATOR PERIODICAL INPUT GRID OR NOT
C                     1 = PERIODICAL
C                     OTHERWISE NON-PERIODICAL
C         *DLAM*    - STEPSIZE BETWEEN LONGITUDES IN INPUT (DEG).
C         *DPHI*    - STEPSIZE BETWEEN LATITUDES  IN INPUT (DEG).
C         *RLATS*   - MOST SOUTHERN LATITUDE OF INPUT (DEG).
C         *RLONL*   - MOST WESTERN LONGITUDE OF INPUT (DEG).
C         *U(I,J)*  - "HORIZONTAL" COMPONENT (INPUT).
C         *V(I,J)*  - "VERTICAL  " COMPONENT (INPUT).
C         *US*      - SPACE INTERPOLATED "HORIZONTAL" COMPONENT.
C         *DS*      - SPACE INTERPOLATED "VERTICAL  " COMPONENT.
C
C     METHOD.
C     -------
C
C       DOUBLE LINEAR INTERPOLATION WITHIN A MESH FOR U AND V.
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
#include "commap.h"
C
#include "comtest.h"
C
C ----------------------------------------------------------------------
C
!
      REAL,DIMENSION(IJS:IJL)  :: US,DS

      DIMENSION U(NC,NR), V(NC,NR)
C
C ----------------------------------------------------------------------
C
C*    1. LOOP OVER POINTS IN WAMODEL BLOCKS.
C        -----------------------------------
C
 1000 CONTINUE
      DO 1001 IJ=IJS,IJL
C
C*    1.1 TRANSFORM WAM COORDINATE TO INPUT GRID.
C         ---------------------------------------
C
         KX = KXLT(IJ,IG)
         XI = AMOWEP + REAL(IXLG(IJ,IG)-1)*ZDELLO(KX)-RLONL
         XI = MOD(XI+720.,360.)
         XI = XI/DLAM+1.00001
         XJ = AMOSOP + REAL(KXLT(IJ,IG)-1)*XDELLA-RLATS
         XJ = XJ/DPHI+1.00001
C
C*    1.2 COMPUTE CORNER POINT INDICES OF INPUT GRID.
C         -------------------------------------------
C
         I1  = XI
         J1  = XJ
         J2  = MIN(KROW,J1+1)
         I2  = I1+1
C
C*    1.3 DISTANCES OF INTERPOLATION POINT FROM CORNER POINTS.
C         ----------------------------------------------------
C
         DI1 = XI-I1
         DI2 = 1.-DI1
         DJ1 = XJ-J1
         DJ2 = 1.-DJ1
C
C*    1.4. CORRECTIONS FOR FIRST AND LAST GRID LINES
C          PERIODIC OR UNPERIODIC GRIDS ARE CONSIDERED.
C          --------------------------------------------
C
         IF (IGPER.EQ.1) THEN
            IF (I1.EQ.KCOL) I2 = 1
            IF (I1.EQ.0   ) I1 = KCOL
         ELSE
            IF (I1.EQ.KCOL) I2 = KCOL
         ENDIF
C
C*    1.5 CHECK WHETHER POINTS ARE IN GRID.
C         ---------------------------------
C
         IF (I1.LT.1 .OR. I1.GT.KCOL .OR. J1.LT.1 .OR. J1.GT.KROW) THEN
            WRITE(IU06,*) ' *******************************************'
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' *          FATAL ERROR IN LOCINT          *'
            WRITE(IU06,*) ' *          =====================          *'
            WRITE(IU06,*) ' * POINT IS OUTSIDE OF INPUT GRID          *'
            WRITE(IU06,*) ' * DIMENSION OF INPUT GRID IS KCOL = ', KCOL
            WRITE(IU06,*) ' *                            KROW = ', KROW
            WRITE(IU06,*) ' * POINT REQUESTED IS         XI   = ', XI
            WRITE(IU06,*) ' *                            XJ   = ', XJ
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS      *'
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' *******************************************'
            CALL ABORT1
         ENDIF
C
C*    1.6 LINEAR INTERPOLATION.
C         ---------------------
C
         U1 = U(I1,J1)*DI2 + U(I2,J1)*DI1
         U2 = U(I1,J2)*DI2 + U(I2,J2)*DI1
         US(IJ) = U1*DJ2 + U2*DJ1
         V1 = V(I1,J1)*DI2 + V(I2,J1)*DI1
         V2 = V(I1,J2)*DI2 + V(I2,J2)*DI1
         DS(IJ) = V1*DJ2 + V2*DJ1
 1001 CONTINUE

      RETURN
      END
