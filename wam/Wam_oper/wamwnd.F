      SUBROUTINE WAMWND (US, DS, KINF, KSUP , KEND, UWND, VWND, NC, NR)

C ----------------------------------------------------------------------
C
C**** *WAMWND* - TRANSFORMS INPUT WINDS TO BLOCKED WAM POINTS.
C
C
C*    PURPOSE.
C     --------
C
C       INTERPOLATE AND CONVERT INPUT WINDS TO WAM WINDS FOR ALL
C       POINTS IN THE GRID.
C
C**   INTERFACE.
C     ----------
C
C       *CALL WAMWND (US, DS, UWND, VWND, NC, NR)*
C          *US*   - INTERPOLATED WINDS IN ALL POINTS AND BLOCKS.
C          *DS*   - INTERPOLATED WIND DIRECTION IN ALL POINTS.
C          *KINF* - LOWER BOUND ON FIRST DIMENSION OF FIELD
C          *KSUP* - UPPER BOUND ON FIRST DIMENSION OF FIELD
C          *KEND* - USED UPPER LIMIT ON FIRST DIMENSION OF FIELD
C          *UWND* - INPUT WIND : U COMPONENT.
C          *VWND* - INPUT WIND : V COMPONENT.
C          *NC*   - DIMENSION OF IN INPUT WIND ARRAYS (COLUMNS).
C          *NR*   - DIMENSION OF IN INPUT WIND ARRAYS (ROWS).
C
C
C     METHOD.
C     -------
C
C       THE INPUT WINDS ARE INTERPOLATED TO THE WAVE MODEL GRID POINTS
C       IN BLOCKED FORMAT. THE INTERPOLATED VALUES ARE TRANSFORMED TO
C       MAGNITUDE AND DIRECTION. INPUT MAY BE WIND IN 10M HEIGHT ,
C       SURFACE WINDS OR FRICTION VELOCITIES. THE INPUT GRID HAS TO BE
C       A LATITUDE/LONGITUDE GRID EITHER PERIODIC OR NON PERIODIC.
C
C     EXTERNALS.
C     ----------
C
C       *LOCINT*    - LOCATES WAM POINT IN INPUT GRID AND INTERPOLATES.
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comgrid.h"
C
#include "commesspass.h"
C
#include "comtest.h"
C
#include "comwndg.h"
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
      PARAMETER (ROAIR = 1.225)
      PARAMETER (XKAPPA = 0.41, ALPHACH = 0.0185)
C
      REAL, DIMENSION (KINF:KSUP,NBLO) :: US,DS
      DIMENSION UWND(NC,NR), VWND(NC,NR)
      LOGICAL LFRSTIME
C
      DATA LFRSTIME /.TRUE./
C
C ----------------------------------------------------------------------
C
C*    1. INITIALISE WIND ARRAYS WITH ZERO
C        --------------------------------
C
 1000 CONTINUE
         US = 0.
         DS = 0.
C
C ----------------------------------------------------------------------
C
C*    2. INTERPOLATE INPUT WINDS TO WAVEMODEL GRID.
C        ------------------------------------------
C
 2000 CONTINUE
      DO 2001 IG=1,IGL
C
C*    2.1 INTERPOLATE ONE BLOCK IN SPACE.
C         ------------------------------
C
        IF(KINF.EQ.1.AND.KSUP.EQ.NIBLO.AND.LMESSPASS) THEN
          CALL LOCINT (IG, 1, KEND, NC, NR, KCOL, KROW, IWPER,
     1         DLAM, DPHI, RLONL, RLATS, UWND, VWND, US(1,IG),
     2         DS(1,IG))
        ELSE
          CALL LOCINT (IG, IJS(IG), IJL(IG), NC, NR, KCOL, KROW, IWPER,
     1         DLAM, DPHI, RLONL, RLATS, UWND, VWND, US(IJS(IG),IG),
     2         DS(IJS(IG),IG))
        ENDIF
 2001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    3. PROCESS WINDS ACCORDING TO TYPE
C        NOTHING TO DO FOR WIND SPEED U10 (ICODE = 3).
C        ---------------------------------------------
C
 3000 CONTINUE
C
C*    3.1 TRANSFORM TO MAGNITUDE AND DIRECTION.
C         -------------------------------------
C
      DO IG = 1,IGL
        IF(KINF.EQ.1.AND.KSUP.EQ.NIBLO.AND.LMESSPASS) THEN
          DO  IJ = 1,KEND
            UU = US(IJ,IG)
            VV = DS(IJ,IG)
            US(IJ,IG) = SQRT(UU**2 + VV**2)
            IF (US(IJ,IG).NE.0.) DS(IJ,IG) = ATAN2(UU,VV)
            IF (DS(IJ,IG).LT.0.) DS(IJ,IG) = DS(IJ,IG) + ZPI
          ENDDO
        ELSE
          DO  IJ = IJS(IG),IJL(IG)
            UU = US(IJ,IG)
            VV = DS(IJ,IG)
            US(IJ,IG) = SQRT(UU**2 + VV**2)
            IF (US(IJ,IG).NE.0.) DS(IJ,IG) = ATAN2(UU,VV)
            IF (DS(IJ,IG).LT.0.) DS(IJ,IG) = DS(IJ,IG) + ZPI
          ENDDO
        ENDIF
      ENDDO 
C
      IF (ICODE.EQ.1) THEN
C
C*    3.2  INPUT IS FRICTION VELOCITY.
C          ---------------------------
C
         DO IG = 1,IGL
           IF(KINF.EQ.1.AND.KSUP.EQ.NIBLO.AND.LMESSPASS) THEN
             DO IJ = 1,KEND
               USTAR = MAX(0.01,US(IJ,IG))
               Z0  = ALPHACH/G*USTAR**2
               CD  = XKAPPA/ALOG(10./Z0)
               US(IJ,IG) = USTAR/CD
             ENDDO
           ELSE
             DO IJ = IJS(IG),IJL(IG)
               USTAR = MAX(0.01,US(IJ,IG))
               Z0  = ALPHACH/G*USTAR**2
               CD  = XKAPPA/ALOG(10./Z0)
               US(IJ,IG) = USTAR/CD
             ENDDO
           ENDIF
         ENDDO
C
      ELSE IF (ICODE.EQ.2) THEN
C
C*    3.3 INPUT WINDS ARE SURFACE STRESSES.
C         ---------------------------------
C
         DO IG = 1,IGL
           IF(KINF.EQ.1.AND.KSUP.EQ.NIBLO.AND.LMESSPASS) THEN
             DO IJ = 1,KEND
               USTAR = MAX (0.01, SQRT(US(IJ,IG)/ROAIR))
               Z0  = ALPHACH/G*USTAR**2
               CD  = XKAPPA/ALOG(10./Z0)
               US(IJ,IG) = USTAR/CD
             ENDDO
           ELSE
             DO IJ = IJS(IG),IJL(IG)
               USTAR = MAX (0.01, SQRT(US(IJ,IG)/ROAIR))
               Z0  = ALPHACH/G*USTAR**2
               CD  = XKAPPA/ALOG(10./Z0)
               US(IJ,IG) = USTAR/CD
             ENDDO
           ENDIF
         ENDDO
      ENDIF
C
C     3.4 SET LOWER THEN 1.0 M/S  WIND SPEED TO 1.0 M/S
C         ---------------------------------------------
      DO IG = 1,IGL
        IF(KINF.EQ.1.AND.KSUP.EQ.NIBLO.AND.LMESSPASS) THEN
          DO IJ = 1,KEND
            US(IJ,IG) = MAX(US(IJ,IG),1.0)
          ENDDO
        ELSE
          DO IJ = IJS(IG),IJL(IG)
            US(IJ,IG) = MAX(US(IJ,IG),1.0)
          ENDDO
        ENDIF
      ENDDO
      IF(LFRSTIME) THEN
      WRITE (IU06,*) ' '
      WRITE (IU06,*) ' ************************************************'
      WRITE (IU06,*) ' SUB. WAMWND : '
      WRITE (IU06,*) ' '
      WRITE (IU06,*) ' WIND SPEEDS LOWER THAN 1 M/S WERE RESET TO 1 M/S'
      WRITE (IU06,*) ' '
      WRITE (IU06,*) ' ************************************************'
      WRITE (IU06,*) ' '
      LFRSTIME=.FALSE.
      ENDIF

C
C*    4. TEST OUTPUT OF WAVE MODEL BLOCKS
C        ---------------------------------
C
 4000 CONTINUE
      IF (ITEST.GE.3) THEN
          WRITE (IU06,*) ' '
          WRITE (IU06,*) '      SUB. WAMWND:',
     1                   ' WINDFIELDS CONVERTED TO BLOCKS'
          WRITE (IU06,*) ' '
          DO 4001 IG = 1,MIN(ITESTB,NBLO)
             IJA = IJS(IG)
             IJE = IJS(IG)+4
             WRITE (IU06,*) ' BLOCK= ',IG
             WRITE (IU06,*) ' '
             WRITE (IU06,*) ' US(IJS - IJS+4): ',(US(I,IG),I=IJA,IJE)
             WRITE (IU06,*) ' DS(IJS - IJS+4): ',(DS(I,IG),I=IJA,IJE)
 4001    CONTINUE
      ENDIF

      RETURN
      END
