      SUBROUTINE ADJUST(WEST, EAST)

C ----------------------------------------------------------------------
C
C**** *ADJUST* - ROUTINE TO CORRECT BORDERS OF INTERVALS.
C
C     H.GUNTHER            ECMWF       04/04/1990
C
C*    PURPOSE.
C     -------
C
C       ADJUSTS INTERVAL BORDERS GIVEN IN DEGREE.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *ADJUST (WEST, EAST)*
C          *WEST*    - LEFT INTERVAL BORDER IN DEGREE.
C          *EAST*    - RIGHT INTERVAL BORDER IN DEGREE.
C
C     METHOD.
C     -------
C
C       THE INTERVAL BORDERS ARE CHANGED TO FULLFILL:
C         0. .LE. EAST  .AND. EAST .LT. 360. .AND. WEST .LE. EAST
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
C* 1. CORRECT BORDERS.
C     ----------------

      WEST = MOD(WEST+720.,360.)
      EAST = MOD(EAST+720.,360.)
      IF (WEST.GT.EAST) WEST = WEST-360.

      RETURN
      END
