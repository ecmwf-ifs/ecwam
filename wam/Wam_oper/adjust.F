      SUBROUTINE ADJUST(WEST, EAST)

! ----------------------------------------------------------------------

!**** *ADJUST* - ROUTINE TO CORRECT BORDERS OF INTERVALS.

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       ADJUSTS INTERVAL BORDERS GIVEN IN DEGREE.

!**   INTERFACE.
!     ----------

!       *CALL* *ADJUST (WEST, EAST)*
!          *WEST*    - LEFT INTERVAL BORDER IN DEGREE.
!          *EAST*    - RIGHT INTERVAL BORDER IN DEGREE.

!     METHOD.
!     -------

!       THE INTERVAL BORDERS ARE CHANGED TO FULLFILL:
!         0. .LE. EAST  .AND. EAST .LT. 360. .AND. WEST .LE. EAST

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

!* 1. CORRECT BORDERS.
!     ----------------

      WEST = MOD(WEST+720.,360.)
      EAST = MOD(EAST+720.,360.)
      IF (WEST.GT.EAST) WEST = WEST-360.

      RETURN
      END SUBROUTINE ADJUST






