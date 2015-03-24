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

      IMPLICIT NONE

      REAL, INTENT(INOUT) :: WEST, EAST 
      REAL :: OLD_WEST, OLD_EAST

!* 1. CORRECT BORDERS.
!     ----------------

      OLD_WEST=WEST
      OLD_EAST=EAST
      WEST = MOD(WEST+720.,360.)
      EAST = MOD(EAST+720.,360.)
      IF(OLD_WEST.NE.OLD_EAST) THEN
        IF (WEST.GE.EAST) WEST = WEST-360.
      ENDIF

      RETURN
      END SUBROUTINE ADJUST






