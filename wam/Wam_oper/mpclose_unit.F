      SUBROUTINE MPCLOSE_UNIT

! ----------------------------------------------------------------------

!**** *MPCLOSE_UNIT* -

!     J. BIDLOT     ECMWF    MAY 1996    MESSAGE PASSING
!                            2007   ONLY USE IT TO FLUSH ALL UNITS.

!*    PURPOSE.
!     --------
!     CLOSE ALL OPENED UNIT

!**   INTERFACE.
!     ----------
!     *CALL MPCLOSE_UNIT*

!     METHOD.
!     -------
!     LOOPS ON ALL UNIT NUMBER AND INQUIRES IF UNIT IS OPENED 
!     AND THEN CLOSES IT

!     EXTERNALS.
!     ----------
!       NONE.
!     REFERENCE.
!     ----------
!       NONE.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IU
      LOGICAL :: LOPENED

      DO IU=1,99
        INQUIRE(UNIT=IU,OPENED=LOPENED)
        IF(LOPENED) THEN 
          CALL FLUSH(IU)
!!!!          CLOSE(IU)
        ENDIF
      ENDDO 
      END SUBROUTINE MPCLOSE_UNIT
