      SUBROUTINE ABORT1 

! ----------------------------------------------------------------------

!**** *ABORT1* -

!     J. BIDLOT     ECMWF    MAY 1996    MESSAGE PASSING


!*    PURPOSE.
!     --------
!     CALL ABORT AFTER CLOSING OPENED UNITS 

!**   INTERFACE.
!     ----------
!     *CALL ABORT1*

!     METHOD.
!     -------
!     CALL MPCLOSE_UNIT THEN CALL ABORT 

!     EXTERNALS.
!     ----------
!       MPCLOSE_UNIT
!       MPEI_ABORT

!       ABORT
!     REFERENCE.
!     ----------
!       NONE.

! ----------------------------------------------------------------------

      USE YOWMESPAS, ONLY : LMESSPASS

! ----------------------------------------------------------------------
      IF(LMESSPASS) THEN
        CALL MPCLOSE_UNIT
        CALL MPEI_ABORT('MPEI_ABORT was called to terminate run')
      ELSE
        CALL ABORT
      ENDIF
      RETURN
      END SUBROUTINE ABORT1
