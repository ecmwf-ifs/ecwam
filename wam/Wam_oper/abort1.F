      SUBROUTINE ABORT1 

! ----------------------------------------------------------------------

!**** *ABORT1* -

!     J. BIDLOT     ECMWF    MAY 1996    MESSAGE PASSING


!*    PURPOSE.
!     --------
!     TO ABORT AFTER CLOSING OPENED UNITS 

!**   INTERFACE.
!     ----------
!     *CALL ABORT1*

!     METHOD.
!     -------
!     CALL MPCLOSE_UNIT THEN CALL ABORT 

!     EXTERNALS.
!     ----------
!       MPABORT
!       ABORT

!     REFERENCE.
!     ----------
!       NONE.

! ----------------------------------------------------------------------

      USE YOWMESPAS, ONLY : LMESSPASS

! ----------------------------------------------------------------------
      IF(LMESSPASS) THEN
        CALL MPABORT
      ELSE
        CALL ABORT
      ENDIF
      RETURN
      END SUBROUTINE ABORT1
