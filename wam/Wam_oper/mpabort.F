      SUBROUTINE MPABORT 

! ----------------------------------------------------------------------

!**** *MPABORT* -

!     J. BIDLOT     ECMWF   2001 


!*    PURPOSE.
!     --------
!     CALL ABORT AFTER CLOSING OPENED UNITS 

!**   INTERFACE.
!     ----------
!     *CALL MPABORT*

!     METHOD.
!     -------
!     CALL MPCLOSE_UNIT THEN CALL ABORT 

!     EXTERNALS.
!     ----------
!       MPCLOSE_UNIT
!       MPL_ABORT

!       ABORT
!     REFERENCE.
!     ----------
!       NONE.

! ----------------------------------------------------------------------

      USE MPL_MODULE

! ----------------------------------------------------------------------
      CALL MPCLOSE_UNIT
      CALL MPL_ABORT('MPL_ABORT was called to terminate run')
      RETURN
      END SUBROUTINE MPABORT
