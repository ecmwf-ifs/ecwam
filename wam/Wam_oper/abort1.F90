      SUBROUTINE ABORT1

! ----------------------------------------------------------------------

!**** *ABORT1* -

!     J. BIDLOT     ECMWF    MAY 1996    MESSAGE PASSING


!*    PURPOSE.
!     --------
!     TO ABORT AFTER FLUSHING OPENED UNITS

!**   INTERFACE.
!     ----------
!     *CALL ABORT1*

!     METHOD.
!     -------
!     CALL WAM_ABORT

!     EXTERNALS.
!     ----------
!       WAM_ABORT

!     REFERENCE.
!     ----------
!       NONE.

! ----------------------------------------------------------------------

      USE YOWTEST, ONLY: IU06
      USE YOWABORT, ONLY : WAM_ABORT
! ----------------------------------------------------------------------

      IMPLICIT NONE

      CALL FLUSH (IU06)
      CALL WAM_ABORT()

      END SUBROUTINE ABORT1
