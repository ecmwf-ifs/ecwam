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
!     CALL MPABORT OR ABORT 

!     EXTERNALS.
!     ----------
!       MPABORT
!       ABORT

!     REFERENCE.
!     ----------
!       NONE.

! ----------------------------------------------------------------------

      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      IMPLICIT NONE

      CHARACTER(LEN=30) :: CDMESSAGE

      CALL FLUSH(IU06)

      CDMESSAGE='!! *** WAVE MODEL HAS ABORTED *** !!!!'
      CALL MPABORT(CDMESSAGE)

      RETURN
      END SUBROUTINE ABORT1
