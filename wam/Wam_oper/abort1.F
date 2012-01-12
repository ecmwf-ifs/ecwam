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

      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------
      CHARACTER*30 :: CDMESSAGE

      CALL FLUSH(IU06)

      CDMESSAGE='!! WAVE MODEL HAS ABORTED !!!!'
      IF(LMESSPASS) THEN
        CALL MPABORT(CDMESSAGE)
      ELSE
        CALL ABORT
      ENDIF
      RETURN
      END SUBROUTINE ABORT1
