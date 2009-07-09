      SUBROUTINE MPABORT(CDMESSAGE) 

! ----------------------------------------------------------------------

!**** *MPABORT* -

!     J. BIDLOT     ECMWF   2001 


!*    PURPOSE.
!     --------
!     CALL ABORT AFTER CLOSING OPENED UNITS 

!**   INTERFACE.
!     ----------
!     *CALL MPABORT(CDMESSAGE)*

!     METHOD.
!     -------
!     CALL MPCLOSE_UNIT THEN CALL ABORT 

!     EXTERNALS.
!     ----------
!       MPCLOSE_UNIT
!       ABOR1

!       ABORT
!     REFERENCE.
!     ----------
!       NONE.

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU  
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE MPL_MODULE
      USE SDL_MOD

! ----------------------------------------------------------------------
      CHARACTER*(*),INTENT(IN) :: CDMESSAGE

      LOGICAL :: LLABORT=.TRUE.

!     FLUSH UNITS
      CALL MPCLOSE_UNIT

      IF(LWCOU) THEN
        CALL ABOR1(CDMESSAGE)
      ELSEIF (NPROC > 1) THEN
        CALL MPL_MESSAGE(CDMESSAGE=CDMESSAGE,LDABORT=LLABORT)
      ELSE
        CALL SDL_TRACEBACK
        CALL FLUSH(0)
        CALL SDL_SRLABORT
      ENDIF

      RETURN
      END SUBROUTINE MPABORT
