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

      INTEGER :: IERROR
      LOGICAL :: LLABORT=.TRUE.

      IERROR=-1

!     FLUSH UNITS
      CALL MPCLOSE_UNIT


      IF(LWCOU) THEN
!       USE THE SAME ABORT AS THE IFS.
!        IF(PRESENT(CDMESSAGE)) THEN
          CALL ABOR1(CDMESSAGE)
!        ELSE
!          CALL ABOR1('!!!! THE WAVE MODEL HAS ABORTED !!!!')
!        ENDIF
      ELSEIF (NPROC > 1) THEN
!        IF(PRESENT(CDMESSAGE)) THEN
          CALL MPL_MESSAGE(KERROR=IERROR,CDMESSAGE=CDMESSAGE,
     &                     LDABORT=LLABORT)
!        ELSE
!          CALL MPL_MESSAGE(KERROR=IERROR,
!     &                     CDMESSAGE='!!! THE MODEL HAS ABORTED !!!',
!     &                     LDABORT=LLABORT)
!        ENDIF
      ELSE
        CALL SDL_TRACEBACK
        CALL FLUSH(0)
        CALL SDL_SRLABORT
      ENDIF

      RETURN
      END SUBROUTINE MPABORT
