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
      USE MPL_MODULE
      USE SDL_MOD

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "mpclose_unit.intfb.h"

      CHARACTER(LEN=*),INTENT(IN) :: CDMESSAGE

      LOGICAL :: LLABORT=.TRUE.

!     FLUSH UNITS
      CALL MPCLOSE_UNIT

      IF (LWCOU) THEN
        CALL ABOR1(CDMESSAGE)
        STOP 2
      ELSE
        STOP 1
      ENDIF

      END SUBROUTINE MPABORT
