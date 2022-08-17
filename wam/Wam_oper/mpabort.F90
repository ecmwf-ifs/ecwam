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
!       WAM_ABORT

!     REFERENCE.
!     ----------
!       NONE.

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU  
      USE YOWABORT, ONLY : WAM_ABORT

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "mpclose_unit.intfb.h"

      CHARACTER(LEN=*),INTENT(IN) :: CDMESSAGE

      LOGICAL :: LLABORT=.TRUE.

      IF (LWCOU) THEN
        CALL WAM_ABORT(CDMESSAGE)
        STOP 2
      ELSE
        CALL MPCLOSE_UNIT
        STOP 1
      ENDIF

      END SUBROUTINE MPABORT
