      SUBROUTINE ABORT1 

C ----------------------------------------------------------------------
C
C**** *ABORT1* -
C
C     J. BIDLOT     ECMWF    MAY 1996    MESSAGE PASSING
C
C
C*    PURPOSE.
C     --------
C     CALL ABORT AFTER CLOSING OPENED UNITS 
C
C**   INTERFACE.
C     ----------
C     *CALL ABORT1*
C
C     METHOD.
C     -------
C     CALL MPCLOSE_UNIT THEN CALL ABORT 
C
C     EXTERNALS.
C     ----------
C       MPCLOSE_UNIT
C       MPEI_ABORT
C
C       ABORT
C     REFERENCE.
C     ----------
C       NONE.
C
C ----------------------------------------------------------------------
C
#include "commesspass.h"
C
C ----------------------------------------------------------------------


      IF(LMESSPASS) THEN
         CALL MPCLOSE_UNIT
         CALL MPEI_ABORT('MPEI_ABORT was called to terminate run')
      ELSE
         CALL ABORT
      ENDIF

      RETURN
      END
