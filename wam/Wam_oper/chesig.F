C ----------------------------------------------------------------------
C
      SUBROUTINE CHESIG (KU06, KTST, KRANK, KPROC, LDSIGSTOP, LDSIGREST)
C
C ----------------------------------------------------------------------
C
C**** *CHESIG* - CHECK SIGNALS.
C
C     PURPOSE.
C     --------
C       CHECK FOR EXTERNAL SIGNALS REQUESTING CREATION OF RESTART FILES.
C
C**   INTERFACE.
C     ----------
C        *CALL* *CHESIG(KU06,KRANK,KPROC,LDSIGSTOP,LDSIGREST)
C
C    I/   *KU06*      - LOGICAL I/O UNIT FOR STDOUT.
C    I/   *KTST*      - TEST OUTPUT LEVEL.
C    I/   *KRANK*     - CURRENT PROCESSOR ELEMENT (NUMERIC).
C    I/   *KPROC*     - TOTAL NUMBER OF PROCESSING ELEMENTS.
C     /O  *LDSIGSTOP* - SET .TRUE. IF STOP SIGNAL RECEIVED.
C     /O  *LDSIGREST* - SET .TRUE. IF RESTART FILE SIGNAL RECEIVED.
C
C     METHOD.
C     -------
C       SIGNALS ARE UNBLOCKED.
C       MASTER PROCESS CHECKS FLAGS.
C       THEY MAY BE SET BY SIGNAL PROCESSOR (see IFSSIG)
C       VALUES OF FLAGS ARE SENT TO ALL OTHER PROCESSES.
C
C     EXTERNALS.  
C     ----------
C       IFSSIGB
C       MPE_SEND
C       MPE_RECV
C       MPE_BARRIER
C       FLUSH
C
C     REFERENCE.
C     ----------
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS.
C
C     AUTHOR.
C     -------
C        BJORN HANSEN *ECMWF* MARCH 1997
C
C     Modifications.
C     --------------
C        ORIGINAL : 96-07-29 BY DAVID DENT.
C
C     ------------------------------------------------------------------
C
      IMPLICIT LOGICAL (L)
C
#include "comstat.h"
C
      INTEGER, PARAMETER :: JPSIGTAG=1000
      INTEGER, PARAMETER :: JPICODE=1
C
C ----------------------------------------------------------------------
C
C*    1. UNBLOCK SIGNALS.
C        ----------------
C
      IF (KTST.GE.2)
     . WRITE(KU06,*) '   CHESIG:  FOR PE ', KRANK, ' START'
C
      CALL IFSSIGB
C
C         check for arrival of signals
C
      LDSIGSTOP=.FALSE.
      LDSIGREST=.FALSE.
C
      IF (KTST.GE.2) THEN
      WRITE(KU06,*) '   CHESIG:  FOR PE ', KRANK, ' STOP    = ', ISIGHUP
      WRITE(KU06,*) '   CHESIG:  FOR PE ', KRANK, ' RESTART = ', ISIGINT
      ENDIF
C
      IF (KRANK .EQ. 1) THEN
        DO JE=2, KPROC
          ITAG=JPSIGTAG+JE
          CALL MPE_SEND(ISIGHUP, 2, JPICODE, JE, ITAG, 0, 0, 0, IERR)
          IF (KTST.GE.2)
     .    WRITE(KU06,*)' SEND STATUS TO PE ', JE, ' R-CODE= ', IERR
        ENDDO
      ELSE
        ITAG=JPSIGTAG+KRANK
        IPROC=1
        CALL MPE_RECV(ISIGHUP,2,JPICODE,IPROC,ITAG,
     .                0,0,0,ILREC,IRECD,ITAG,IERR)
        IF (KTST.GE.2)
     .  WRITE(KU06,*)' RECEIVED STATUS ON PE ',KRANK, 'R-CODE= ', IERR,
     .               'STOP= ',ISIGHUP, 'RESTART= ',ISIGINT
      ENDIF
      IF (ISIGHUP.NE.0) LDSIGSTOP=.TRUE.
      IF (ISIGINT.NE.0) LDSIGREST=.TRUE.
      ISIGINT=0
      CALL MPE_BARRIER(IERR)
      IF(IERR.LT.0) CALL MPEI_ABORT('MPE_BARRIER ERROR in CHESIG')

      IF (KTST.GE.2) THEN
        WRITE(KU06,*) '   CHESIG:  FOR PE ', KRANK, ' DONE'
        WRITE(KU06,*) '   CHESIG:  LDSIGSTOP:', LDSIGSTOP,
     .                ' LDSIGREST: ', LDSIGREST
        CALL FLUSH(KU06)
      ENDIF
C
C     ------------------------------------------------------------------
C
      RETURN
      END

