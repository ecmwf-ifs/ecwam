! ----------------------------------------------------------------------

      SUBROUTINE CHESIG (KU06, KTST, KRANK, KPROC, LDSIGSTOP, LDSIGREST)

! ----------------------------------------------------------------------

!**** *CHESIG* - CHECK SIGNALS.

!     PURPOSE.
!     --------
!       CHECK FOR EXTERNAL SIGNALS REQUESTING CREATION OF RESTART FILES.

!**   INTERFACE.
!     ----------
!        *CALL* *CHESIG(KU06,KRANK,KPROC,LDSIGSTOP,LDSIGREST)

!    I/   *KU06*      - LOGICAL I/O UNIT FOR STDOUT.
!    I/   *KTST*      - TEST OUTPUT LEVEL.
!    I/   *KRANK*     - CURRENT PROCESSOR ELEMENT (NUMERIC).
!    I/   *KPROC*     - TOTAL NUMBER OF PROCESSING ELEMENTS.
!     /O  *LDSIGSTOP* - SET .TRUE. IF STOP SIGNAL RECEIVED.
!     /O  *LDSIGREST* - SET .TRUE. IF RESTART FILE SIGNAL RECEIVED.

!     METHOD.
!     -------
!       SIGNALS ARE UNBLOCKED.
!       MASTER PROCESS CHECKS FLAGS.
!       THEY MAY BE SET BY SIGNAL PROCESSOR (see IFSSIG)
!       VALUES OF FLAGS ARE SENT TO ALL OTHER PROCESSES.

!     EXTERNALS.  
!     ----------
!       IFSSIGB
!       MPE_SEND
!       MPE_RECV
!       MPE_BARRIER
!       FLUSH

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS.

!     AUTHOR.
!     -------
!        BJORN HANSEN *ECMWF* MARCH 1997

!     Modifications.
!     --------------
!        ORIGINAL : 96-07-29 BY DAVID DENT.

!     ------------------------------------------------------------------


      USE YOWSTAT  , ONLY : ISIGHUP  ,ISIGINT


      IMPLICIT LOGICAL (L)

      INTEGER, PARAMETER :: JPSIGTAG=1000
      INTEGER, PARAMETER :: JPICODE=1

! ----------------------------------------------------------------------

!*    1. UNBLOCK SIGNALS.
!        ----------------

      IF (KTST.GE.2)
     & WRITE(KU06,*) '   CHESIG:  FOR PE ', KRANK, ' START'

      CALL IFSSIGB

!         check for arrival of signals

      LDSIGSTOP=.FALSE.
      LDSIGREST=.FALSE.

      IF (KTST.GE.2) THEN
        WRITE(KU06,*) '   CHESIG:  FOR PE ', KRANK, ' STOP    = ', 
     &   ISIGHUP
        WRITE(KU06,*) '   CHESIG:  FOR PE ', KRANK, ' RESTART = ', 
     &   ISIGINT
      ENDIF

      IF (KRANK .EQ. 1) THEN
        DO JE=2, KPROC
          ITAG=JPSIGTAG+JE
          CALL MPE_SEND(ISIGHUP, 2, JPICODE, JE, ITAG, 0, 0, 0, IERR)
          IF (KTST.GE.2)
     &    WRITE(KU06,*)' SEND STATUS TO PE ', JE, ' R-CODE= ', IERR
        ENDDO
      ELSE
        ITAG=JPSIGTAG+KRANK
        IPROC=1
        CALL MPE_RECV(ISIGHUP,2,JPICODE,IPROC,ITAG,
     &                0,0,0,ILREC,IRECD,ITAG,IERR)
        IF (KTST.GE.2)
     &  WRITE(KU06,*)' RECEIVED STATUS ON PE ',KRANK, 'R-CODE= ', IERR,
     &               'STOP= ',ISIGHUP, 'RESTART= ',ISIGINT
      ENDIF
      IF (ISIGHUP.NE.0) LDSIGSTOP=.TRUE.
      IF (ISIGINT.NE.0) LDSIGREST=.TRUE.
      ISIGINT=0
      IF (KPROC .GT. 1) THEN
        CALL MPE_BARRIER(IERR)
        IF(IERR.LT.0) CALL MPEI_ABORT('MPE_BARRIER ERROR in CHESIG')
      ENDIF

      IF (KTST.GE.2) THEN
        WRITE(KU06,*) '   CHESIG:  FOR PE ', KRANK, ' DONE'
        WRITE(KU06,*) '   CHESIG:  LDSIGSTOP:', LDSIGSTOP,
     &                ' LDSIGREST: ', LDSIGREST
        CALL FLUSH(KU06)
      ENDIF

!     ------------------------------------------------------------------

      RETURN
      END SUBROUTINE CHESIG

