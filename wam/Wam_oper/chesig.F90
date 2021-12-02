! ----------------------------------------------------------------------

      SUBROUTINE CHESIG (KU06, KRANK, KPROC, LDSIGSTOP, LDSIGREST)

! ----------------------------------------------------------------------

!**** *CHESIG* - CHECK SIGNALS.

!     PURPOSE.
!     --------
!       CHECK FOR EXTERNAL SIGNALS REQUESTING CREATION OF RESTART FILES.

!**   INTERFACE.
!     ----------
!        *CALL* *CHESIG(KU06,KRANK,KPROC,LDSIGSTOP,LDSIGREST)

!    I/   *KU06*      - LOGICAL I/O UNIT FOR STDOUT.
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
!       MPL_SEND
!       MPL_RECV
!       MPL_BARRIER
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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWSTAT  , ONLY : ISIGHUP  ,ISIGINT
      USE MPL_MODULE

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KU06, KRANK, KPROC
      LOGICAL, INTENT(OUT) :: LDSIGSTOP, LDSIGREST

      INTEGER(KIND=JWIM) :: JE, ITAG, IPROC
      INTEGER(KIND=JWIM) :: IBUF(2)
      INTEGER(KIND=JWIM), PARAMETER :: JPSIGTAG=1000

! ----------------------------------------------------------------------

!*    1. UNBLOCK SIGNALS.
!        ----------------

      CALL IFSSIGB

!         check for arrival of signals

      LDSIGSTOP=.FALSE.
      LDSIGREST=.FALSE.

      IF (KRANK == 1) THEN
        IF (ISIGINT == 1) THEN
          WRITE(*,'(A)')' SIGNAL RECEIVED: WRITE RESTART FILES'
        ENDIF
        IF (ISIGHUP == 1) THEN
          WRITE(*,'(A)')' SIGNAL RECEIVED: STOP'
        ENDIF

        DO JE=2,KPROC
          ITAG=JPSIGTAG+JE
          IBUF(1)=ISIGHUP
          IBUF(2)=ISIGINT
          CALL MPL_SEND(IBUF,KDEST=JE,KTAG=ITAG,CDSTRING='CHESIG:')
        ENDDO
      ELSE
        ITAG=JPSIGTAG+KRANK
        IPROC=1
        CALL MPL_RECV(IBUF,KSOURCE=IPROC,KTAG=ITAG,CDSTRING='CHESIG:')
        ISIGHUP=IBUF(1)
        ISIGINT=IBUF(2)
      ENDIF
      IF (ISIGHUP /= 0) LDSIGSTOP=.TRUE.
      IF (ISIGINT /= 0) LDSIGREST=.TRUE.
      ISIGINT=0
      IF (KPROC > 1) THEN
        CALL MPL_BARRIER(CDSTRING='CHESIG:')
      ENDIF

!     ------------------------------------------------------------------
      END SUBROUTINE CHESIG
