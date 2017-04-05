      SUBROUTINE MPDISTRIBFL(ISEND, ITAG, NBLKS, NBLKE, KINF, KSUP,
     &                       MINF, MSUP, FL)

! ----------------------------------------------------------------------

!****  *MPDISTRIBFL* - DISTRIBUTE FL ACROSS PROCESSORS 

!     J. BIDLOT    ECMWF   SEPTEMBER 1997 

!     PURPOSE.
!     --------
!     SEND THEIR RESPECTIVE CONTRIBUTION OF ARRAY FL 
!     FROM PROCESSOR ISEND TO THE OTHER PE's.

!*    INTERFACE.
!     ----------

!      CALL *MPDISTRIBFL(ISEND, ITAG, NBLKS, NBLKE, KINF, KSUP,
!    &                   MINF, MSUP, FL)

!     *ISEND*     RANK OF THE PROCESS ONTO WHICH FIELD IS COLLECTED 
!     *ITAG*      TAG ASSOCIATED WITH AS A PARTICULAR CALL TO SUBROUTINE
!                 THIS IS NECESSARY TO DIFFERENTIATE THE DIFFERENT CALLS 
!     *NBLKS*     INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NBLKE*     INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *KINF*      INDEX OF THE FIRST DIRECTION OF FL TO BE DISTRIBUTED
!     *KSUP*      INDEX OF THE LAST DIRECTION OF FL TO BE DISTRIBUTED
!     *MINF*      INDEX OF THE FIRST FREQUENCY OF FL TO BE DISTRIBUTED
!     *MSUP*      INDEX OF THE LAST FREQUENCY OF FL TO BE DISTRIBUTED
!     *FL*        INPUT/OUTPUT ARRAY CONTAINING THE PART OF THE SPECTRUM 

!     METHOD.
!     -------
!     MPL SEND OF ARRAY FL FROM PROCESSOR CORRESPONDING TO ISEND TO 
!     ALL OTHER PROCESSORS.

!     EXTERNALS.
!     ----------
!     MPL PACKAGE :
!         MPL_SEND
!         MPL_RECV
!         MPL_ABORT 

!     REFERENCES.
!     -----------
!         NONE
! -------------------------------------------------------------------

      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,MPMAXLENGTH
      USE YOWPARAM , ONLY : NIBLO
      USE YOWUNPOOL, ONLY : LLUNSTR
      USE YOWPD,     ONLY : RANK, MYRANK, NP,  MNP=>NPA, EXCHANGE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE MPL_MODULE

!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: ISEND, ITAG, KINF, KSUP, MINF, MSUP
      INTEGER :: IR, IP, M, K, IJ, MPLENGTH, KCOUNT, KRTAG, KRCOUNT
      INTEGER,DIMENSION(NPROC) :: NBLKS, NBLKE
      INTEGER,DIMENSION(NPROC) :: ISENDREQ

      REAL :: ZHOOK_HANDLE
      REAL,DIMENSION(0:NIBLO,KINF:KSUP,MINF:MSUP) :: FL
      REAL*8, ALLOCATABLE :: AC(:,:,:)
      REAL,ALLOCATABLE :: ZCOMBUFS(:,:),ZCOMBUFR(:)

!----------------------------------------------------------------------

      IF((ISEND.EQ.0).OR.(NPROC.EQ.1)) RETURN

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('MPDISTRIBFL',0,ZHOOK_HANDLE)
#endif

      MPLENGTH=MPMAXLENGTH*(KSUP-KINF+1)*(MSUP-MINF+1)
      ALLOCATE(ZCOMBUFR(MPLENGTH))

      IF(IRANK.EQ.ISEND) THEN
!     1.1 SEND NON BLOCKING TO ALL PROCESSORS
!         -----------------------------------
        ALLOCATE(ZCOMBUFS(MPLENGTH,NPROC))

        DO IP=1,NPROC
          KCOUNT=0
          DO M=MINF,MSUP
            DO K=KINF,KSUP
              DO IJ=NBLKS(IP),NBLKE(IP)
                KCOUNT=KCOUNT+1
                ZCOMBUFS(KCOUNT,IP)=FL(IJ,K,M)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

        DO IP=1,NPROC
          CALL MPL_SEND(ZCOMBUFS(1:MPLENGTH,IP),KDEST=IP,KTAG=ITAG,
     &     KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(IP),
     &     CDSTRING='MPDISTRIBFL SEND:')
        ENDDO
      ENDIF

!     1.2 RECEIVE CONTRIBUTION TO THE FIELD FROM PE ISEND 
!         ------------------------------------------------ 

      CALL MPL_RECV(ZCOMBUFR(1:MPLENGTH),KSOURCE=ISEND,KTAG=ITAG,
     &   KOUNT=KRCOUNT,KRECVTAG=KRTAG,CDSTRING='MPDISTRIBFL:')
      IF(KRCOUNT.NE.MPLENGTH) CALL MPL_ABORT
     &   ('MPL_RECV ERROR in MPDISTRIBFL:MISMATCHED MESSAGE LENGTH')
      IF(KRTAG.NE.ITAG) CALL MPL_ABORT
     &   ('MPL_RECV ERROR in MPDISTRIBFL MISMATCHED TAGS' )

      IF (LLUNSTR) THEN
        ALLOCATE(AC(MNP,KINF:KSUP,MINF:MSUP))
        DO M=MINF,MSUP
          DO K=KINF,KSUP
            DO IP = 1, RANK(IRANK)%NP
              AC(IP,K,M) = DBLE(ZCOMBUFR(IP))
            ENDDO
          ENDDO 
        ENDDO
      ELSE
        KCOUNT=0
        DO M=MINF,MSUP
          DO K=KINF,KSUP
            DO IJ=NBLKS(IRANK),NBLKE(IRANK)
              KCOUNT=KCOUNT+1
              FL(IJ,K,M)=ZCOMBUFR(KCOUNT)
            ENDDO
          ENDDO 
        ENDDO 
      ENDIF

!     1.3 WAIT ANY OUTSTANDING SENDS TO COMPLETE
!         --------------------------------------

      IF(IRANK.EQ.ISEND) THEN
        CALL MPL_WAIT(KREQUEST=ISENDREQ,
     &                CDSTRING='MPDISTRIBFL:')
        DEALLOCATE(ZCOMBUFS)
      ENDIF

      DEALLOCATE(ZCOMBUFR)

!!! it's not very prety but when unstructured, you also need to have the values on the halo
!!! it should be combined with the previous exchange
      IF (LLUNSTR) THEN
        DO M=MINF,MSUP
          DO K=KINF,KSUP
            CALL EXCHANGE(AC(:,K,M))
            DO IP = 1 , MNP
              FL(IP,K,M) = AC(IP,K,M)
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(AC)
      ENDIF


#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('MPDISTRIBFL',1,ZHOOK_HANDLE)
#endif

      RETURN
      END SUBROUTINE MPDISTRIBFL
