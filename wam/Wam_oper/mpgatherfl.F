      SUBROUTINE MPGATHERFL(IRECV,ITAG,NSTART,NEND,MINF,MSUP,MI,MS,FL) 

!****  *MPGATHERFL* - GATHER FL ONTO A SINGLE PROCESS 

!     J. BIDLOT    ECMWF   APRIL 1996  MESSAGE PASSING
!     J. BIDLOT    ECMWF   MARCH 1997  add use of MINF AND MSUP
!     J. BIDLOT    ECMWF   JANUARY 2003 use of MPL_GATHERV

!     PURPOSE.
!     --------

!     GATHER ARRAY FL DISTRIBUTED ACROSS THE
!     DIFFERENT PROCESSES ONTO THE SINGLE PROCESS IRECV

!*    INTERFACE.
!     ----------

!     CALL *MPGATHERFL*(IRECV,ITAG,NSTART,NEND,MINF,MSUP,MI,MS,FL) 

!     *IRECV*     RANK OF THE PROCESS ONTO WHICH FIELD IS COLLECTED 
!     *ITAG*      TAG ASSOCIATED WITH AS A PARTICULAR CALL TO SUBROUTINE
!                 THIS IS NECESSARY TO DIFFERENTIATE THE DIFFERENT CALLS 
!     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *MINF*      INDEX OF THE FIRST USEFUL FREQUENCY OF FL
!     *MSUP*      INDEX OF THE LAST USEFUL FREQUENCY OF FL
!     *MI*        INDEX OF THE FIRST FREQUENCY OF FL
!     *MS*        INDEX OF THE LAST FREQUENCY OF FL
!     *FL*        INPUT/OUTPUT ARRAY CONTAINING THE PART OF THE SPECTRUM 

!     METHOD.
!     -------
!     MPL SEND OF ARRAY FL TO PROCESS CORRESPONDING TO IRECV FOR
!     ALL PROCESS EXCEPT FOR THE PROCESS CORRESPONDING TO IRECV
!     WHERE IT IS RECEIVED.

!     EXTERNALS.
!     ----------
!     MPL PACKAGE :
!         MPL_GATHERV
!         or
!         MPL_SEND
!         MPL_RECV
!         MPL_ABORT 

!     REFERENCES.
!     -----------
!         NONE
! -------------------------------------------------------------------

      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLO
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,MPMAXLENGTH
      USE MPL_MODULE

!----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(0:NIBLO,NANG,MI:MS) :: FL

!----------------------------------------------------------------------

      REAL,ALLOCATABLE,DIMENSION(:) :: ZCOMBUF,PSENDBUF,PRECVBUF

      INTEGER,DIMENSION(NPROC) :: KRECVCOUNTS

      IRTAG=ITAG+NFRE

      IF((IRECV.EQ.0).OR.(NPROC.EQ.1)) THEN
        RETURN

      ELSEIF(MINF.EQ.MSUP .AND. MI.EQ.MS) THEN
!     if the FL array is essentially 2-dimensional then use 
!     a call to MPL_GATHERV which should be more efficient then 
!     a send and a loop on receive

        ILEN=(NEND(IRANK)-NSTART(IRANK)+1)*NANG
        ALLOCATE(PSENDBUF(ILEN))
        IF(IRANK.NE.IRECV) THEN
          KCOUNT=0
          DO K=1,NANG
            DO IJ=NSTART(IRANK),NEND(IRANK)
              KCOUNT=KCOUNT+1
              PSENDBUF(KCOUNT)=FL(IJ,K,MINF)
            ENDDO
          ENDDO
        ENDIF

        ILENR=0
        DO IP=1,NPROC
         KRECVCOUNTS(IP)=(NEND(IP)-NSTART(IP)+1)*NANG
         ILENR=ILENR+KRECVCOUNTS(IP)
        ENDDO
        ALLOCATE(PRECVBUF(ILENR))

        CALL MPL_GATHERV(PSENDBUF=PSENDBUF,
     &                   KRECVCOUNTS=KRECVCOUNTS,
     &                   PRECVBUF=PRECVBUF,
     &                   KROOT=IRECV,
     &                   CDSTRING='MPGATHERFL:')

        IF(IRANK.EQ.IRECV) THEN
          KCOUNT=0
          DO IP=1,NPROC
            IF(IP.NE.IRECV) THEN
              DO K=1,NANG
                DO IJ=NSTART(IP),NEND(IP)
                  KCOUNT=KCOUNT+1
                  FL(IJ,K,MINF)=PRECVBUF(KCOUNT)
                ENDDO
              ENDDO 
            ELSE
              KCOUNT=KCOUNT+(NEND(IP)-NSTART(IP)+1)*NANG
            ENDIF
          ENDDO 
        ENDIF

        DEALLOCATE(PSENDBUF)
        DEALLOCATE(PRECVBUF)

!! the old way is below
      ELSEIF(IRANK.NE.IRECV) THEN
!     1.1 SEND TO THE PROCESS THAT GATHERS THE WHOLE FIELD
!         ------------------------------------------------
        MPLENGTH=MPMAXLENGTH*NANG*(MSUP-MINF+1)
        ALLOCATE(ZCOMBUF(MPLENGTH))

        KCOUNT=0
        DO M=MINF,MSUP
          DO K=1,NANG
            DO IJ=NSTART(IRANK),NEND(IRANK)
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=FL(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO

!       synchronisation

        CALL GSTATS(677,0)
        CALL MPL_BARRIER(CDSTRING='MPGATHERFL:')

! SEND CONTRIBUTION TO RECEIVING PE

        CALL MPL_SEND(ZCOMBUF(1:MPLENGTH),KDEST=IRECV,KTAG=ITAG,
     &   CDSTRING='MPGATHERFL:')
        CALL GSTATS(677,1)

        DEALLOCATE(ZCOMBUF)
      ELSE
!     1.2 RECEIVE CONTRIBUTION TO THE FIELD FROM THE OTHER PROCESSES
!         ---------------------------------------------------------- 

!       synchronisation

        CALL GSTATS(677,0)
        CALL MPL_BARRIER(CDSTRING='MPGATHERFL:')
        CALL GSTATS(677,1)

! GET CONTRIBUTION FROM OTHER PE'S
        MPLENGTH=MPMAXLENGTH*NANG*(MSUP-MINF+1)
        ALLOCATE(ZCOMBUF(MPLENGTH))

        DO IP=1,NPROC-1
          CALL GSTATS(677,0)
          CALL MPL_RECV(ZCOMBUF(1:MPLENGTH),KTAG=ITAG,KFROM=KRFROM,
     &     KOUNT=KRCOUNT,KRECVTAG=KRTAG,CDSTRING='MPGATHERFL:')
          IF(KRCOUNT.NE.MPLENGTH) CALL MPL_ABORT
     &     ('MPL_RECV ERROR in MPGATHERFL:MISMATCHED MESSAGE LENGTH')
          IF(KRTAG.NE.ITAG) CALL MPL_ABORT
     &     ('MPL_RECV ERROR in MPGATHERFL MISMATCHED TAGS' )
          CALL GSTATS(677,1)

          KCOUNT=0
          DO M=MINF,MSUP
            DO K=1,NANG
              DO IJ=NSTART(KRFROM),NEND(KRFROM)
                KCOUNT=KCOUNT+1
                FL(IJ,K,M)=ZCOMBUF(KCOUNT)
              ENDDO
            ENDDO 
          ENDDO 

        ENDDO

        DEALLOCATE(ZCOMBUF)
      ENDIF


      RETURN
      END SUBROUTINE MPGATHERFL
