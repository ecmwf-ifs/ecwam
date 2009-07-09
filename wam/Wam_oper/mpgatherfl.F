      SUBROUTINE MPGATHERFL(IRECV,NSTART,NEND,KINF,KSUP,
     &                      MINF,MSUP,FL)

!****  *MPGATHERFL* - GATHER FL ONTO A SINGLE PROCESS 

!     J. BIDLOT    ECMWF   APRIL 1996  MESSAGE PASSING
!     J. BIDLOT    ECMWF   MARCH 1997  add use of MINF AND MSUP
!     J. BIDLOT    ECMWF   JANUARY 2003 use of MPL_GATHERV
!     J. BIDLOT    ECMWF   JANUARY 2009 use of KINF AND KSUP

!     PURPOSE.
!     --------

!     GATHER ARRAY FL DISTRIBUTED ACROSS THE
!     DIFFERENT PROCESSES ONTO THE SINGLE PROCESS IRECV

!*    INTERFACE.
!     ----------

!     CALL *MPGATHERFL*(IRECV,NSTART,NEND,KINF,KSUP,
!                       MINF,MSUP,FL) 

!     *IRECV*     RANK OF THE PROCESS ONTO WHICH FIELD IS COLLECTED 
!     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *KINF*      INDEX OF THE FIRST DIRECTION OF FL
!     *KSUP*      INDEX OF THE LAST DIRECTION OF FL
!     *MINF*      INDEX OF THE FIRST FREQUENCY OF FL
!     *MSUP*      INDEX OF THE LAST FREQUENCY OF FL
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

!     REFERENCES.
!     -----------
!         NONE
! -------------------------------------------------------------------

      USE YOWPARAM , ONLY : NIBLO
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE MPL_MODULE

!----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(0:NIBLO,KINF:KSUP,MINF:MSUP) :: FL

!----------------------------------------------------------------------

      REAL,ALLOCATABLE,DIMENSION(:) :: PSENDBUF,PRECVBUF

      INTEGER,DIMENSION(NPROC) :: KRECVCOUNTS

      N23=(KSUP-KINF+1)*(MSUP-MINF+1)

      IF((IRECV.EQ.0).OR.(NPROC.EQ.1)) THEN
        RETURN
      ELSE
        ILEN=(NEND(IRANK)-NSTART(IRANK)+1)*N23
        ALLOCATE(PSENDBUF(ILEN))
        KCOUNT=0
        DO M=MINF,MSUP
          DO K=KINF,KSUP
            DO IJ=NSTART(IRANK),NEND(IRANK)
              KCOUNT=KCOUNT+1
              PSENDBUF(KCOUNT)=FL(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO

        ILENR=0
        DO IP=1,NPROC
          KRECVCOUNTS(IP)=(NEND(IP)-NSTART(IP)+1)*N23
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
              DO M=MINF,MSUP
                DO K=KINF,KSUP
                  DO IJ=NSTART(IP),NEND(IP)
                    KCOUNT=KCOUNT+1
                    FL(IJ,K,M)=PRECVBUF(KCOUNT)
                  ENDDO
                ENDDO 
              ENDDO 
          ENDDO 
        ENDIF

        DEALLOCATE(PSENDBUF)
        DEALLOCATE(PRECVBUF)

      ENDIF

      RETURN
      END SUBROUTINE MPGATHERFL
