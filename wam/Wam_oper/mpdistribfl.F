      SUBROUTINE MPDISTRIBFL(ISEND,ITAG,NSTART,NEND,MINF,MSUP,MI,MS,FL)
C
C ----------------------------------------------------------------------
C
C****  *MPDISTRIBFL* - DISTRIBUTE FL ACROSS PROCESSORS 
C
C     J. BIDLOT    ECMWF   SEPTEMBER 1997 
C
C     PURPOSE.
C     --------
C     SEND THEIR RESPECTIVE CONTRIBUTION OF ARRAY FL FROM PROCESSOR ISEND
C     TO THE OTHER PE's.
C
C*    INTERFACE.
C     ----------
C
C     CALL *MPDISTRIBFL*(ISEND,ITAG,NSTART,NEND,MINF,MSUP,MI,MS,FL)
C
C     *ISEND*     RANK OF THE PROCESS ONTO WHICH FIELD IS COLLECTED 
C     *ITAG*      TAG ASSOCIATED WITH AS A PARTICULAR CALL TO SUBROUTINE
C                 THIS IS NECESSARY TO DIFFERENTIATE THE DIFFERENT CALLS 
C     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C     *MINF*      INDEX OF THE FIRST FREQUENCY OF FL WHICH WILL BE DISTRIBUTED
C     *MSUP*      INDEX OF THE LAST FREQUENCY OF FL WHICH WILL BE DISTRIBUTED
C     *MI*        INDEX OF THE FIRST FREQUENCY OF FL
C     *MS*        INDEX OF THE LAST FREQUENCY OF FL
C     *FL*        INPUT/OUTPUT ARRAY CONTAINING THE PART OF THE SPECTRUM 
C
C     METHOD.
C     -------
C     MPE SEND OF ARRAY FL FROM PROCESSOR CORRESPONDING TO ISEND TO 
C     ALL OTHER PROCESSORS.
C     
C     EXTERNALS.
C     ----------
C     MPE PACKAGE :
C         MPE_SEND
C         MPE_RECV
C         MPEI_ABORT 
C
C     REFERENCES.
C     -----------
C         NONE
C -------------------------------------------------------------------
C
#include "param.h"
C
#include "comtest.h"
C
#include "commpp.h"
C
#include "txtmpp.h"

C----------------------------------------------------------------------
C     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(0:NIBLO,NANG,MI:MS) :: FL
C----------------------------------------------------------------------

      REAL,ALLOCATABLE :: ZCOMBUF(:)

C
      IF((ISEND.EQ.0).OR.(NPROC.EQ.1)) THEN
         RETURN
      ELSE IF(IRANK.EQ.ISEND) THEN
C     1.1 SEND TO ALL PROCESSORS OTHER THAN ISEND
C         ------------------------------------------------
         MPLENGTH=MPMAXLENGTH*NANG*(MSUP-MINF+1)
         ALLOCATE(ZCOMBUF(MPLENGTH))

         DO IP=1,ISEND-1
            KCOUNT=0
            DO M=MINF,MSUP
              DO K=1,NANG
                 DO IJ=NSTART(IP),NEND(IP)
                    KCOUNT=KCOUNT+1
                    ZCOMBUF(KCOUNT)=FL(IJ,K,M)
                 END DO
              END DO
            END DO
            CALL MPE_SEND(ZCOMBUF,MPLENGTH,2,IP,ITAG,0,0,0,IERR)
            IF(IERR.LT.0) CALL MPEI_ABORT
     1                    ('MPE_SEND ERROR 1 in MPDISTRIBFL' )
         ENDDO

         DO IP=ISEND+1,NPROC
            KCOUNT=0
            DO M=MINF,MSUP
              DO K=1,NANG
                 DO IJ=NSTART(IP),NEND(IP)
                    KCOUNT=KCOUNT+1
                    ZCOMBUF(KCOUNT)=FL(IJ,K,M)
                 END DO
              END DO
            END DO
            CALL MPE_SEND(ZCOMBUF,MPLENGTH,2,IP,ITAG,0,0,0,IERR)
            IF(IERR.LT.0) CALL MPEI_ABORT
     1                    ('MPE_SEND ERROR 2 in MPDISTRIBFL' )
         ENDDO

       ELSE

C     1.2 RECEIVE CONTRIBUTION TO THE FIELD FROM PE ISEND 
C         ------------------------------------------------ 
C
         MPLENGTH=MPMAXLENGTH*NANG*(MSUP-MINF+1)
         ALLOCATE(ZCOMBUF(MPLENGTH))

            CALL MPE_RECV(ZCOMBUF,MPLENGTH,2,ISEND,ITAG,0,0,0,
     $                    KRCOUNT,KRFROM,KRTAG,IERR)
            IF(IERR.LT.0) CALL MPEI_ABORT
     1                       ('MPE_RECV ERROR in MPDISTRIBFL' )
            IF(KRCOUNT.NE.MPLENGTH) CALL MPEI_ABORT
     1      ('MPE_RECV ERROR in MPDISTRIBFL:MISMATCHED MESSAGE LENGTH')
            IF(KRTAG.NE.ITAG) CALL MPEI_ABORT
     1      ('MPE_RECV ERROR in MPDISTRIBFL MISMATCHED TAGS' )

            KCOUNT=0
            DO M=MINF,MSUP
              DO K=1,NANG
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  FL(IJ,K,M)=ZCOMBUF(KCOUNT)
                END DO
              END DO 
            END DO 

      ENDIF

      DEALLOCATE(ZCOMBUF)
c
      RETURN
      END
