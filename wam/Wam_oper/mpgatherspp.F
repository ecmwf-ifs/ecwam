      SUBROUTINE MPGATHERSPP(IRECV,ITAG,NSTART,NEND,NSPFLD,NSCFLD,
     1                      FL3,FL1,FLPTS,
     2                      EMEAN,FMEAN,THQ,SMEAN,ESWELL,FSWELL,THSWELL, 
     3                      THWISEA,U10,THW,US) 
C
C****  *MPGATHERSSPP* - GATHERS	SPECTRUM AND SCALAR FIELDS ONTO A SINGLE
C****                   PROCESS FOR OUTPUT IN OUTSPP AT SELECTED GRID
C****                   POINTS
C
!     J. BIDLOT    ECMWF   APRIL 1996  MESSAGE PASSING
C
C     PURPOSE.
C     --------
C
C     GATHER NSPFLD SPECTRUM AND NSCFLD SCALAR FIELDS DISTRIBUTED ACROSS THE
C     DIFFERENT PROCESSES ONTO THE SINGLE PROCESS IRECV FOR SELECTED GRID
C     POINTS.
C
C*    INTERFACE.
C     ----------
C
C     CALL *MPGATHERSPP*(IRECV,ITAG,NSTART,NEND,NSPFLD,NSCFLD, 
C    1                   FL1,FL3,FLPTS,
C    2                   EMEAN,FMEAN,THQ,SMEAN,ESWELL,FSWELL,THSWELL, 
C    3                   THWISEA,U10,THW,US)
C
C
C     *IRECV*     RANK OF THE PROCESS ONTO WHICH FIELD IS COLLECTED 
C     *ITAG*      TAG ASSOCIATED WITH A PARTICULAR CALL TO SUBROUTINE
C                 THIS IS NECESSARY TO DIFFERENTIATE THE DIFFERENT CALLS 
C                 THE USE OF A SEPARATE TAG FOR EACH SPECIFIC CALL SHOULD
C                 ENSURE THAT NO CONFLICT ARISES BETWEEN DIFFERENT MESSAGES
C     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C     *NSPFLD*    NUMBER OF SPECTRUM FIELDS
C     *NSCFLD*    NUMBER OF SCALAR FILEDS
C     *FL3*       BLOCK OF SPECTRA.
C     *FL1*       SWELL SPECTRUM.
C     *FLPTS*     SPECTRA AT THE SELECTED GRID POINTS (OUTPUT) 
C     *EMEAN*     TOTAL ENERGY.
C     *FMEAN*     MEAN FREQUENCY (1./MEAN PERIOD).
C     *THQ*       MEAN WAVE DIRECTION IN RADIANS.
C     *SMEAN*     MEAN SQUARE SLOPE.
C     *ESWELL*    TOTAL SWELL ENERGY
C     *FSWELL*    MEAN SWELL FREQUENCY
C     *THSWELL*   MEAN SWELL DIRECTION
C     *THWISEA*   MEAN WIND INTERACTING SEA DIRECTION
C     *U10*    NEW WIND SPEED IN M/S.
C     *THW*    WIND DIRECTION IN RADIANS
C     *US*     NEW FRICTION VELOCITY IN M/S.
C
C     METHOD.
C     -------
C     MPE SEND OF A MESSAGE BUFFER CONTAINING ALL THE FIELDS AT ALL SELECTED
C     GRID POINTS HELD BY THAT PROCESS TO THE PROCESS CORRESPONDING TO IRECV
C     FOR ALL PROCESS EXCEPT FOR THE PROCESS CORRESPONDING TO IRECV
C     WHERE THE MESSAGE BUFFERS ARE RECEIVED AND REORDERED.
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
#include "comcout.h"
C
#include "comtest.h"
C
#include "commpp.h"
C
#include "txtmpp.h"

!----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL1, FL3
      REAL,DIMENSION(NSPFLD,NGOUT,NANG,NFRE) :: FLPTS
      REAL,DIMENSION(NIBLO) :: EMEAN,FMEAN,THQ,SMEAN,ESWELL,FSWELL,
     1                         THSWELL,THWISEA
      REAL,DIMENSION(NIBLO) :: U10,THW,US  
!----------------------------------------------------------------------

      REAL,ALLOCATABLE :: ZCOMBUF(:)

c
      MAXLENGTH=(NSCFLD+NSPFLD*NANG*NFRE)*NGOUT
      ALLOCATE(ZCOMBUF(MAXLENGTH))

C     1.0 DEFAULT ACTION IF NO FIELD GATHERING 
C         ------------------------------------
      IF((IRECV.EQ.0).OR.(NPROC.EQ.1)) THEN
         DO NGOU=1,NGOUT
           IJ=IJAR(NGOU)
           DO  M = 1, NFRE
             DO K = 1, NANG
                FLPTS(1,NGOU,K,M)=FL3(IJ,K,M)
             END DO
           END DO
           DO  M = 1, NFRE
             DO K = 1, NANG
               FLPTS(2,NGOU,K,M)=FL1(IJ,K,M)
             END DO
           END DO
         END DO

      ELSE IF(IRANK.NE.IRECV) THEN
C     1.1 SEND TO THE PROCESS THAT GATHERS THE WHOLE FIELD
C         ------------------------------------------------
C
C     LOAD COMMUNICATION BUFFER
         KCOUNT=0
         DO NGOU=1,NGOUT

           IJ=IJAR(NGOU)
           IF((IJ.GE.NSTART(IRANK)).AND.(IJ.LE.NEND(IRANK))) THEN
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=EMEAN(IJ)
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=FMEAN(IJ)
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=SMEAN(IJ)
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=THQ(IJ)
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=ESWELL(IJ)
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=FSWELL(IJ)
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=THSWELL(IJ)
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=THWISEA(IJ)
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=U10(IJ)
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=THW(IJ)
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=US(IJ)

              DO  M = 1, NFRE
                DO K = 1, NANG
                  KCOUNT=KCOUNT+1
                  ZCOMBUF(KCOUNT)=FL3(IJ,K,M)
                END DO
              END DO

              DO  M = 1, NFRE
                DO K = 1, NANG
                  KCOUNT=KCOUNT+1
                  ZCOMBUF(KCOUNT)=FL1(IJ,K,M)
                END DO
              END DO
           ENDIF

         END DO

C     SEND BUFFER
         CALL MPE_SEND(ZCOMBUF,KCOUNT,2,IRECV,ITAG,0,0,0,IERR)
         IF(IERR.LT.0) CALL MPEI_ABORT
     1                    ('MPE_SEND ERROR in MPGATHERSPP' )


       ELSE
C     1.2.1  RECEIVE CONTRIBUTION TO THE FIELD FROM THE OTHER PROCESSES
C            ---------------------------------------------------------- 
         DO IP=1,NPROC-1
            CALL MPE_RECV(ZCOMBUF,MAXLENGTH,2,-1,ITAG,0,0,0,
     $                    KRCOUNT,KRFROM,KRTAG,IERR)
            IF(IERR.LT.0) CALL MPEI_ABORT
     1                       ('MPE_RECV ERROR in MPGATHERSPP' )
            IF(KRTAG.NE.ITAG) CALL MPEI_ABORT
     1        ('MPE_RECV ERROR in MPGATHERSPP:  MISMATCHED TAGS' )

C  DECODE BUFFER
            KCOUNT=0

            DO NGOU=1,NGOUT
              IJ=IJAR(NGOU)
              IF((IJ.GE.NSTART(KRFROM)).AND.(IJ.LE.NEND(KRFROM))) THEN
                KCOUNT=KCOUNT+1
                EMEAN(IJ)=ZCOMBUF(KCOUNT)
                KCOUNT=KCOUNT+1
                FMEAN(IJ)=ZCOMBUF(KCOUNT)
                KCOUNT=KCOUNT+1
                SMEAN(IJ)=ZCOMBUF(KCOUNT)
                KCOUNT=KCOUNT+1
                THQ(IJ)=ZCOMBUF(KCOUNT)
                KCOUNT=KCOUNT+1
                ESWELL(IJ)=ZCOMBUF(KCOUNT)
                KCOUNT=KCOUNT+1
                FSWELL(IJ)=ZCOMBUF(KCOUNT)
                KCOUNT=KCOUNT+1
                THSWELL(IJ)=ZCOMBUF(KCOUNT)
                KCOUNT=KCOUNT+1
                THWISEA(IJ)=ZCOMBUF(KCOUNT)
                KCOUNT=KCOUNT+1
                U10(IJ)=ZCOMBUF(KCOUNT)
                KCOUNT=KCOUNT+1
                THW(IJ)=ZCOMBUF(KCOUNT)
                KCOUNT=KCOUNT+1
                US(IJ)=ZCOMBUF(KCOUNT)

                DO  M = 1, NFRE
                  DO K = 1, NANG
                    KCOUNT=KCOUNT+1
                    FLPTS(1,NGOU,K,M)=ZCOMBUF(KCOUNT)
                  END DO
                END DO

                DO  M = 1, NFRE
                  DO K = 1, NANG
                    KCOUNT=KCOUNT+1
                    FLPTS(2,NGOU,K,M)=ZCOMBUF(KCOUNT)
                  END DO
                END DO

              ENDIF
            END DO
                IF(KRCOUNT.NE.KCOUNT) CALL MPEI_ABORT
     1        ('MPE_RECV ERROR in MPGATHERSPP: WRONG LENGTH OF MESSAGE')

         END DO
C     1.2.2  CONTRIBUTION FROM RECEIVING PROCESS 
C            ----------------------------------- 

         DO NGOU=1,NGOUT
           IJ=IJAR(NGOU)
           IF((IJ.GE.NSTART(IRANK)).AND.(IJ.LE.NEND(IRANK)))THEN
             DO  M = 1, NFRE
               DO K = 1, NANG
                 FLPTS(1,NGOU,K,M)=FL3(IJ,K,M)
               END DO
             END DO
             DO  M = 1, NFRE
               DO K = 1, NANG
                 FLPTS(2,NGOU,K,M)=FL1(IJ,K,M)
               END DO
             END DO
           ENDIF
         END DO

      ENDIF

      DEALLOCATE(ZCOMBUF)
c
      RETURN
      END

