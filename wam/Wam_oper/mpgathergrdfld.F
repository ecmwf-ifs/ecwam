      SUBROUTINE MPGATHERGRDFLD(ISEND,IRECV,ITAG,FIELD) 
C
C****  *MPGATHERGRDFLD* - GATHER GRID DATA FIELD FROM PROCESS ISEND
C****                     ONTO PROCESS IRECV 
C
!     J. BIDLOT    ECMWF   APRIL 1996  MESSAGE PASSING
C
C     PURPOSE.
C     --------
C
C     TRANSMIT GRID (2-D) DATA CONTAINED IN THE ARRAY FIELD THAT MIGHT BE 
C     ON A DIFFERENT PROCESS TO THE RECEIVING PROCESS
C
C*    INTERFACE.
C     ----------
C
C     CALL *MPGATHERGRDFLD*(ISEND,IRECV,ITAG,FIELD) 
C
C     *ISEND*     PROCESS RANK SENDING THE GRID FIELD 
C     *IRECV*     PROCESS RANK RECEIVING THE GRID FIELD 
C     *ITAG*      TAG USED TO DIFFERENTIATE THE CALLS TO THIS SUBROUTINE 
C     *FIELD*     INPUT/OUTPUT ARRAY CONTAINING THE FIELD
C
C     METHOD.
C     -------
C     MPE SEND to PROCESS IRECV FROM PROCESS ISEND 
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

!----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      REAL,DIMENSION(NGX,NGY) :: FIELD
!----------------------------------------------------------------------

      REAL,DIMENSION(NGX*NGY):: ZCOMBUF

c
      KCOUNT=NGX*NGY 

      IF((ISEND.EQ.IRECV).OR.(NPROC.EQ.1)) THEN
         RETURN

      ELSE IF(IRANK.EQ.ISEND) THEN 
C     1.1 SEND FIELD FROM PROCESS ISEND TO PROCESS IRECV
C         -----------------------------------------------
C
C        LOAD COMMUNICATION BUFFER
          IJ=0
          DO J=1,NGY
          DO I=1,NGX
             IJ=IJ+1
             ZCOMBUF(IJ)=FIELD(I,J)
          END DO
          END DO

          CALL MPE_SEND(ZCOMBUF,KCOUNT,2,IRECV,ITAG,0,0,0,IERR)
          IF(IERR.LT.0) CALL MPEI_ABORT
     1                  ('MPE_SEND ERROR in MPGATHERGRDFLD' )

       ELSE IF(IRANK.EQ.IRECV) THEN
C     1.2 RECEIVE FIELD FROM PROCESS ISEND ON PROCESS IRECV
C         ------------------------------------------------- 
          CALL MPE_RECV(ZCOMBUF,KCOUNT,2,ISEND,ITAG,0,0,0,
     $                  KRCOUNT,KRFROM,KRTAG,IERR)
          IF(IERR.LT.0) CALL MPEI_ABORT
     1                 ('MPE_RECV ERROR in MPGATHERGRDFLD' )
          IF(KRCOUNT.NE.KCOUNT) CALL MPEI_ABORT
     1    ('MPE_RECV ERROR in MPGATHERGRDFLD:MISMATCHED MESSAGE LENGTH')
          IF(KRTAG.NE.ITAG) CALL MPEI_ABORT
     1    ('MPE_RECV ERROR in MPGATHERGRDFLD: MISMATCHED TAGS' )
C
C      DECODE COMMUNICATION BUFFER
          IJ=0
          DO J=1,NGY
          DO I=1,NGX
             IJ=IJ+1
             FIELD(I,J)=ZCOMBUF(IJ)
          ENDDO
          ENDDO
      ENDIF

      RETURN
      END

