      SUBROUTINE MPDISTRIBINTFLD(ISEND,ITAG,NSTART,NEND,IFIELD)
C
C ----------------------------------------------------------------------
C
C****  *MPDISTRIBINTFLD* - DISTRIBUTE IFIELD ACROSS PROCESSORS 
C
C     J. BIDLOT    ECMWF   SEPTEMBER 1997 
C
C     PURPOSE.
C     --------
C     SEND THE RESPECTIVE CONTRIBUTION OF ARRAY IFIELD FROM PROCESSOR ISEND
C     TO THE OTHER PE's.
C
C*    INTERFACE.
C     ----------
C
C     CALL *MPDISTRIBFL*(ISEND,ITAG,NSTART,NEND,MI,MS,IFIELD) 
C
C     *ISEND*     RANK OF THE PROCESS ONTO WHICH IFIELD IS COLLECTED 
C     *ITAG*      TAG ASSOCIATED WITH AS A PARTICULAR CALL TO SUBROUTINE
C                 THIS IS NECESSARY TO DIFFERENTIATE THE DIFFERENT CALLS 
C     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C     *IFIELD*     INPUT/OUTPUT ARRAY
C
C     METHOD.
C     -------
C     MPE SEND OF ARRAY IFIELD FROM PROCESSOR CORRESPONDING TO ISEND TO 
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
      INTEGER,DIMENSION(NIBLO) :: IFIELD
C----------------------------------------------------------------------

      INTEGER,ALLOCATABLE :: ICOMBUF(:)

c
      IF((ISEND.EQ.0).OR.(NPROC.EQ.1)) THEN
         RETURN
      ELSE IF(IRANK.EQ.ISEND) THEN
C     1.1 SEND TO ALL PROCESSORS OTHER THAN ISEND
C         ------------------------------------------------
         MPLENGTH=MPMAXLENGTH
         ALLOCATE(ICOMBUF(MPLENGTH))

         DO IP=1,ISEND-1
            KCOUNT=0
            DO IJ=NSTART(IP),NEND(IP)
               KCOUNT=KCOUNT+1
               ICOMBUF(KCOUNT)=IFIELD(IJ)
            END DO
            CALL MPE_SEND(ICOMBUF,MPLENGTH,1,IP,ITAG,0,0,0,IERR)
            IF(IERR.LT.0) CALL MPEI_ABORT
     1                    ('MPE_SEND ERROR 1 in MPDISTRIBINTFLD' )
         ENDDO

         DO IP=ISEND+1,NPROC
            KCOUNT=0
            DO IJ=NSTART(IP),NEND(IP)
               KCOUNT=KCOUNT+1
               ICOMBUF(KCOUNT)=IFIELD(IJ)
            END DO
            CALL MPE_SEND(ICOMBUF,MPLENGTH,1,IP,ITAG,0,0,0,IERR)
            IF(IERR.LT.0) CALL MPEI_ABORT
     1                    ('MPE_SEND ERROR 2 in MPDISTRIBINTFLD' )
         ENDDO

       ELSE

C     1.2 RECEIVE CONTRIBUTION TO THE IFIELD FROM PE ISEND 
C         ------------------------------------------------ 
C
         MPLENGTH=MPMAXLENGTH
         ALLOCATE(ICOMBUF(MPLENGTH))

            CALL MPE_RECV(ICOMBUF,MPLENGTH,1,ISEND,ITAG,0,0,0,
     $                    KRCOUNT,KRFROM,KRTAG,IERR)
            IF(IERR.LT.0) CALL MPEI_ABORT
     1                       ('MPE_RECV ERROR in MPDISTRIBINTFLD' )
            IF(KRCOUNT.NE.MPLENGTH) CALL MPEI_ABORT
     1      ('MPE_RECV ERROR in MPDISTRIBINTFLD:MISMATCHED MSG LENGTH')
            IF(KRTAG.NE.ITAG) CALL MPEI_ABORT
     1      ('MPE_RECV ERROR in MPDISTRIBINTFLD MISMATCHED TAGS' )

            KCOUNT=0
            DO IJ=NSTART(IRANK),NEND(IRANK)
               KCOUNT=KCOUNT+1
               IFIELD(IJ)=ICOMBUF(KCOUNT)
            END DO

      ENDIF

      DEALLOCATE(ICOMBUF)
c
      RETURN
      END
