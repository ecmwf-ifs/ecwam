      SUBROUTINE MPBCASTSCFLD(ITAG,NSTART,NEND,FIELD) 
C
C****  *MPBCASTSCFLD* - BROADCAST SCALAR BLOCK DATA FIELD ONTO ALL
C****                   PROCESS 
C
!     J. BIDLOT    ECMWF   SEPTEMBER 1996  MESSAGE PASSING
C
C     PURPOSE.
C     --------
C
C     BROADCAST SCALAR BLOCK DATA CONTAINED IN ARRAY FIELD DISTRIBUTED 
C     ACROSS THE DIFFERENT PROCESSES ONTO ALL PE's
C
C*    INTERFACE.
C     ----------
C
C     CALL *MPGATHERSCFLD*(ITAG,NSTART,NEND,FIELD) 
C
C     *ITAG*      TAG ASSOCIATED WITH A PARTICULAR CALL TO SUBROUTINE
C                 THIS IS NECESSARY TO DIFFERENTIATE THE DIFFERENT CALLS 
C                 THE USE OF A SEPARATE TAG FOR EACH SPECIFIC CALL SHOULD
C                 ENSURE THAT NO CONFLICT ARISES BETWEEN DIFFERENT MESSAGES
C     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C     *FIELD*     INPUT/OUTPUT ARRAY CONTAINING THE FIELD
C
C     METHOD.
C     -------
C     MPE BROADCAST OF ARRAY FIELD TO ALL PE's
C     FOLLOWED BY RECEIVE OF THE DIFFERENT CONTRIBUTION.
C     
C     EXTERNALS.
C     ----------
C     MPE PACKAGE :
C         MPE_BROADCAST
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

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(NIBLO) :: FIELD
!----------------------------------------------------------------------

      REAL,DIMENSION(MPMAXLENGTH) :: ZCOMBUF

c
      IF(NPROC.EQ.1) THEN
         RETURN
      ELSE 
C     1.1 SEND TO ALL PE'S THE RESPECTIVE CONTRIBUTION 
C         ---------------------------------------------
         KCOUNT=0
         DO IJ=NSTART(IRANK),NEND(IRANK)
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=FIELD(IJ)
         END DO
         CALL MPE_BROADCAST
     1              (ZCOMBUF,MPMAXLENGTH,2,IRANK,ITAG,0,0,0,IERR)
         IF(IERR.LT.0) CALL MPEI_ABORT
     1                    ('MPE_BROADCAST ERROR in MPBCASTSCFLD' )

C     1.2 RECEIVE CONTRIBUTION TO THE FIELD FROM THE OTHER PROCESSES
C         ---------------------------------------------------------- 
         DO IP=1,NPROC-1
            CALL MPE_RECV(ZCOMBUF,MPMAXLENGTH,2,-1,ITAG,0,0,0,
     $                    KRCOUNT,KRFROM,KRTAG,IERR)
            IF(IERR.LT.0) CALL MPEI_ABORT
     1                       ('MPE_RECV ERROR in MPBCASTSCFLD' )
            IF(KRTAG.NE.ITAG) CALL MPEI_ABORT
     1        ('MPE_RECV ERROR in MPBCASTSCFLD MISMATCHED TAGS' )

            KCOUNT=0
            DO IJ=NSTART(KRFROM),NEND(KRFROM)
               KCOUNT=KCOUNT+1
               FIELD(IJ)=ZCOMBUF(KCOUNT)
            END DO 
         ENDDO
      ENDIF
c
      RETURN
      END
