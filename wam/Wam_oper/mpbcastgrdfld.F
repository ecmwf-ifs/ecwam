      SUBROUTINE MPBCASTGRDFLD(ISEND,ITAG,NXF,NYF,FIELD) 
C
C****  *MPBCASTSCFLD* - BROADCAST 2-D DATA FIELD FROM 1 PE TO ALL OTHERS
C
C     J. BIDLOT    ECMWF   SEPTEMBER 1997
C
C     PURPOSE.
C     --------
C
C     BROADCAST DATA CONTAINED IN 2-D ARRAY FIELD FROM PE ISEND 
C     TO ALL OTHER PE's.
C
C*    INTERFACE.
C     ----------
C
C     CALL *MPGATHERSCFLD*(ITAG,NSTART,NEND,FIELD) 
C
C     *ISEND*     RANK OF THE PROCESS FROM WHICH THE FIELD IS BROADCASTED
C     *ITAG*      TAG ASSOCIATED WITH A PARTICULAR CALL TO SUBROUTINE
C     *NXF*       FIRST DIMENSION OF FIELD
C     *NYF*       SECOND DIMENSION OF FIELD
C     *FIELD*     INPUT/OUTPUT ARRAY CONTAINING THE FIELD
C
C     METHOD.
C     -------
C     MPE BROADCAST OF ARRAY FIELD FROM ISEND
C     AND RECEIVE OF THE CONTRIBUTION ON THE OTHER PE's.
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

C----------------------------------------------------------------------

      REAL,DIMENSION(NXF,NYF) :: FIELD

      REAL,ALLOCATABLE :: ZCOMBUF(:)
      MPLENGTH=NXF*NYF

c
      IF((ISEND.EQ.0).OR.(NPROC.EQ.1)) THEN
         RETURN
      ELSE IF(IRANK.EQ.ISEND) THEN
C     1.1 SEND TO ALL PROCESSORS OTHER THAN ISEND
C         ----------------------------------------
         ALLOCATE(ZCOMBUF(MPLENGTH))

         KCOUNT=0
         DO J=1,NYF
            DO I=1,NXF
               KCOUNT=KCOUNT+1
               ZCOMBUF(KCOUNT)=FIELD(I,J)
            ENDDO
         ENDDO

         CALL MPE_BROADCAST
     1              (ZCOMBUF,MPLENGTH,2,IRANK,ITAG,0,0,0,IERR)
         IF(IERR.LT.0) CALL MPEI_ABORT
     1                    ('MPE_BROADCAST ERROR in MPBCASTGRDFLD' )

         DEALLOCATE(ZCOMBUF)

      ELSE

C     1.2 RECEIVE CONTRIBUTION TO THE FIELD FROM PE ISEND
C         ------------------------------------------------
C
         ALLOCATE(ZCOMBUF(MPLENGTH))

         CALL MPE_RECV(ZCOMBUF,MPLENGTH,2,ISEND,ITAG,0,0,0,
     $                 KRCOUNT,KRFROM,KRTAG,IERR)
         IF(IERR.LT.0) CALL MPEI_ABORT
     1                    ('MPE_RECV ERROR in MPBCASTGRDFLD' )
         IF(KRTAG.NE.ITAG) CALL MPEI_ABORT
     1     ('MPE_RECV ERROR in MPBCASTGRDFLD MISMATCHED TAGS' )

         KCOUNT=0
         DO J=1,NYF
            DO I=1,NXF
               KCOUNT=KCOUNT+1
               FIELD(I,J)=ZCOMBUF(KCOUNT)
            ENDDO
         ENDDO

         DEALLOCATE(ZCOMBUF)

      ENDIF
c
      RETURN
      END
