      SUBROUTINE MPCRTBL 

C ----------------------------------------------------------------------
C
C**** *MPCRTBL* -
C
C     J. BIDLOT     ECMWF    MARCH 1996    MESSAGE PASSING
C
C
C*    PURPOSE.
C     --------
C     CREATE A TABLE THAT ASSOCIATES EACH INTEGRATED OUTPUT PARAMETER
C     WITH A SPECIFIC PROCESS. THIS TABLE IS USED TO SPREAD THE RETRIEVAL
C     OF THOSE PARAMETER FIELD IN GRID FORMAT ON AS MANY PROCESSOR AS
C     POSSIBLE
C
C**   INTERFACE.
C     ----------
C
C
C     METHOD.
C     -------
C     SYSTEMATICALY ASSIGN A PROCESS TO A PARAMETER AND START AGAIN WHEN
C     RUN OUT OF PROCESS.
C
C     EXTERNALS.
C     ----------
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comcout.h"
C
#include "comtest.h"
C
#include "comtext.h"
C
#include "commpp.h"
C
#include "comstat.h"
C
C ----------------------------------------------------------------------
C
C
C     1. CREATE TABLE MAPPING OUTPUT INTEGRATED PARAMETER WITH PROCESS RANK
C        ------------------------------------------------------------------

        IR=1
      DO IFLAG=1,JPPFLAG
        IF((IASSI.EQ.1).AND.
     1     ((IFLAG.EQ.1).OR.(IFLAG.EQ.7))) THEN
          IPFGTBL(IFLAG)=IR
          IR=IR+1
          IF(IR.GT.NPROC) IR=1
        ELSE IF(PFLAG(IFLAG).OR.FFLAG(IFLAG).OR.GFLAG(IFLAG)) THEN
          IPFGTBL(IFLAG)=IR
          IR=IR+1
          IF(IR.GT.NPROC) IR=1
        ELSE
          IPFGTBL(IFLAG)=0
        ENDIF
      END DO
C
C FIND P.E. FOR THE GATHERING OF THE WIND STRESS FIELD
C WHENEVER ASSIMILATION IS DONE BUT NOT OUTPUT OF IT  WAS REQUESTED.
C
      IF((IASSI.EQ.1).AND.(IPFGTBL(4).EQ.0)) THEN 
         IFLAG=JPPFLAG
         DO WHILE ((IPFGTBL(IFLAG).EQ.0).AND.(IFLAG.GT.1))
           IFLAG=IFLAG-1
         END DO
         IPFGTBL(4)=IPFGTBL(IFLAG)+1
         IF(IPFGTBL(4).GT.NPROC)IPFGTBL(4)=1
      ENDIF
C
C FIND P.E. FOR OUTPUT OF RESTART FILES
C
      IFLAG=JPPFLAG
cc
ccc set to 1 in order to simplify ouput to FDB
      IPFGTBL(JPPFLAG+1)=1

CCC      DO WHILE ((IPFGTBL(IFLAG).EQ.0).AND.(IFLAG.GT.1))
CCC         IFLAG=IFLAG-1
CCC      END DO
CCC      IPFGTBL(JPPFLAG+1)=IPFGTBL(IFLAG)+1
CCC      IF(IPFGTBL(JPPFLAG+1).GT.NPROC)IPFGTBL(JPPFLAG+1)=1

      END
