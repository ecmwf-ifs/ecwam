      SUBROUTINE OUTBETA (IJS, IJL, IG, IGL,USNEW,Z0NEW,NSTART,NEND,
     1                    BETAG, NLONW, NLATW)
C
C ----------------------------------------------------------------------
C
C**** *OUTBETA* - OUTPUT OF CHARNOCK CONSTANT FOR ATMOSPHERIC
C****             MODEL IN GRIDDED FORMAT.
C
C     P.JANSSEN      KNMI/ECMWF  JANUARY 1992
!     J.BIDLOT       ECMWF       FEBRUARY 1996  MESSAGE PASSING
C
C
C*    PURPOSE.
C     --------
C
C       OUTPUT OF CHARNOCK CONSTANT FOR ATMOSPHERIC MODEL
C
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *OUTBETA (IJS, IJL, IG, IGL,USNEW,Z0NEW,NSTART,NEND)*
C         *IJS*    - INDEX OF FIRST GRIDPOINT.
C         *IJL*    - INDEX OF LAST GRIDPOINT.
C         *IG*     - BLOCK NUMBER.
C         *USNEW*     NEW FRICTION VELOCITY IN M/S.
C         *Z0NEW*     ROUGHNESS LENGTH IN M.
C         *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C         *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C
C     EXTERNALS.
C     ----------
C
C       *MAKEGRID*  - GENERATES GRIDDED FIELDS.
C       *MPBCASTSCFLD DISTRIBUTE TO ALL PE's THE CONTRIBUTION OF A BLOCK
C                     FILE
C
C     METHOD.
C     -------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.
C
#include "param.h"
C
#include "parcons.h"
C
#include "comcout.h"
C
#include "comice.h"
C
#include "comintp.h"
C
#include "commesspass.h"
C
#include "commpp.h"
C
#include "comwind.h"
C
C
C ----------------------------------------------------------------------
C
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(NINF:NSUP) :: USNEW,Z0NEW

C ----------------------------------------------------------------------
      DIMENSION BETAG(NLONW,NLATW)
C ----------------------------------------------------------------------
C
      DIMENSION BETA(NIBLO)
      CHARACTER TITL*100
C
C ----------------------------------------------------------------------
C
C*    1. COMPUTE CHARNOCK 'CONSTANT' BETA.
C        ---------------------------------
C
      DO IJ = IJS,IJL
         BETA(IJ) = G*Z0NEW(IJ)/(USNEW(IJ)**2 +0.0001)
      ENDDO
C      
      IF (LICERUN) THEN
         DO IJ = IJS,IJL
            ZICE     = FLOAT(ICEMASK(IJ,IG))
            BETA(IJ) = ZICE*BETA(IJ)
         ENDDO
      ENDIF   
C
C*    1.1 STORE PARAMETERS FOR ACTUAL BLOCK IN GRID ARRAYS.
C         -------------------------------------------------
C
      IF(LMESSPASS) THEN

        CALL MPBCASTSCFLD(20,NSTART,NEND,BETA)

        CALL MAKEGRID (BETA , BETAG  , 1, NEND(NPROC), IG)

      ELSE
        CALL MAKEGRID (BETA , BETAG  , IJS, IJL, IG)
      ENDIF
C
C
C*    1.2 SET LANDVALUES TO CHARNOCK CONSTANT.
C         ------------------------------------
C
      IF (IG.EQ.IGL) THEN
        DO I=1,NGX
          DO J=1,NGY
            IF (BETAG(I,J).LE.0.0) THEN
               BETAG(I,J) = 0.018
            END IF
          ENDDO
        ENDDO
      END IF
C
C ----------------------------------------------------------------------
      RETURN
      END
