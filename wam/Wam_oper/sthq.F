      SUBROUTINE STHQ (F3, IJS, IJL)

C ----------------------------------------------------------------------
C
C**** *STHQ* - COMPUTATION OF MEAN WAVE DIRECTION AT EACH GRID POINT.
C
C     S.D. HASSELMANN
C     OPTIMIZED BY L. ZAMBRESKY
C
C*    PURPOSE.
C     --------
C
C       TO COMPUTE MEAN WAVE DIRECTION AT EACH GRID POINT.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *STHQ (F3, IJS, IJL)*
C          *F3*  - SPECTRUM.
C          *IJS* - INDEX OF FIRST GRIDPOINT
C          *IJL* - INDEX OF LAST GRIDPOINT
C
C     METHOD.
C     -------
C
C       INTEGRATION OF SPECTRUM TIMES SIN AND COS OVER DIRECTION.
C
C     EXTERNALS.
C     ----------
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
#include "param.h"
C
#include "comfred.h"
C
#include "commean.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
C ----------------------------------------------------------------------
C
      DIMENSION TEMP(IJS:IJL)
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
      DIMENSION SI(IJS:IJL), CI(IJS:IJL), F3(NINF-1:NSUP,NANG,NFRE)
C
C ----------------------------------------------------------------------
C
C*    1. INITIALISE SIN AND COS ARRAYS.
C        ------------------------------
C
 1000 CONTINUE
      DO 1001 IJ=IJS,IJL
        SI(IJ) = 0.
        CI(IJ) = 0.
 1001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
C        ------------------------------------------
C
 2000 CONTINUE
      DO 2001 K=1,NANG
         CIH = COSTH(K)
         SIH = SINTH(K)
         DO 2002 IJ=IJS,IJL
            TEMP(IJ) = 0.
 2002    CONTINUE
         DO 2003 M=1,NFRE
            DO 2004 IJ=IJS,IJL
               TEMP(IJ) = TEMP(IJ)+F3(IJ,K,M)*DFIM(M)
 2004       CONTINUE
 2003    CONTINUE
         DO 2005 IJ=IJS,IJL
            SI(IJ) = SI(IJ)+SIH*TEMP(IJ)
            CI(IJ) = CI(IJ)+CIH*TEMP(IJ)
 2005    CONTINUE
 2001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    3. COMPUTE MEAN DIRECTION.
C        -----------------------
C
 3000 CONTINUE
      DO 3001 IJ=IJS,IJL
         IF (CI(IJ).EQ.0.) CI(IJ) = EPSMIN
 3001 CONTINUE
      DO 3002 IJ=IJS,IJL
         THQ(IJ) = ATAN2(SI(IJ),CI(IJ))
 3002 CONTINUE
      DO 3003 IJ=IJS,IJL
         IF (THQ(IJ).LT.0.) THQ(IJ) = THQ(IJ) + ZPI
 3003 CONTINUE

      RETURN
      END
