      SUBROUTINE FEMEAN (F, IJS, IJL)

C ----------------------------------------------------------------------
C
C**** *FEMEAN* - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT.
C
C     S.D. HASSELMANN
C     MODIFIED : P.JANSSEN (INTEGRATION OF F**-4 TAIL)
C     OPTIMIZED BY : L. ZAMBRESKY AND H. GUENTHER
C
C
C*    PURPOSE.
C     --------
C
C       COMPUTE MEAN FREQUENCY AT EACH GRID POINT.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *FEMEAN (F, IJS, IJL)*
C              *F*   - SPECTRUM.
C              *IJS* - INDEX OF FIRST GRIDPOINT
C              *IJL* - INDEX OF LAST GRIDPOINT
C
C     METHOD.
C     -------
C
C       NONE.
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
#include "comstat.h"
C
#include "comshal.h"
C
CSHALLOW
C
C ----------------------------------------------------------------------
C
      DIMENSION F(NINF-1:NSUP,NANG,NFRE)
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
      DIMENSION TEMP1(IJS:IJL), TEMP2(IJS:IJL)
C
C ----------------------------------------------------------------------
C
C*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
C        ------------------------------------------------
C
 1000 CONTINUE
      DO 1001 IJ=IJS,IJL
         FMEAN(IJ) = EPSMIN
 1001 CONTINUE
      DELT25 = 0.20*DELTH
C
C ----------------------------------------------------------------------
C
C*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
C        ------------------------------------------
C
 2000 CONTINUE
      IF (ISHALLO.EQ.1) THEN
C
C*    2.1 DEEP WATER INTEGRATION.
C         -----------------------
C
 2100 CONTINUE
         DO 2101 M=1,NFRE
            FD = DFIM(M)/FR(M)
            DO 2102 IJ=IJS,IJL
               TEMP2(IJ) = 0.
 2102       CONTINUE
            DO 2103 K=1,NANG
               DO 2104 IJ=IJS,IJL
                  TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
 2104          CONTINUE
 2103       CONTINUE
            DO 2105 IJ=IJS,IJL
               FMEAN(IJ) = FMEAN(IJ)+FD*TEMP2(IJ)
 2105       CONTINUE
 2101    CONTINUE
CSHALLOW
      ELSE
C
C*    2.2 SHALLOW WATER INTEGRATION.
C         --------------------------
C
 2200 CONTINUE
         DO 2201 IJ=IJS,IJL
            AKMEAN(IJ) = EPSMIN
 2201    CONTINUE
         DO 2202 M=1,NFRE
            FD=DFIM(M)/FR(M)
            DO 2203 IJ=IJS,IJL
               TEMP1(IJ) = DFIM(M)/SQRT(TFAK(INDEP(IJ),M))
               TEMP2(IJ) = 0.
 2203       CONTINUE
            DO 2204 K=1,NANG
               DO 2205 IJ=IJS,IJL
                  TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
 2205          CONTINUE
 2204       CONTINUE
            DO 2206 IJ=IJS,IJL
               FMEAN(IJ) = FMEAN(IJ)+FD*TEMP2(IJ)
               AKMEAN(IJ) = AKMEAN(IJ)+TEMP1(IJ)*TEMP2(IJ)
 2206       CONTINUE
 2202    CONTINUE
C
C        ADD TAIL TO MEAN WAVENUMBER AND NORMALIZE WITH TOTAL ENERGY.
C
         DEL2 = DELT25*SQRT(G)/ZPI
         DO 2207 IJ=IJS,IJL
            AKMEAN(IJ) = AKMEAN(IJ)+DEL2*TEMP2(IJ)
            AKMEAN(IJ) = (EMEAN(IJ)/AKMEAN(IJ))**2
 2207    CONTINUE
      ENDIF
CSHALLOW
C
C*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
C*       NORMALIZE WITH TOTAL ENERGY.
C        ------------------------------------------
C
 3000 CONTINUE
      DO 3001 IJ=IJS,IJL
         FMEAN(IJ) = FMEAN(IJ)+DELT25*TEMP2(IJ)
         FMEAN(IJ) = EMEAN(IJ)/FMEAN(IJ)
 3001 CONTINUE

      RETURN
      END
