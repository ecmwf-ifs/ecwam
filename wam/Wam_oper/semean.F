      SUBROUTINE SEMEAN (F3, IJS, IJL)

C ----------------------------------------------------------------------
C
C**** *SEMEAN* - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.
C
C     S.D. HASSELMANN.
C     OPTIMIZED BY: L. ZAMBRESKY AND H. GUENTHER
C
C*    PURPOSE.
C     --------
C
C       TO COMPUTE TOTAL ENERGY AT EACH GRID POINT.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *SEMEAN(F3, IJS, IJL)*
C          *F3*  - SPECTRUM.
C          *IJS* - INDEX OF FIRST GRIDPOINT
C          *IJL* - INDEX OF LAST GRIDPOINT
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
C ----------------------------------------------------------------------
C
      DIMENSION F3(NINF-1:NSUP,NANG,NFRE), TEMP(IJS:IJL)
C
#include "parcons.h"
C
C ----------------------------------------------------------------------
C
C*    1. INITIALISE ENERGY ARRAY.
C        ------------------------
C
 1000 CONTINUE
      DO 1001 IJ=IJS,IJL
         EMEAN(IJ) = EPSMIN
 1001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2. INTEGRATE OVER FREQUENCIES AND DIRECTION.
C        -----------------------------------------
C
 2000 CONTINUE
      DO 2001 M=1,NFRE
         DO 2002 IJ=IJS,IJL
            TEMP(IJ) = 0.
 2002   CONTINUE
        DO 2003 K=1,NANG
           DO 2004 IJ=IJS,IJL
              TEMP(IJ) = TEMP(IJ)+F3(IJ,K,M)
 2004      CONTINUE
 2003   CONTINUE
        DO 2005 IJ=IJS,IJL
           EMEAN(IJ) = EMEAN(IJ)+DFIM(M)*TEMP(IJ)
 2005   CONTINUE
 2001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    3. ADD TAIL ENERGY.
C        ----------------
C
 3000 CONTINUE
      DELT25 = FR(NFRE)/4.*DELTH
      DO 3001 IJ=IJS,IJL
         EMEAN(IJ) = EMEAN(IJ)+DELT25*TEMP(IJ)
 3001 CONTINUE

      RETURN
      END
