      SUBROUTINE SBOTTOM (F, FL, IJS, IJL, IG, SL)

CSHALLOW
C ----------------------------------------------------------------------
C
C**** *SBOTTOM* - COMPUTATION OF BOTTOM FRICTION.
C
C     G.J.KOMEN AND Q.D.GAO
C     OPTIMIZED BY L.F. ZAMBRESKY
C     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
C
C*    PURPOSE.
C     --------
C
C       COMPUTATION OF BOTTOM FRICTION DISSIPATION
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *SBOTTOM (F, FL, IJS, IJL, IG, SL)*
C          *F*   - SPECTRUM.
C          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
C          *IJS* - INDEX OF FIRST GRIDPOINT
C          *IJL* - INDEX OF LAST GRIDPOINT
C          *IG*  - BLOCK NUMBER
C          *SL*  - TOTAL SOURCE FUNCTION ARRAY
C
C     METHOD.
C     -------
C
C       SEE REFERENCES.
C
C     REFERENCES.
C     -----------
C
C       HASSELMANN ET AL, D. HYDR. Z SUPPL A12(1973) (JONSWAP)
C       BOUWS AND KOMEN, JPO 13(1983)1653-1658
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comshal.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
      PARAMETER (CONST = -2.0*0.038/G)
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F,FL,SL
      DIMENSION SBO(IJS:IJL)
C
C ----------------------------------------------------------------------
C
      DO 1050 M=1,NFRE
         DO 1010 IJ=IJS,IJL
            ARG = 2.* DEPTH(IJ,IG)*TFAK(INDEP(IJ),M)
            ARG = MIN(ARG,50.)
            SBO(IJ) = CONST*TFAK(INDEP(IJ),M)/SINH(ARG)
 1010    CONTINUE

         DO 1040 K=1,NANG
            DO 1030 IJ=IJS,IJL
               SL(IJ,K,M) = SL(IJ,K,M)+SBO(IJ)*F(IJ,K,M)
               FL(IJ,K,M) = FL(IJ,K,M)+SBO(IJ)
 1030      CONTINUE
 1040    CONTINUE
 1050 CONTINUE

      RETURN
      END
