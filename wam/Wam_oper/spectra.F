      SUBROUTINE SPECTRA (IJS, IJL, IG, FL1)

C ----------------------------------------------------------------------
C
C**** *SPECTRA*  - COMPUTATION OF 2-D SPECTRA FOR ONE BLOCK.
C
C     S. HASSELMANN  - JULY 1990
C     H. GUNTHER     - DECEMBER 1990   MODIFIED FOR CYCLE_4.
C     J. BIDLOT      - FEBRUARY 1996   MESSAGE PASSING
C
C*    PURPOSE.
C     --------
C
C       INITIALISATION OF A BLOCK BY 2-D SPECTRA.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *SPECTRA (IJS, IJL, IG, FL1)*
C          *IJS*     INTEGER  FIRST POINT IN BLOCK.
C          *IJL*     INTEGER  LAST  POINT IN BLOCK.
C          *IG*      INTEGER  BLOCK NUMBER.
C          *FL1*      REAL      2-D SPECTRUM FOR EACH GRID POINT
C
C     METHOD.
C     -------
C       1-D JONSWAP SPECTRA AND COSINE**2 SPREADING FUNCTIONS ARE
C       COMPUTED FROM GIVEN WINDS AND PARAMETERS AT EACH GRID POINT.
C       THE 1-D SPECTRA ARE SPREAD OVER THE DIRECTIONS BY MULTIPLICATION
C       WITH THE SPREADING FUNCTION.
C
C     EXTERNALS.
C     ----------
C
C       *SPT*       - COMPUTATION OF COS**2 SPREADING FUNCTION.
C
C
C     REFERENCES.
C     -----------
C
C       K.HASSELMAN,D.B.ROSS,P.MUELLER AND W.SWELL
C          A PARAMETRIC WAVE PREDICTION MODEL
C          JOURNAL OF PHYSICAL OCEANOGRAPHY, VOL. 6, NO. 2, MARCH 1976
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comfred.h"
C
#include "comjons.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comwind.h"
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
      DIMENSION ST(IJS:IJL,NANG), ET(IJS:IJL,NFRE), STH(NANG)
C

!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL1

C ----------------------------------------------------------------------
      PIRHF(F,ALF)=ALF*G**2/(ZPI**4*F**5)
C
C ----------------------------------------------------------------------
C
C*    1. COMPUTE JONSWAP SPECTRUM.
C        -------------------------
C
 1000 CONTINUE
C
      DO 1001 M=1,NFRE
         DO 1002 IJ=IJS,IJL
            IF (FP(IJ).NE.0.) THEN
               FRH = FR(M)
               SIGMA = SA
               IF (FRH.GT.FP(IJ)) SIGMA = SB
               EARG = .5*((FRH-FP(IJ)) / (SIGMA*FP(IJ)))**2
               FJON = 1
               IF (EARG.LT.99.) FJON = GAMMA**EXP(-EARG)
               FMPF = 1.25*(FP(IJ)/FRH)**4
               FJONH = 0.
               IF (FMPF.LT.99.) FJONH = EXP(-FMPF)
               ET(IJ,M) = PIRHF(FRH,ALPHJ(IJ))*FJONH*FJON
            ELSE
               ET(IJ,M) = 0.
            ENDIF
 1002    CONTINUE
 1001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2. COMPUTATION OF SPREADING FUNCTION.
C        ----------------------------------
C
 2000 CONTINUE
C
      DO 2001 IJ=IJS,IJL
         CALL SPR (NANG, NANG, THES(IJ), TH, STH)
         DO 2002 K=1,NANG
            ST(IJ,K) = STH(K)
 2002    CONTINUE
 2001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    3. COMPUTATION OF 2-D SPECTRUM.
C        ----------------------------
C
 3000 CONTINUE
C
      DO 3001 M=1,NFRE
         DO 3002 K=1,NANG
            DO 3003 IJ=IJS,IJL
               FL1(IJ,K,M) = ET(IJ,M)*ST(IJ,K)
 3003       CONTINUE
 3002    CONTINUE
 3001 CONTINUE
C
      RETURN
      END
