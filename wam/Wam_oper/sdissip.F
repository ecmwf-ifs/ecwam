      SUBROUTINE SDISSIP (F, FL, IJS, IJL, SL)

C ----------------------------------------------------------------------
C
C**** *SDISSIP* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.
C
C     S.D.HASSELMANN.
C     MODIFIED TO SHALLOW WATER : G. KOMEN , P. JANSSEN
C     OPTIMIZATION : L. ZAMBRESKY
C     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
C
C*    PURPOSE.
C     --------
C       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
C       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
C       OF DISSIPATION SOURCE FUNCTION.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *SDISSIP (F, FL, IJS, IJL, SL)*
C          *F*   - SPECTRUM.
C          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
C          *IJS* - INDEX OF FIRST GRIDPOINT
C          *IJL* - INDEX OF LAST GRIDPOINT
C          *SL*  - TOTAL SOURCE FUNCTION ARRAY
C
C     METHOD.
C     -------
C
C       SEE REFERENCES.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       G.KOMEN, S. HASSELMANN AND K. HASSELMANN, ON THE EXISTENCE
C          OF A FULLY DEVELOPED WINDSEA SPECTRUM, JGR, 1984.
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
CSHALLOW
#include "comshal.h"
C
CSHALLOW
C
#include "comstat.h"
C
C ----------------------------------------------------------------------
C
      DIMENSION TEMP1(IJS:IJL), SDS(IJS:IJL)
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F,FL,SL
      PARAMETER (CDIS = 4.5)
      PARAMETER (CONSD = -0.5*CDIS*ZPI**9/G**4)
      PARAMETER (CONSS = -0.5*CDIS*ZPI)

C ----------------------------------------------------------------------
C
C*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
C*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
C        --------------------------------------------------------------
C
 2000 CONTINUE
      IF (ISHALLO.EQ.1) THEN
         DO 2001 IJ=IJS,IJL
            SDS(IJ) = CONSD*EMEAN(IJ)**2*FMEAN(IJ)**9
 2001    CONTINUE

         DO 2002 M=1,NFRE
            DO 2003 IJ=IJS,IJL
               X         = (FR(M)/FMEAN(IJ))**2
               TEMP1(IJ) = SDS(IJ)*( X + X**2)
 2003       CONTINUE
            DO 2004 K=1,NANG
               DO 2005 IJ=IJS,IJL
                  SL(IJ,K,M) = SL(IJ,K,M)+TEMP1(IJ)*F(IJ,K,M)
                  FL(IJ,K,M) = FL(IJ,K,M)+TEMP1(IJ)
 2005          CONTINUE
 2004       CONTINUE
 2002    CONTINUE
      ELSE
CSHALLOW
         DO 2101 IJ=IJS,IJL
            SDS(IJ) = CONSS*FMEAN(IJ)*EMEAN(IJ)**2*AKMEAN(IJ)**4
 2101    CONTINUE

         DO 2102 M=1,NFRE
            DO 2103 IJ=IJS,IJL
               X         = TFAK(INDEP(IJ),M)/AKMEAN(IJ)
               TEMP1(IJ) = SDS(IJ)*( X + X**2)
 2103       CONTINUE
            DO 2104 K=1,NANG
               DO 2105 IJ=IJS,IJL
                  SL(IJ,K,M) = SL(IJ,K,M)+TEMP1(IJ)*F(IJ,K,M)
                  FL(IJ,K,M) = FL(IJ,K,M)+TEMP1(IJ)
 2105          CONTINUE
 2104       CONTINUE
 2102    CONTINUE
CSHALLOW
      ENDIF
C
      RETURN
      END
