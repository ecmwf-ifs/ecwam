      SUBROUTINE MEANSQS(F, IJS, IJL, USNEW)

C ----------------------------------------------------------------------
C
C**** *MEANSQS* - COMPUTATION OF MEAN SQUARE SLOPE.
C
C     P.A.E.M. JANSSEN
C     J. BIDLOT  ECMWF  FEBRUARY 1996  MESSAGE PASSING
C
C
C*    PURPOSE.
C     --------
C
C       COMPUTE MEAN SQUARE SLOPE AT EACH GRID POINT.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *MEANSQS (F, IJS, IJL, USNEW)*
C              *F*   - SPECTRUM.
C              *IJS* - INDEX OF FIRST GRIDPOINT
C              *IJL* - INDEX OF LAST GRIDPOINT
C              *USNEW*     NEW FRICTION VELOCITY IN M/S.
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
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.
C
#include "param.h"
C
#include "parcons.h"
C
#include "comfred.h"
C
#include "commean.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comshal.h"
C
#include "comstat.h"
C
#include "comwind.h"
C
C ----------------------------------------------------------------------
C
!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F

      REAL,DIMENSION(NINF:NSUP) :: USNEW

C ----------------------------------------------------------------------
C
      PARAMETER (XLAMBDAC=0.005,SURFT=0.000075)
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
         SMEAN(IJ) = EPSMIN
 1001 CONTINUE  

      FS     = FR(NFRE)      
      XKC    = ZPI/XLAMBDAC
      FC     = SQRT(G*XKC+SURFT*XKC**3)/ZPI
      CONST2 = ZPI**4*FS**5/G**2
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
            FD = DFIM(M)*(ZPI*FR(M))**4/G**2
            DO 2102 IJ=IJS,IJL
               TEMP2(IJ) = 0.
 2102       CONTINUE
            DO 2103 K=1,NANG
               DO 2104 IJ=IJS,IJL
                  TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
 2104          CONTINUE
 2103       CONTINUE
            DO 2105 IJ=IJS,IJL
               SMEAN(IJ) = SMEAN(IJ)+FD*TEMP2(IJ)
 2105       CONTINUE
 2101    CONTINUE
CSHALLOW
      ELSE
C
C*    2.2 SHALLOW WATER INTEGRATION.
C         --------------------------
C
 2200 CONTINUE
         DO 2202 M=1,NFRE
            DO 2203 IJ=IJS,IJL
               TEMP1(IJ) = DFIM(M)*TFAK(INDEP(IJ),M)**2
               TEMP2(IJ) = 0.
 2203       CONTINUE
            DO 2204 K=1,NANG
               DO 2205 IJ=IJS,IJL
                  TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
 2205          CONTINUE
 2204       CONTINUE
            DO 2206 IJ=IJS,IJL
               SMEAN(IJ) = SMEAN(IJ)+TEMP1(IJ)*TEMP2(IJ)
 2206       CONTINUE
 2202    CONTINUE
      ENDIF
CSHALLOW
C
C*    3. ADD TAIL CORRECTION TO MEAN SQUARE SLOPE.
C        ------------------------------------------
C
 3000 CONTINUE
      DO 3001 IJ=IJS,IJL      
         CONST1    = LOG(FC/FS)  
         CP        = G/(ZPI*FMEAN(IJ))
         XI        = CP/USNEW(IJ)
         ALPHAP    = MAX(0.21/XI,0.0010)
         CONST3    = MAX(CONST2*DELTH*TEMP2(IJ),ALPHAP)
         CONST3    = CONST2*DELTH*TEMP2(IJ)+ALPHAP 
         SMEAN(IJ) = SMEAN(IJ)+CONST1*CONST3
 3001 CONTINUE

      RETURN
      END
