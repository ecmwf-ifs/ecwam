      SUBROUTINE SINPUT (F, FL, IJS, IJL, IG,THWOLD,USOLD,Z0OLD,
     1                   THWNEW,USNEW,Z0NEW,SL)
C ----------------------------------------------------------------------
C
C**** *SINPUT* - COMPUTATION OF INPUT SOURCE FUNCTION.
C
C     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
C
C     OPTIMIZED BY : H. GUENTHER
C
C     MODIFIED BY : J-R BIDLOT NOVEMBER 1995
C                   J-R BIDLOT FEBRUARY 1996-97
C
C*    PURPOSE.
C     ---------
C
C       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
C       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
C       INPUT SOURCE FUNCTION.
C
C**   INTERFACE.
C     ----------
C
C     *CALL* *SINPUT (F, FL, IJS, IJL, IG,THWOLD,USOLD,Z0OLD,
C    1                   THWNEW,USNEW,Z0NEW,SL)
C          *F*   - SPECTRUM.
C          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
C          *IJS* - INDEX OF FIRST GRIDPOINT.
C          *IJL* - INDEX OF LAST GRIDPOINT.
C          *IG*  - BLOCK NUMBER.
C      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
C                  NOTATION (POINTING ANGLE OF WIND VECTOR,
C                  CLOCKWISE FROM NORTH).
C      *THWOLD*    INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
C                  WIND VELOCITY.
C      *USNEW*     NEW FRICTION VELOCITY IN M/S.
C      *USOLD*     INTERMEDIATE STORAGE OF MODULUS OF FRICTION
C                  VELOCITY.
C      *Z0NEW*     ROUGHNESS LENGTH IN M.
C      *Z0OLD*     INTERMEDIATE STORAGE OF ROUGHNESS LENGTH IN
C                  M.
C      *SL*        REAL      TOTAL SOURCE FUNCTION ARRAY.
C
C     METHOD.
C     -------
C
C       SEE REFERENCE.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     MODIFICATION
C     ------------
C
C     REMOVAL OF CALL TO CRAY SPECIFIC FUNCTIONS EXPHF AND ALOGHF
C     BY THEIR STANDARD FORTRAN EQUIVALENT EXP and ALOGHF
C
C     REFERENCE.
C     ----------
C
C       P. JANSSEN, J.P.O., 1989.
C       P. JANSSEN, J.P.O., 1991
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comcoup.h"
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
#include "comwind.h"
C
C ----------------------------------------------------------------------
C
C     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F,FL,SL
      REAL,DIMENSION(NINF:NSUP,NBLO) :: THWOLD,USOLD,Z0OLD
      REAL,DIMENSION(NINF:NSUP) :: THWNEW,USNEW,Z0NEW
C
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
      PARAMETER (ROAIR = 1.225, ROWATER = 1000.)
      PARAMETER (XEPS = ROAIR/ROWATER, XINVEPS = 1./XEPS)
C
      DIMENSION TEMP(IJS:IJL,NANG), TEMP1(IJS:IJL,NANG),
     1          UCO(IJS:IJL), ZCO(IJS:IJL), UCN(IJS:IJL), ZCN(IJS:IJL),
     2          UFAC1(IJS:IJL), UFAC2(IJS:IJL), CM(IJS:IJL)
C
C ----------------------------------------------------------------------
C
C*    1. PRECALCULATED ANGULAR DEPENDENCE.
C        ---------------------------------
C
 1000 CONTINUE
      DO 1001 K=1,NANG
         TKD=TH(K)
         DO 1002 IJ=IJS,IJL
            TEMP (IJ,K) = COS(TKD-THWOLD(IJ,IG))
            TEMP1(IJ,K) = COS(TKD-THWNEW(IJ))
 1002    CONTINUE
 1001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2. LOOP OVER FREQUENCIES.
C        ----------------------
C
      CONST1  = XEPS*BETAMAX/XKAPPA**2
      DO 2001 M=1,NFRE
         FAC = ZPI*FR(M)
         CONST=FAC*CONST1
C
C*      INVERSE OF PHASE VELOCITIES.
C       ----------------------------
C
         IF (ISHALLO.EQ.1) THEN
            DO 2002 IJ=IJS,IJL
               CM(IJ) = FAC/G
 2002       CONTINUE
         ELSE
            DO 2003 IJ=IJS,IJL
               CM(IJ) = TFAK(INDEP(IJ),M)/FAC
 2003       CONTINUE
         END IF
C
C*      PRECALCULATE FREQUENCY DEPENDENCE.
C       ----------------------------------
C
         DO 2004 IJ=IJS,IJL
            UCO(IJ) = USOLD(IJ,IG)*CM(IJ) + ZALP
            ZCO(IJ) = ALOG(G*Z0OLD(IJ,IG)*CM(IJ)**2)
            UCN(IJ) = USNEW(IJ)*CM(IJ) + ZALP
            ZCN(IJ) = ALOG(G*Z0NEW(IJ)*CM(IJ)**2)
 2004    CONTINUE
C
C*    2.1 LOOP OVER DIRECTIONS.
C         ---------------------
C
         DO 2101 K=1,NANG
            DO 2102 IJ=IJS,IJL
               UFAC1(IJ) = 0.
               UFAC2(IJ) = 0.
 2102       CONTINUE
            DO 2103 IJ=IJS,IJL
               IF (TEMP(IJ,K).GT.0.01) THEN
                  X    = TEMP(IJ,K)*UCO(IJ)
                  ZARG = XKAPPA/X
                  ZLOG = ZCO(IJ) + ZARG
                  IF (ZLOG.LT.0.) THEN
                     ZLOG2X=ZLOG*ZLOG*X
                     UFAC1(IJ) = CONST*EXP(ZLOG)*ZLOG2X*ZLOG2X
                  ENDIF
               ENDIF
 2103       CONTINUE
C
            DO 2104 IJ=IJS,IJL
               IF (TEMP1(IJ,K).GT.0.01) THEN
                  X    = TEMP1(IJ,K)*UCN(IJ)
                  ZARG = XKAPPA/X
                  ZLOG = ZCN(IJ) + ZARG
                  IF (ZLOG.LT.0.) THEN
                     ZLOG2X=ZLOG*ZLOG*X
                     UFAC2(IJ) = CONST*EXP(ZLOG)*ZLOG2X*ZLOG2X
                  ENDIF
               ENDIF
 2104       CONTINUE
C
C*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
C         ------------------------------------------------
C
            DO 2201 IJ=IJS,IJL
               SL(IJ,K,M) = 0.5*(UFAC1(IJ)+UFAC2(IJ))*F(IJ,K,M)
               FL(IJ,K,M) = UFAC2(IJ)
 2201       CONTINUE
 2101    CONTINUE
 2001 CONTINUE

      RETURN
      END
