      SUBROUTINE STRESSO (F,IJS, IJL, IG,THWNEW,USNEW,Z0NEW,TAUW,SL)

C ----------------------------------------------------------------------
C
C**** *STRESSO* - COMPUTATION OF WAVE STRESS.
C
C     H. GUNTHER      GKSS/ECMWF  NOVEMBER  1989 CODE MOVED FROM SINPUT.
C     P.A.E.M. JANSSEN      KNMI  AUGUST    1990
C     J. BIDLOT             ECMWF FEBRUARY  1996-97
C
C*    PURPOSE.
C     --------
C
C       COMPUTE NORMALIZED WAVE STRESS FROM INPUT SOURCE FUNCTION
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *STRESSO (F,IJS, IJL, IG,THWNEW,USNEW,Z0NEW,TAUW,SL)*
C          *F*   - WAVE SPECTRUM.
C          *IJS* - INDEX OF FIRST GRIDPOINT.
C          *IJL* - INDEX OF LAST GRIDPOINT.
C          *IG*  - ACTUAL BLOCK NUMBER.
C      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
C                  NOTATION (POINTING ANGLE OF WIND VECTOR,
C                  CLOCKWISE FROM NORTH).
C      *USNEW*     NEW FRICTION VELOCITY IN M/S.
C      *Z0NEW*     ROUGHNESS LENGTH IN M.
C      *TAUW*      WAVE STRESS IN (M/S)**2
C      *SL*        REAL      TOTAL SOURCE FUNCTION ARRAY.
C
C     METHOD.
C     -------
C
C       THE INPUT SOURCE FUNCTION IS INTEGRATED OVER FREQUENCY
C       AND DIRECTIONS.
C       BECAUSE ARRAY *SL* IS USED, ONLY THE INPUT SOURCE
C       HAS TO BE STORED IN *SL* (CALL FIRST SINPUT, THEN
C       STRESSO, AND THEN THE REST OF THE SOURCE FUNCTIONS)
C
C     EXTERNALS.
C     -----------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       R SNYDER ET AL,1981.
C       G. KOMEN, S. HASSELMANN AND K. HASSELMANN, JPO, 1984.
C       P. JANSSEN, JPO, 1985
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comcoup.h"
C
#include "comfred.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comtabl.h"
C
#include "comwind.h"
C
C ----------------------------------------------------------------------
C
C     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F,SL
      REAL,DIMENSION(NINF:NSUP,NBLO) :: TAUW
      REAL,DIMENSION(NINF:NSUP) :: THWNEW,USNEW,Z0NEW
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
      PARAMETER (ROAIR = 1.225, ROWATER = 1000.)
      PARAMETER (XEPS = ROAIR/ROWATER, XINVEPS = 1./XEPS)
C
      DIMENSION CONSTF(NFRE), TAUHF(IJS:IJL), TEMP(IJS:IJL),
     1          XSTRESS(IJS:IJL), YSTRESS(IJS:IJL)
C ----------------------------------------------------------------------
C
C*    1. PRECOMPUTE FREQUENCY SCALING.
C        -----------------------------
C
      CONST0  = DELTH*(ZPI)**4*FR(NFRE)**5/G**2
C
      DO 1000 M=1,NFRE
         CONSTF(M) =ZPI*XINVEPS*FR(M)*DFIM(M)
 1000 CONTINUE
C
C*    2. COMPUTE WAVE STRESS OF ACTUEL BLOCK.
C        ------------------------------------
C
C*    2.1 PRESET STRESS ARRAYS.
C         ---------------------
C
      DO 2100 IJ=IJS,IJL
         XSTRESS(IJ)=0.
         YSTRESS(IJ)=0.
 2100 CONTINUE
C
C*    2.2 INTEGRATE INPUT SOURCE FUNCTION OVER FREQUENCY AND DIRECTIONS.
C         --------------------------------------------------------------
C
      DO 2200 K=1,NANG
         DO 2210 M=1,NFRE
            CONST1=CONSTF(M)*SINTH(K)
            CONST2=CONSTF(M)*COSTH(K)
            DO 2220 IJ=IJS,IJL
               XSTRESS(IJ)=XSTRESS(IJ)+SL(IJ,K,M)*CONST1
               YSTRESS(IJ)=YSTRESS(IJ)+SL(IJ,K,M)*CONST2
 2220       CONTINUE
 2210    CONTINUE
 2200 CONTINUE
C
C*    2.3 CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
C     ----------------------------------------------------
C
      DO 2300 IJ=IJS,IJL
         TEMP(IJ) = 0.
 2300 CONTINUE
C
      DO 2310 K=1,NANG
         TKD = TH(K)
         DO 2320 IJ=IJS,IJL
            COSW     = MAX(COS(TKD-THWNEW(IJ)),0.)
            TEMP(IJ) = TEMP(IJ)+F(IJ,K,NFRE)*COSW**3
 2320    CONTINUE
 2310 CONTINUE
C
      DO 2330 IJ=IJS,IJL
         UST   = MAX(USNEW(IJ),0.000001)
         XI    = UST / DELUST
         XI    = MIN(FLOAT(IUSTAR),XI)
         I     = MIN (IUSTAR-1, INT(XI))
         I     = MAX (0, I)
         DELI1 = MIN (1. ,XI-FLOAT(I))
         DELI2   = 1. - DELI1
C
         XJ    = (G*Z0NEW(IJ)/UST**2-ALPHA) / DELALP
         XJ    = MIN(FLOAT(IALPHA),XJ)
         J     = MIN (IALPHA-1, INT(XJ))
         J     = MAX (0, J)
         DELJ1 = MIN (1. ,XJ-FLOAT(J))
         DELJ2   = 1. - DELJ1
C
         TAU1 = (TAUHFT(I,J  )*DELI2 + TAUHFT(I+1,J  )*DELI1)*DELJ2
     1        + (TAUHFT(I,J+1)*DELI2 + TAUHFT(I+1,J+1)*DELI1)*DELJ1
C
         TAUHF(IJ) = CONST0*TEMP(IJ)*UST**2*TAU1
 2330 CONTINUE
C
      DO 2340 IJ=IJS,IJL
         XSTRESS(IJ) = XSTRESS(IJ)+TAUHF(IJ)*SIN(THWNEW(IJ))
         YSTRESS(IJ) = YSTRESS(IJ)+TAUHF(IJ)*COS(THWNEW(IJ))
         TAUW(IJ,IG) = SQRT(XSTRESS(IJ)**2+YSTRESS(IJ)**2)
 2340 CONTINUE

      RETURN
      END
