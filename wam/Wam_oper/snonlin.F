      SUBROUTINE SNONLIN (F, FL, IJS, IJL, IG, SL)

C ----------------------------------------------------------------------
C
C**** *SNONLIN* - COMPUTATION OF NONLINEAR TRANSFER RATE AND ITS
C****             FUNCTIONAL DERIVATIVE (DIAGONAL TERMS ONLY) AND
C****             ADDITION TO CORRESPONDING NET EXPRESSIONS.
C
C     S.D. HASSELMANN.  MPI
C
C     G. KOMEN, P. JANSSEN   KNMI             MODIFIED TO SHALLOW WATER
C     H. GUENTHER, L. ZAMBRESKY               OPTIMIZED
C     H. GUENTHER       GKSS/ECMWF  JUNE 1991 INTERACTIONS BETWEEN DIAG-
C                                             AND PROGNOSTIC PART.
C     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
C
C*    PURPOSE.
C     --------
C
C       SEE ABOVE.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *SNONLIN (F, FL, IJS, IJL, IG, SL)*
C          *F*   - SPECTRUM.
C          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
C          *IJS* - INDEX OF FIRST GRIDPOINT
C          *IJL* - INDEX OF LAST GRIDPOINT
C          *IG*  - BLOCK NUMBER.
C          *SL*  - TOTAL SOURCE FUNCTION ARRAY.
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
#include "comindn.h"
C
#include "commean.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
CSHALLOW
#include "comshal.h"
CSHALLOW
C
#include "comstat.h"
C
C ----------------------------------------------------------------------
C
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F,FL,SL
C
      DIMENSION FTEMP(IJS:IJL), AD(IJS:IJL), DELAD(IJS:IJL), 
     1          DELAP(IJS:IJL), DELAM(IJS:IJL)
CSHALLOW
      DIMENSION ENH(IJS:IJL)
C
C ----------------------------------------------------------------------
C
C*    1. SHALLOW WATER INITIALISATION.
C        -----------------------------
C
      IF (ISHALLO.NE.1) THEN
        DO 1000 IJ=IJS,IJL
          ENH(IJ) = 0.75*DEPTH(IJ,IG)*AKMEAN(IJ)
          ENH(IJ) = MAX(ENH(IJ),.5)
          ENH(IJ) = 1.+(5.5/ENH(IJ))*(1.-.833*ENH(IJ)) *
     1              EXP(-1.25*ENH(IJ))
 1000   CONTINUE
      END IF
CSHALLOW
C
C ----------------------------------------------------------------------
C
C*    2. FREQUENCY LOOP.
C        ---------------
C
      DO 2000 MC=1,NFRE+4
        MP  = IKP (MC)
        MP1 = IKP1(MC)
        MM  = IKM (MC)
        MM1 = IKM1(MC)
        FFACP  = 1.
        FFACP1 = 1.
        FFACM1 = 1.
        FTAIL  = 1.
        IC  = MC
        IP  = MP
        IP1 = MP1
        IM  = MM
        IM1 = MM1
        IF (IP1.GT.NFRE) THEN
          FFACP1 = FRH(IP1-NFRE+1)
          IP1 = NFRE
          IF (IP .GT.NFRE) THEN
            FFACP  = FRH(IP -NFRE+1)
            IP  = NFRE
            IF (IC .GT.NFRE) THEN
              FTAIL  = FRH(IC -NFRE+1)
              IC  = NFRE
              IF (IM1.GT.NFRE) THEN
                FFACM1 = FRH(IM1-NFRE+1)
                IM1 = NFRE
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        FKLAMP  = FKLAP(MC)
        FKLAMP1 = FKLAP1(MC)
        GW2 = FKLAMP1*FFACP*DAL1
        GW1 = GW2*CL11
        GW2 = GW2*ACL1
        GW4 = FKLAMP*FFACP1*DAL1
        GW3 = GW4*CL11
        GW4 = GW4*ACL1
        FKLAMPA = FKLAMP*CL11
        FKLAMPB = FKLAMP*ACL1
        FKLAMP2 = FKLAMP1*ACL1
        FKLAMP1 = FKLAMP1*CL11
        FKLAPA2 = FKLAMPA**2
        FKLAPB2 = FKLAMPB**2
        FKLAP12 = FKLAMP1**2
        FKLAP22 = FKLAMP2**2

        FKLAMM  = FKLAM(MC)
        FKLAMM1 = FKLAM1(MC)
        GW6 = FKLAMM1*DAL2
        GW5 = GW6*CL21
        GW6 = GW6*ACL2
        GW8 = FKLAMM*FFACM1*DAL2
        GW7 = GW8*CL21
        GW8 = GW8*ACL2
        FKLAMMA = FKLAMM*CL21
        FKLAMMB = FKLAMM*ACL2
        FKLAMM2 = FKLAMM1*ACL2
        FKLAMM1 = FKLAMM1*CL21
        FKLAMA2 = FKLAMMA**2
        FKLAMB2 = FKLAMMB**2
        FKLAM12 = FKLAMM1**2
        FKLAM22 = FKLAMM2**2

        IF (ISHALLO.EQ.1) THEN
          DO 2001 IJ=IJS,IJL
            FTEMP(IJ) = AF11(MC)
 2001     CONTINUE
        ELSE
CSHALLOW
          DO 2002 IJ=IJS,IJL
            FTEMP(IJ) = AF11(MC)*ENH(IJ)
 2002     CONTINUE
CSHALLOW
        ENDIF
C
C*    2.1 LOOP FOR ANLULAR SYMMETRY.
C         -------------------------
C
        DO 2100 KH=1,2
C
C*    2.1.1   ANGULAR LOOP.
C             -------------
C
          DO 2110 K=1,NANG
            K1  = K1W (K,KH)
            K2  = K2W (K,KH)
            K11 = K11W(K,KH)
            K21 = K21W(K,KH)
C
C*    2.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND
C*            DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
C             ----------------------------------------------
C
            IF (MC.GT.4) THEN
              DO 2111 IJ=IJS,IJL
                SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )
     1              + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
                SAM = GW5*F(IJ,K2 ,IM ) + GW6*F(IJ,K21,IM )
     1              + GW7*F(IJ,K2 ,IM1) + GW8*F(IJ,K21,IM1)
                FIJ = F(IJ,K  ,IC )*FTAIL
                FAD1 = FIJ*(SAP+SAM)
                FAD2 = FAD1-2.*SAP*SAM
                FAD1 = FAD1+FAD2
                FCEN = FTEMP(IJ)*FIJ
                AD(IJ) = FAD2*FCEN
                DELAD(IJ) = FAD1*FTEMP(IJ)
                DELAP(IJ) = (FIJ-2.*SAM)*DAL1*FCEN
                DELAM(IJ) = (FIJ-2.*SAP)*DAL2*FCEN
 2111         CONTINUE

              DO 2112 IJ=IJS,IJL
                SL(IJ,K2 ,MM ) = SL(IJ,K2 ,MM ) + AD(IJ)*FKLAMM1
                SL(IJ,K21,MM ) = SL(IJ,K21,MM ) + AD(IJ)*FKLAMM2
C
                FL(IJ,K2 ,MM ) = FL(IJ,K2 ,MM ) + DELAM(IJ)*FKLAM12
                FL(IJ,K21,MM ) = FL(IJ,K21,MM ) + DELAM(IJ)*FKLAM22
 2112         CONTINUE

              IF (MM1.LE.NFRE) THEN
                DO 2113 IJ=IJS,IJL
                  SL(IJ,K2 ,MM1) = SL(IJ,K2 ,MM1) + AD(IJ)*FKLAMMA
                  SL(IJ,K21,MM1) = SL(IJ,K21,MM1) + AD(IJ)*FKLAMMB
C
                  FL(IJ,K2 ,MM1) = FL(IJ,K2 ,MM1) + DELAM(IJ)*FKLAMA2
                  FL(IJ,K21,MM1) = FL(IJ,K21,MM1) + DELAM(IJ)*FKLAMB2
 2113           CONTINUE

                IF (MC .LE.NFRE) THEN
                  DO 2114 IJ=IJS,IJL
                    SL(IJ,K  ,MC ) = SL(IJ,K  ,MC ) - 2.*AD(IJ)
C
                    FL(IJ,K  ,MC ) = FL(IJ,K  ,MC ) - 2.*DELAD(IJ)
 2114             CONTINUE

                  IF (MP .LE.NFRE) THEN
                    DO 2115 IJ=IJS,IJL
                      SL(IJ,K1 ,MP ) = SL(IJ,K1 ,MP ) + AD(IJ)*FKLAMP1
                      SL(IJ,K11,MP ) = SL(IJ,K11,MP ) + AD(IJ)*FKLAMP2
C
                      FL(IJ,K1 ,MP ) = FL(IJ,K1 ,MP )
     1                               + DELAP(IJ)*FKLAP12
                      FL(IJ,K11,MP ) = FL(IJ,K11,MP )
     1                               + DELAP(IJ)*FKLAP22
 2115               CONTINUE

                    IF (MP1.LE.NFRE) THEN
                      DO 2116 IJ=IJS,IJL
                        SL(IJ,K1 ,MP1) = SL(IJ,K1 ,MP1)
     1                                 + AD(IJ)*FKLAMPA
                        SL(IJ,K11,MP1) = SL(IJ,K11,MP1)
     1                                 + AD(IJ)*FKLAMPB
C
                        FL(IJ,K1 ,MP1) = FL(IJ,K1 ,MP1)
     1                                 + DELAP(IJ)*FKLAPA2
                        FL(IJ,K11,MP1) = FL(IJ,K11,MP1)
     1                                 + DELAP(IJ)*FKLAPB2
 2116                 CONTINUE
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF

            ELSE

              DO 2117 IJ=IJS,IJL
                SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )
     1              + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
                FIJ = F(IJ,K,IM)
                FAD2 = FIJ*SAP
                FAD1 = 2.*FAD2
                FCEN = FTEMP(IJ)*FIJ
                AD(IJ) = FAD2*FCEN
                DELAD(IJ) = FAD1*FTEMP(IJ)
                DELAP(IJ) = FIJ*DAL1*FCEN
 2117       CONTINUE
C
            DO 2118 IJ=IJS,IJL
                SL(IJ,K  ,MC ) = SL(IJ,K  ,MC ) - 2.*AD(IJ)
                SL(IJ,K1 ,MP ) = SL(IJ,K1 ,MP ) + AD(IJ)*FKLAMP1
                SL(IJ,K11,MP ) = SL(IJ,K11,MP ) + AD(IJ)*FKLAMP2
                SL(IJ,K1 ,MP1) = SL(IJ,K1 ,MP1) + AD(IJ)*FKLAMPA
                SL(IJ,K11,MP1) = SL(IJ,K11,MP1) + AD(IJ)*FKLAMPB
C
                FL(IJ,K  ,MC ) = FL(IJ,K  ,MC ) - 2.*DELAD(IJ)
                FL(IJ,K1 ,MP ) = FL(IJ,K1 ,MP ) + DELAP(IJ)*FKLAP12
                FL(IJ,K11,MP ) = FL(IJ,K11,MP ) + DELAP(IJ)*FKLAP22
                FL(IJ,K1 ,MP1) = FL(IJ,K1 ,MP1) + DELAP(IJ)*FKLAPA2
                FL(IJ,K11,MP1) = FL(IJ,K11,MP1) + DELAP(IJ)*FKLAPB2
 2118         CONTINUE
            ENDIF
C
C*    BRANCH BACK TO 2.1.1 FOR NEXT DIRECTION.
C
 2110     CONTINUE
C
C*    BRANCH BACK TO 2.1 FOR MIRROR INTERACTIONS.
C
 2100   CONTINUE
C
C*    BRANCH BACK TO 2. FOR NEXT FREQUENCY.
C
 2000 CONTINUE

      RETURN
      END
