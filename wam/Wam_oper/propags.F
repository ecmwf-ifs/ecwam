      SUBROUTINE PROPAGS (F1, F3, MIJS, MIJL, IG, SL, FCONST)

C ----------------------------------------------------------------------
C
C**** *PROPAGS* - COMPUTATION OF A PROPAGATION TIME STEP.
C
C     S.D. HASSELMANN.
C     OPTIMIZED BY: L. ZAMBRESKY AND H. GUENTHER
C
C     MODIFIED BY   H. GUNTHER   01/06/90    -   LAND POINTS ARE TAKEN
C                             OUT OF BLOCKS AND REFRACTION INTEGRATION
C                             CORRECTED FOR N-S AND S-N PROPAGATION.
C
C     K.P. HUBBERT                /07/89    -   DEPTH AND CURRENT
C     S. HASSELMANN   MPIFM       /04/90        REFRACTION SHALLOW
C
C     H. GUNTHER   GKSS/ECMWF   17/01/91    -   MODIFIED FOR CYCLE_4
C
C     J. BIDLOT    ECMWF   FEBRUARY 1997    - ADD SL AND FCONST
C
C*    PURPOSE.
C     --------
C
C       COMPUTATION OF A PROPAGATION TIME STEP.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *PROPAGS(F1, F3, MIJS, MIJL, IG, SL, FCONST)*
C          *F1*   - SPECTRUM AT TIME T.
C          *F3*   - SPECTRUM AT TIME T+DELT.
C          *MIJS* - INDEX OF FIRST POINT
C          *MIJL* - INDEX OF LAST POINT
C          *IG*   - BLOCK NUMBER.
C          *SL*   - SIGMA DOT 
C          *FCONST* - SIGMA/SINH(2KD)
C
C     METHOD.
C     -------
C
C       FIRST ORDER FLUX SCHEME.
C
C     EXTERNALS.
C     ----------
C
C       *DOTDC*     - READ DOT TERMS FOR REFRACTION AND SCATTER TABLE.
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
#include "comcurr.h"
C
#include "comfred.h"
C
#include "comgrid.h"
C
#include "commap.h"
C
#include "commpp.h"
C
#include "comrefd.h"
C
#include "comshal.h"
C
#include "comstat.h"
C
#include "comubuf.h"
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F1,F3,SL
      REAL,DIMENSION(NIBLO) :: DELLA0,DCO,DP1,DP2,DPN,DPS,DLE,DLW
      REAL,DIMENSION(0:NIBLO) :: DPH,DLA
      REAL,DIMENSION(NIBLO) :: DTP,DTM,DRGP,DRGM,DTC
      REAL,DIMENSION(0:NIBLO,NFRE):: CGOND
      REAL,DIMENSION(NIBLO,NFRE):: FCONST

      REAL,DIMENSION(NIBLC) ::  DOP,DOM,DRCP,DRCM
      REAL,DIMENSION(NIBLD) ::  DRDP,DRDM
C
C ----------------------------------------------------------------------
C
C*    0. SPECTRUM AT LAND TO ZERO.
C        -------------------------
C
      DO M=1,NFRE
        DO K=1,NANG
          F3(NINF-1,K,M) = 0.
        ENDDO
      ENDDO
C
      DELPRO = FLOAT(IDELPRO)   
C
C*    0.0 GROUP VELOCITIES. (if shallow water on)
C         -----------------
      IF (ISHALLO.NE.1) THEN
        DO 110 M=1,NFRE
          CGOND(NINF-1,M) = TCGOND(NDEPTH,M)
           DO 113 IJ=NINF,IJLT(IG)
               CGOND(IJ,M) = TCGOND(INDEP(IJ),M)
 113       CONTINUE
 110    CONTINUE
      ENDIF

C
C*    0.1 READ REFRACTION DOT TERMS.
C         --------------------------
C
      IF (IREFRA.NE.0)
     1   CALL DOTDC (MIJS, MIJL, IG, IGL, ISHALLO, IREFRA, SL, FCONST)
C
C*    0.2 SPHERICAL OR CARTESIAN GRID?
C         ----------------------------
C
      IF (ICASE.EQ.1) THEN
C
C*    0.2.1 SPHERICAL GRID.
C           ---------------
C
C*    0.2.1.1 COSINE OF LATITUDE.
C             -------------------
C
         DO 211 IJ = NINF,IJLT(IG)
            JH = KXLT(IJ,IG)
            DCO(IJ) = 1./COSPH(JH)
            DELLA0(IJ) = DELPRO/DELLAM(JH)
  211    CONTINUE
C
C*    0.2.1.2 COMPUTE COS PHI FACTOR FOR ADJOINING GRID POINT.
C             ------------------------------------------------
C
         DO 212 IJ = MIJS,MIJL
            JH = KLAT(IJ,1)
            IF (JH.EQ.NINF-1) THEN
               DP1(IJ) = 1.
            ELSE
               DP1(IJ) = DCO(IJ)/DCO(JH)
            ENDIF
            JH = KLAT(IJ,2)
            IF (JH.EQ.NINF-1) THEN
               DP2(IJ) = 1.
            ELSE
               DP2(IJ) = DCO(IJ)/DCO(JH)
            ENDIF
  212    CONTINUE
         IF (IREFRA.NE.2) THEN
C
C*       BRANCH TO 3. IF WITHOUT REFRACTION OR DEPTH.
C        --------------------------------------------
C
            GOTO 3000
         ELSE
C
C*       BRANCH TO 4. IF DEPTH AND CURRENT REFRACTION.
C        ---------------------------------------------
C
            GOTO 4000
         ENDIF
      ELSE
C
C*    0.2.2 CARTESIAN GRID.
C           ---------------
C
         DO 221 IJ = NINF,IJLT(IG)
            JH = KXLT(IJ,IG)
            DELLA0(IJ) = DELPRO/DELLAM(JH)
  221    CONTINUE
C
C*    0.2.2.1 BRANCH TO 2. IF DEPTH AND CURRENT REFRACTION.
C             ---------------------------------------------
C
         IF (IREFRA.EQ.2) GOTO 2000
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    1. PROPAGATION FOR CARTESIAN GRID
C*       WITHOUT REFRACTION OR DEPTH REFRATION.
C        --------------------------------------
C
 1000 CONTINUE
C
      DELPH0 = DELPRO/DELPHI
      DELTH0 = 0.25*DELPRO/DELTH
C
C*    1.1 LOOP OVER DIRECTIONS.
C         ---------------------
C
      DO 1100 K=1,NANG
         SD = SINTH(K)
         CD = COSTH(K)
C
C*    1.1.1 INDEX FOR ADJOINING POINTS.
C           ---------------------------
C
         IF (SD.LT.0) THEN
            IJLA = 2
         ELSE
            IJLA = 1
         ENDIF
         IF (CD.LT.0) THEN
            IJPH = 2
         ELSE
            IJPH = 1
         ENDIF
C
         IF (ISHALLO.EQ.1) THEN
C
C*    1.1.2 DEEP WATER.
C           -----------
C
            SD = ABS(SD)
            CD = ABS(CD)
C
C*    1.1.2.1 LOOP OVER FREQUENCIES.
C             ----------------------
C
            DO 1120 M=1,NFRE
C
C*    1.1.2.1.1 LOOP OVER GRIDPOINTS.
C               ---------------------
C
               DNO = CD*GOM(M)
               DEA = SD*GOM(M)
               DO 1121 IJ = MIJS,MIJL
                  DTT = 1.-(ABS(SD)*DELLA0(IJ)+ABS(CD)*DELPH0)*GOM(M)
                  DNO = DNO*DELPH0
                  DEA = DEA*DELLA0(IJ)
                  F3(IJ,K,M) = DTT * F1(IJ,K,M )
     1                       + DNO * F1(KLAT(IJ,IJPH),K  ,M)
     2                       + DEA * F1(KLON(IJ,IJLA),K  ,M)
 1121          CONTINUE
C
C*    BRANCH BACK TO 1.1.2.1 FOR NEXT FREQUENCY.
C
 1120       CONTINUE
         ELSE
CSHALLOW
C
C*    1.1.3 SHALLOW WATER.
C           --------------
C
            SD = 0.5*SD
            CD = 0.5*CD
C
C*    1.1.3.1 DEPTH REFRACTION.
C             -----------------
C
            IF(IREFRA.EQ.1) THEN
               KP1 = K+1
               IF (KP1.GT.NANG) KP1 = 1
               KM1 = K-1
               IF (KM1.LT.1) KM1 = NANG
               DO 1131 IJ = MIJS,MIJL
                  DRDP(IJ) = (THDD(IJ,K) + THDD(IJ,KP1))*DELTH0
                  DRDM(IJ) = (THDD(IJ,K) + THDD(IJ,KM1))*DELTH0
 1131          CONTINUE
            ENDIF
C
C*    1.1.3.2 LOOP OVER FREQUENCIES.
C             ----------------------
C
            DO 1130 M=1,NFRE
C
C*    1.1.3.2.2 WEIGHTS IN INTEGRATION SCHEME.
C               ------------------------------
C
               IF (SD.GE.0.) THEN
                  DO 1133 IJ=MIJS,MIJL
                     SDD = SD*DELLA0(IJ)
                     DLA(IJ) = SDD*(CGOND(KLON(IJ,1),M) + CGOND(IJ,M))
                     DTC(IJ) = SDD*(CGOND(KLON(IJ,2),M) + CGOND(IJ,M))
 1133             CONTINUE
               ELSE
                  DO 1134 IJ=MIJS,MIJL
                     SDD = SD*DELLA0(IJ)
                     DLA(IJ) =-SDD*(CGOND(KLON(IJ,2),M) + CGOND(IJ,M))
                     DTC(IJ) =-SDD*(CGOND(KLON(IJ,1),M) + CGOND(IJ,M))
 1134             CONTINUE
               ENDIF

               IF (CD.GE.0.) THEN
                  DO 1135 IJ=MIJS,MIJL
                     CDD = CD*DELPH0
                     DPH(IJ) = CDD*(CGOND(KLAT(IJ,1),M) + CGOND(IJ,M))
                     DTC(IJ) = DTC(IJ)
     1                       + CDD*(CGOND(KLAT(IJ,2),M) + CGOND(IJ,M))
 1135             CONTINUE
               ELSE
                  DO 1136 IJ=MIJS,MIJL
                     CDD = CD*DELPH0
                     DPH(IJ) =-CDD*(CGOND(KLAT(IJ,2),M) + CGOND(IJ,M))
                     DTC(IJ) = DTC(IJ)
     1                        -CDD*(CGOND(KLAT(IJ,1),M) + CGOND(IJ,M))
 1136             CONTINUE
               ENDIF
               IF (IREFRA.EQ.1) THEN
                  DO 1137 IJ = MIJS,MIJL
                     DTHP = FCONST(IJ,M)*DRDP(IJ)
                     DTHM = FCONST(IJ,M)*DRDM(IJ)
                     DTC(IJ) = DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
                     DTP(IJ) = -DTHP+ABS(DTHP)
                     DTM(IJ) =  DTHM+ABS(DTHM)
 1137             CONTINUE
               ENDIF
C
C*    1.1.3.2.3 LOOP OVER GRIDPOINTS.
C               ---------------------
C
               DO 1138 IJ = MIJS,MIJL
                  F3(IJ,K,M) = (1.-DTC(IJ))*F1(IJ,K,M )
     1                       + DPH(IJ) * F1(KLAT(IJ,IJPH),K  ,M)
     2                       + DLA(IJ) * F1(KLON(IJ,IJLA),K  ,M)
 1138          CONTINUE
               IF (IREFRA.EQ.1) THEN
                  DO 1139 IJ = MIJS,MIJL
                     F3(IJ,K,M) = F3(IJ,K,M )
     1                          + DTP(IJ) * F1(IJ,KP1,M)
     2                          + DTM(IJ) * F1(IJ,KM1,M)
 1139             CONTINUE
               ENDIF
C
C*    BRANCH BACK TO 1.1.3.2 FOR NEXT FREQUENCY.
C
 1130       CONTINUE
CSHALLOW
         ENDIF
C
C*    BRANCH BACK TO 1.1 FOR NEXT DIRECTION.
C
 1100 CONTINUE
C
C*    1.2 END OF PROPAGATION FOR CARTESIAN GRID
C*        WITHOUT REFRACTION OR DEPTH REFRACTION, RETURN.
C         -----------------------------------------------
C
      RETURN
C
C ----------------------------------------------------------------------
C
C*    2. PROPAGATION FOR CARTESIAN GRID
C*       WITH DEPTH AND CURRENT REFRACTION.
C        ----------------------------------
C
 2000 CONTINUE
C
      DELPH0 = 0.25*DELPRO/DELPHI
      DELTH0 = 0.25*DELPRO/DELTH
      DELLA0 = 0.25*DELLA0
      DELFR0 = 0.25*DELPRO/(0.1*ZPI)
C
C*    2.1 LOOP OVER DIRECTIONS.
C         ---------------------
C
      DO 2100 K=1,NANG
         KP1 = K+1
         IF (KP1.GT.NANG) KP1 = 1
         KM1 = K-1
         IF (KM1.LT.1) KM1 = NANG
         SD = SINTH(K)
         CD = COSTH(K)
C
C*    2.1.1 DEPTH REFRACTION IF SHALLOW WATER.
C           ----------------------------------
C
         IF (ISHALLO.NE.1) THEN
            DO 2101 IJ = MIJS,MIJL
               DRDP(IJ) = (THDD(IJ,K) + THDD(IJ,KP1))*DELTH0
               DRDM(IJ) = (THDD(IJ,K) + THDD(IJ,KM1))*DELTH0
 2101       CONTINUE
         ENDIF
C
C*    2.1.2 CURRENT REFRACTION.
C           -------------------
C
         DO 2102 IJ = MIJS,MIJL
            DRCP(IJ) = (THDC(IJ,K) + THDC(IJ,KP1))*DELTH0
            DRCM(IJ) = (THDC(IJ,K) + THDC(IJ,KM1))*DELTH0
 2102    CONTINUE
C
C*    2.1.3 LOOP OVER FREQUENCIES.
C           ----------------------
C
         DO 2130 M=1,NFRE
            IF (ISHALLO.EQ.1) THEN
C
C*    2.1.3.1 DEEP WATER.
C             -----------
C
               MP1 = MIN(NFRE,M+1)
               MM1 = MAX(1,M-1)
               DFP = PI*2.1*DELFR0
C
C*    2.1.3.1.1 GROUP VELOCITIES.
C               -----------------
C
               CGS = GOM(M)*SD
               CGC = GOM(M)*CD
C
C*    2.1.3.1.2 WEIGHTS IN INTEGRATION SCHEME.
C               ------------------------------
C
               DLA(NINF-1) = CGS*DELLA0(1)
               DPH(NINF-1) = CGC*DELPH0
               DO 2131 IJ=NINF,IJLT(IG)
                  DLA(IJ) = (U(IJ,IG) + CGS)*DELLA0(IJ)
                  DPH(IJ) = (V(IJ,IG) + CGC)*DELPH0
 2131          CONTINUE
               DO 2132 IJ=MIJS,MIJL
                  DLWE = DLA(IJ) + DLA(KLON(IJ,1))
                  DLEA = DLA(IJ) + DLA(KLON(IJ,2))
                  DLE(IJ) = -DLEA+ABS(DLEA)
                  DLW(IJ) =  DLWE+ABS(DLWE)
                  DTC(IJ) =  DLEA+ABS(DLEA)-DLWE+ABS(DLWE)

                  DPSO = DPH(IJ) + DPH(KLAT(IJ,1))
                  DPNO = DPH(IJ) + DPH(KLAT(IJ,2))
                  DPN(IJ) = -DPNO+ABS(DPNO)
                  DPS(IJ) =  DPSO+ABS(DPSO)
                  DTC(IJ) =  DTC(IJ) + DPNO+ABS(DPNO)-DPSO+ABS(DPSO)

                  DTHP = DRCP(IJ)
                  DTHM = DRCM(IJ)
                  DTP(IJ) = -DTHP+ABS(DTHP)
                  DTM(IJ) =  DTHM+ABS(DTHM)
                  DTC(IJ) =  DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)

                  DTHP    = SL(IJ,K,NFRE) * DFP
                  DTC(IJ) = DTC(IJ) + 2.* ABS(DTHP)
                  DOP(IJ) = (-DTHP+ABS(DTHP))/1.1
                  DOM(IJ) = ( DTHP+ABS(DTHP))*1.1
 2132          CONTINUE
            ELSE
CSHALLOW
C
C*    2.1.3.2 SHALLOW WATER.
C             --------------
C
               MP1 = MIN(NFRE,M+1)
               MM1 = MAX(1,M-1)
               DFP = DELFR0/FR(M)
               DFM = DELFR0/FR(MM1)
C
C
C*    2.1.3.2.2 WEIGHTS IN INTEGRATION SCHEME.
C               ------------------------------
C
               DLA(NINF-1) = SD*CGOND(NINF-1,M)*DELLA0(1)
               DPH(NINF-1) = CD*CGOND(NINF-1,M)*DELPH0
               DO 2134 IJ=NINF,IJLT(IG)
                  DLA(IJ) = (U(IJ,IG) + SD*CGOND(IJ,M))*DELLA0(IJ)
                  DPH(IJ) = (V(IJ,IG) + CD*CGOND(IJ,M))*DELPH0
 2134          CONTINUE
               DO 2135 IJ=MIJS,MIJL
                  DLWE = DLA(IJ) + DLA(KLON(IJ,1))
                  DLEA = DLA(IJ) + DLA(KLON(IJ,2))
                  DLE(IJ) = -DLEA+ABS(DLEA)
                  DLW(IJ) =  DLWE+ABS(DLWE)
                  DTC(IJ) = DLEA+ABS(DLEA)-DLWE+ABS(DLWE)

                  DPSO = DPH(IJ) + DPH(KLAT(IJ,1))
                  DPNO = DPH(IJ) + DPH(KLAT(IJ,2))
                  DPN(IJ) = -DPNO+ABS(DPNO)
                  DPS(IJ) =  DPSO+ABS(DPSO)
                  DTC(IJ) = DTC(IJ) + DPNO+ABS(DPNO)-DPSO+ABS(DPSO)

                  DTHP = FCONST(IJ,M)*DRDP(IJ) + DRCP(IJ)
                  DTHM = FCONST(IJ,M)*DRDM(IJ) + DRCM(IJ)
                  DTC(IJ) = DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
                  DTP(IJ) = -DTHP+ABS(DTHP)
                  DTM(IJ) =  DTHM+ABS(DTHM)

                  DTHP = (SL(IJ,K,M) + SL(IJ,K,MP1))*DFP
                  DTHM = (SL(IJ,K,M) + SL(IJ,K,MM1))*DFM
                  DTC(IJ) =  DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
                  DOP(IJ) = (-DTHP+ABS(DTHP))/1.1
                  DOM(IJ) = ( DTHM+ABS(DTHM))*1.1
 2135          CONTINUE
CSHALLOW
            ENDIF
C
C*    2.1.3.3 LOOP OVER GRIDPOINTS.
C             ---------------------
C
            DO 2136 IJ = MIJS,MIJL
               F3(IJ,K,M) = (1.-DTC(IJ))*F1(IJ,K,M )
     1                    + DPN(IJ) * F1(KLAT(IJ,2),K  ,M)
     2                    + DPS(IJ) * F1(KLAT(IJ,1),K  ,M)
     3                    + DLE(IJ) * F1(KLON(IJ,2),K  ,M)
     4                    + DLW(IJ) * F1(KLON(IJ,1),K  ,M)
     5                    + DTP(IJ) * F1(IJ        ,KP1,M)
     6                    + DTM(IJ) * F1(IJ        ,KM1,M)
     7                    + DOP(IJ) * F1(IJ        ,K  ,MP1)
     8                    + DOM(IJ) * F1(IJ        ,K  ,MM1)
 2136       CONTINUE
C
C*    BRANCH BACK TO 2.1.3 FOR NEXT FREQUENCY.
C
 2130    CONTINUE
C
C*    BRANCH BACK TO 2.1 FOR NEXT DIRECTION.
C
 2100 CONTINUE
C
C*    2.2 END OF PROPAGATION FOR CARTESIAN GRID
C*        WITH DEPTH AND CURRENT REFRACTION, RETURN.
C         ------------------------------------------
C
      RETURN
C
C ----------------------------------------------------------------------
C
C*    3. PROPAGATION FOR SPHERICAL LATITUDE/LONGITUDE GRID
C*       WITHOUT OR DEPTH REFRACTION.
C        -------------------------------------------------
C
 3000 CONTINUE
C
      DELTH0 = 0.25*DELPRO/DELTH
      DELPH0 = 0.5*DELPRO/DELPHI
      IF (ISHALLO.NE.1) THEN
         DO 3010 IJ = 1,IJLT(IG)
            DELLA0(IJ) = 0.5*DELLA0(IJ)
 3010    CONTINUE
      ENDIF
C
C*    3.1 LOOP OVER DIRECTIONS.
C         ---------------------
C
      DO 3100 K=1,NANG
         KP1 = K+1
         IF (KP1.GT.NANG) KP1 = 1
         KM1 = K-1
         IF (KM1.LT.1) KM1 = NANG
         SD = SINTH(K)
         CD = COSTH(K)
         SDA = ABS(SD)
         CDA = ABS(CD)
C
C*    3.1.1 COMPUTE GRID REFRACTION.
C           ------------------------
C
         SP  = DELTH0*(SINTH(K)+SINTH(KP1))/R
         SM  = DELTH0*(SINTH(K)+SINTH(KM1))/R
         DO 3101 IJ = MIJS,MIJL
            JH = KXLT(IJ,IG)
            TANPH = SINPH(JH)*DCO(IJ)
            DRGP(IJ) = TANPH*SP
            DRGM(IJ) = TANPH*SM
 3101    CONTINUE
C
C*    3.1.2 INDEX FOR ADJOINING POINTS.
C           ---------------------------
C
         IF (SD.LT.0) THEN
            IJLA = 2
         ELSE
            IJLA = 1
         ENDIF
         IF (CD.LT.0) THEN
            IJPH = 2
         ELSE
            IJPH = 1
         ENDIF
C
         IF (ISHALLO.EQ.1) THEN
C
C*    3.1.3 DEEP WATER.
C           -----------
C
C*    3.1.3.1 LAT / LONG WEIGHTS IN INTEGRATION SCHEME.
C             -----------------------------------------
C
            DO 3131 IJ=MIJS,MIJL
               DLE(IJ) = DCO(IJ)*SDA*DELLA0(IJ)
 3131       CONTINUE
            IF (CD.GT.0.) THEN
               DO 3132 IJ=MIJS,MIJL
                  DTC(IJ) = DLE(IJ) + DELPH0*CDA*(DP2(IJ) + 1.)
                  DPN(IJ) = DELPH0*CDA*(DP1(IJ) + 1.)
 3132          CONTINUE
            ELSE
               DO 3133 IJ=MIJS,MIJL
                  DTC(IJ) = DLE(IJ) + DELPH0*CDA*(DP1(IJ) + 1.)
                  DPN(IJ) = DELPH0*CDA*(DP2(IJ) + 1.)
 3133          CONTINUE
            ENDIF
C
C*    3.1.3.2 REFRACTION WEIGHTS IN INTEGRATION SCHEME.
C             -----------------------------------------
C
            DO 3134 IJ=MIJS,MIJL
               DTHP = DRGP(IJ)
               DTHM = DRGM(IJ)
               DTC(IJ) = DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
               DTP(IJ) = -DTHP+ABS(DTHP)
               DTM(IJ) =  DTHM+ABS(DTHM)
 3134       CONTINUE
C
C*    3.1.3.3 LOOP OVER FREQUENCIES.
C             ----------------------
C
            DO 3135 M=1,NFRE
C
C*    3.1.3.3.1 LOOP OVER GRIDPOINTS.
C               ---------------------
C
               DO 3136 IJ = MIJS,MIJL
                  DTT = 1. - DTC(IJ)*GOM(M)
                  F3(IJ,K,M) = DTT*F1(IJ,K,M ) + GOM(M) *
     1                        (DPN(IJ) * F1(KLAT(IJ,IJPH),K  ,M)
     2                       + DLE(IJ) * F1(KLON(IJ,IJLA),K  ,M)
     3                       + DTP(IJ) * F1(IJ           ,KP1,M)
     4                       + DTM(IJ) * F1(IJ           ,KM1,M))
 3136          CONTINUE
C
C*    BRANCH BACK TO 3.1.3.3 FOR NEXT FREQUENCY.
C
 3135       CONTINUE
         ELSE
CSHALLOW
C
C*    3.1.4 SHALLOW WATER.
C           --------------
C
C
C*    3.1.4.1 COMPUTE DEPTH REFRACTION.
C             -------------------------
C
         IF (IREFRA.EQ.1) THEN
            DO 3141 IJ = MIJS,MIJL
               DRDP(IJ) = (THDD(IJ,K) + THDD(IJ,KP1))*DELTH0
               DRDM(IJ) = (THDD(IJ,K) + THDD(IJ,KM1))*DELTH0
 3141       CONTINUE
         ENDIF
C
C*    3.1.4.2 LOOP OVER FREQUENCIES.
C             ----------------------
C
            DO 3142 M=1,NFRE
C
C
C*    3.1.4.3.2 LAT / LONG WEIGHTS IN INTEGRATION SCHEME.
C               -----------------------------------------
C
               IF (SD.GT.0.) THEN
                  DO 3144 IJ=MIJS,MIJL
                     DTC(IJ) = 1. - DCO(IJ)*SDA*DELLA0(IJ)*
     1                         (CGOND(KLON(IJ,2),M) + CGOND(IJ,M))
                     DLE(IJ) = DCO(IJ)*SDA*DELLA0(IJ)*
     1                         (CGOND(KLON(IJ,1),M) + CGOND(IJ,M))
 3144             CONTINUE
               ELSE
                  DO 3145 IJ=MIJS,MIJL
                     DTC(IJ) = 1. - DCO(IJ)*SDA*DELLA0(IJ)*
     1                         (CGOND(KLON(IJ,1),M) + CGOND(IJ,M))
                     DLE(IJ) = DCO(IJ)*SDA*DELLA0(IJ)*
     1                         (CGOND(KLON(IJ,2),M) + CGOND(IJ,M))
 3145             CONTINUE
               ENDIF
               IF (CD.GT.0.) THEN
                  DO 3146 IJ=MIJS,MIJL
                     DTC(IJ) = DTC(IJ) - CDA*DELPH0*
     1                      (CGOND(KLAT(IJ,2),M)*DP2(IJ) + CGOND(IJ,M))
                     DPN(IJ) = CDA*DELPH0*
     1                      (CGOND(KLAT(IJ,1),M)*DP1(IJ) + CGOND(IJ,M))
 3146             CONTINUE
               ELSE
                  DO 3147 IJ=MIJS,MIJL
                     DTC(IJ) = DTC(IJ) - CDA*DELPH0*
     1                      (CGOND(KLAT(IJ,1),M)*DP1(IJ) + CGOND(IJ,M))
                     DPN(IJ) = CDA*DELPH0*
     1                      (CGOND(KLAT(IJ,2),M)*DP2(IJ) + CGOND(IJ,M))
 3147             CONTINUE
               ENDIF
C
C*    3.1.4.2.3 REFRACTION WEIGHTS IN INTEGRATION SCHEME.
C               -----------------------------------------
C
               IF (IREFRA.EQ.0) THEN
                  DO 3148 IJ=MIJS,MIJL
                     DTHP = DRGP(IJ)*CGOND(IJ,M)
                     DTHM = DRGM(IJ)*CGOND(IJ,M)
                     DTC(IJ) = DTC(IJ) - DTHP-ABS(DTHP)+DTHM-ABS(DTHM)
                     DTP(IJ) = -DTHP+ABS(DTHP)
                     DTM(IJ) =  DTHM+ABS(DTHM)
 3148             CONTINUE
               ELSE
                  DO 3149 IJ=MIJS,MIJL
                     DTHP = DRGP(IJ)*CGOND(IJ,M)+FCONST(IJ,M)*DRDP(IJ)
                     DTHM = DRGM(IJ)*CGOND(IJ,M)+FCONST(IJ,M)*DRDM(IJ)
                     DTC(IJ) = DTC(IJ) - DTHP-ABS(DTHP)+DTHM-ABS(DTHM)
                     DTP(IJ) = -DTHP+ABS(DTHP)
                     DTM(IJ) =  DTHM+ABS(DTHM)
 3149             CONTINUE
                ENDIF
C
C*    3.1.4.2.4 LOOP OVER GRIDPOINTS.
C               ---------------------
C
               DO 3150 IJ = MIJS,MIJL
                  F3(IJ,K,M) = DTC(IJ)*F1(IJ,K,M )
     1                       + DPN(IJ) * F1(KLAT(IJ,IJPH),K  ,M)
     2                       + DLE(IJ) * F1(KLON(IJ,IJLA),K  ,M)
     3                       + DTP(IJ) * F1(IJ           ,KP1,M)
     4                       + DTM(IJ) * F1(IJ           ,KM1,M)
 3150          CONTINUE
C
C*    BRANCH BACK TO 3.1.4.2 FOR NEXT FREQUENCY.
C
 3142       CONTINUE
CSHALLOW
         ENDIF
C
C*    BRANCH BACK TO 3.1 FOR NEXT DIRECTION.
C
 3100 CONTINUE
C
C*    3.2 END OF PROPAGATION FOR SPHERICAL GRID
C*        WITHOUT REFRACTION OR DEPTH REFRACTION, RETURN.
C         -----------------------------------------------
C
      RETURN
C
C ----------------------------------------------------------------------
C
C*    4. PROPAGATION FOR SPHERICAL LATITUDE/LONGITUDE GRID
C*       WITH DEPTH AND CURRENT REFRACTION.
C        -------------------------------------------------
C
 4000 CONTINUE
C
      DELPH0 = 0.25*DELPRO/DELPHI
      DELTH0 = 0.25*DELPRO/DELTH
      DELLA0 = 0.25*DELLA0
      DELFR0 = 0.25*DELPRO/(0.1*ZPI)
C
C*    4.1 LOOP OVER DIRECTIONS.
C         ---------------------
C
      DO 4100 K=1,NANG
         KP1 = K+1
         IF (KP1.GT.NANG) KP1 = 1
         KM1 = K-1
         IF (KM1.LT.1) KM1 = NANG
         SD = SINTH(K)
         CD = COSTH(K)
C
C*    4.1.1 COMPUTE GRID REFRACTION.
C           ------------------------
C
         SP = DELTH0*(SINTH(K)+SINTH(KP1))/R
         SM = DELTH0*(SINTH(K)+SINTH(KM1))/R
         DO 4111 IJ = MIJS,MIJL
            JH = KXLT(IJ,IG)
            TANPH = SINPH(JH)*DCO(IJ)
            DRGP(IJ) = TANPH*SP
            DRGM(IJ) = TANPH*SM
 4111    CONTINUE
C
C*    4.1.2 COMPUTE DEPTH REFRACTION.
C           -------------------------
C
         IF (ISHALLO.NE.1) THEN
            DO 4121 IJ = MIJS,MIJL
               DRDP(IJ) = (THDD(IJ,K) + THDD(IJ,KP1))*DELTH0
               DRDM(IJ) = (THDD(IJ,K) + THDD(IJ,KM1))*DELTH0
 4121       CONTINUE
         ENDIF
C
C*    4.1.3 COMPUTE CURRENT REFRACTION.
C           ---------------------------
C
         DO 4131 IJ = MIJS,MIJL
            DRCP(IJ) = (THDC(IJ,K) + THDC(IJ,KP1))*DELTH0
            DRCM(IJ) = (THDC(IJ,K) + THDC(IJ,KM1))*DELTH0
 4131    CONTINUE
C
C*    4.1.4 LOOP OVER FREQUENCIES.
C           ----------------------
C
         DO 4140 M=1,NFRE
            MP1 = MIN(NFRE,M+1)
            MM1 = MAX(1,M-1)
            IF (ISHALLO.EQ.1) THEN
C
C*    4.1.4.1 DEEP WATER.
C             -----------
C
C*    4.1.4.1.1 GROUP VELOCITIES.
C               -----------------
C
               DFP = PI*2.1*DELFR0
               CGS = GOM(M)*SD
               CGC = GOM(M)*CD
C
C*    4.1.4.1.2 WEIGHTS IN INTEGRATION SCHEME.
C               ------------------------------
C
               DLA(NINF-1) = CGS*DELLA0(1)
               DPH(NINF-1) = CGC*DELPH0
               DO 4141 IJ=NINF,IJLT(IG)
                  DLA(IJ) = (U(IJ,IG) + CGS)*DELLA0(IJ)*DCO(IJ)
                  DPH(IJ) = (V(IJ,IG) + CGC)*DELPH0
 4141          CONTINUE
               DO 4142 IJ=MIJS,MIJL
                  DLWE = DLA(IJ) + DLA(KLON(IJ,1))
                  DLEA = DLA(IJ) + DLA(KLON(IJ,2))
                  DLE(IJ) = -DLEA+ABS(DLEA)
                  DLW(IJ) =  DLWE+ABS(DLWE)
                  DTC(IJ) =  DLEA+ABS(DLEA)-DLWE+ABS(DLWE)

                  DPSO = DPH(IJ) + DPH(KLAT(IJ,1))*DP1(IJ)
                  DPNO = DPH(IJ) + DPH(KLAT(IJ,2))*DP2(IJ)
                  DPN(IJ) = -DPNO+ABS(DPNO)
                  DPS(IJ) =  DPSO+ABS(DPSO)
                  DTC(IJ) = DTC(IJ) + DPNO+ABS(DPNO)-DPSO+ABS(DPSO)

                  DTHP = DRGP(IJ)*GOM(M) + DRCP(IJ)
                  DTHM = DRGM(IJ)*GOM(M) + DRCM(IJ)
                  DTC(IJ) =  DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
                  DTP(IJ) = -DTHP+ABS(DTHP)
                  DTM(IJ) =  DTHM+ABS(DTHM)

                  DTHP =  SL(IJ,K,NFRE) * DFP
                  DTC(IJ) =  DTC(IJ) + 2. * ABS(DTHP)
                  DOP(IJ) = (-DTHP+ABS(DTHP))/1.1
                  DOM(IJ) = ( DTHP+ABS(DTHP))*1.1
 4142          CONTINUE
            ELSE
CSHALLOW
C
C*    4.1.4.2 SHALLOW WATER.
C             --------------
C
              DFP = DELFR0/FR(M)
              DFM = DELFR0/FR(MM1)
C
C*    4.1.4.2.2 LON/LAT/DIR WEIGHTS IN INTEGRATION SCHEME.
C               ------------------------------------------
C
               DLA(0) = SD*CGOND(0,M)*DELLA0(1)
               DPH(0) = CD*CGOND(0,M)*DELPH0
               DLA(NINF-1) = SD*CGOND(NINF-1,M)*DELLA0(1)
               DPH(NINF-1) = CD*CGOND(NINF-1,M)*DELPH0
               DO 4144 IJ=NINF,IJLT(IG)
                  DLA(IJ) =(U(IJ,IG)+SD*CGOND(IJ,M))*DELLA0(IJ)*DCO(IJ)
                  DPH(IJ) =(V(IJ,IG)+CD*CGOND(IJ,M))*DELPH0
 4144          CONTINUE
               DO 4145 IJ=MIJS,MIJL
                  DLWE = DLA(IJ) + DLA(KLON(IJ,1))
                  DLEA = DLA(IJ) + DLA(KLON(IJ,2))
                  DLE(IJ) = -DLEA+ABS(DLEA)
                  DLW(IJ) =  DLWE+ABS(DLWE)
                  DTC(IJ) =  DLEA+ABS(DLEA)-DLWE+ABS(DLWE)

                  DPSO = DPH(IJ) + DPH(KLAT(IJ,1))*DP1(IJ)
                  DPNO = DPH(IJ) + DPH(KLAT(IJ,2))*DP2(IJ)
                  DPN(IJ) = -DPNO+ABS(DPNO)
                  DPS(IJ) =  DPSO+ABS(DPSO)
                  DTC(IJ) = DTC(IJ) + DPNO+ABS(DPNO)-DPSO+ABS(DPSO)

                  DTHP=DRGP(IJ)*CGOND(IJ,M)
     1                 +FCONST(IJ,M)*DRDP(IJ)+DRCP(IJ)
                  DTHM=DRGM(IJ)*CGOND(IJ,M)
     1                 +FCONST(IJ,M)*DRDM(IJ)+DRCM(IJ)
                  DTC(IJ) =  DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
                  DTP(IJ) = -DTHP+ABS(DTHP)
                  DTM(IJ) =  DTHM+ABS(DTHM)

                  DTHP = (SL(IJ,K,M) + SL(IJ,K,MP1))*DFP
                  DTHM = (SL(IJ,K,M) + SL(IJ,K,MM1))*DFM
                  DTC(IJ) =  DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
                  DOP(IJ) = (-DTHP+ABS(DTHP))/1.1
                  DOM(IJ) = ( DTHM+ABS(DTHM))*1.1
 4145          CONTINUE
CSHALLOW
            ENDIF
C
C*    4.1.4.3 LOOP OVER GRIDPOINTS.
C             ---------------------
C
            DO 4146 IJ = MIJS,MIJL
               F3(IJ,K,M) = (1.-DTC(IJ))*F1(IJ,K,M )
     1                    + DPN(IJ) * F1(KLAT(IJ,2),K  ,M)
     2                    + DPS(IJ) * F1(KLAT(IJ,1),K  ,M)
     3                    + DLE(IJ) * F1(KLON(IJ,2),K  ,M)
     4                    + DLW(IJ) * F1(KLON(IJ,1),K  ,M)
     5                    + DTP(IJ) * F1(IJ        ,KP1,M)
     6                    + DTM(IJ) * F1(IJ        ,KM1,M)
     7                    + DOP(IJ) * F1(IJ        ,K  ,MP1)
     8                    + DOM(IJ) * F1(IJ        ,K  ,MM1)
 4146          CONTINUE
C
C*    BRANCH BACK TO 4.1.4 FOR NEXT FREQUENCY.
C
 4140    CONTINUE
C
C*    BRANCH BACK TO 4.2 FOR NEXT DIRECTION.
C
 4100 CONTINUE
C
C*    4.4 END OF PROPAGATION FOR SPHERICAL GRID
C*        WITH DEPTH AND CURRENT REFRACTION, RETURN.
C         ------------------------------------------
C
      RETURN
      END
