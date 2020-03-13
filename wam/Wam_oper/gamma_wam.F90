C----------------------------------------------------------------------
C
C
      SUBROUTINE COUPLER(FRSTIME)
C
C----------------------------------------------------------------
C
C**** *COUPLER* - COMPUTATION OF LOW- AND HIGH-FREQUENCY STRESS.
C
C     P.A.E.M. JANSSEN
C
C
C
C     PURPOSE.
C     ---------
C
C          COMPUTE LOW AND HIGH-FREQUENCY WAVE STRESS
C          COMPUTE TOTAL STRESS
C
C**   INTERFACE.
C     ----------
C
C          *CALL* *COUPLER*
C
C     METHOD.
C     -------
C
C          SEE REFERENCE FOR WAVE STRESS CALCULATION.
C          A STEADY STATE WIND PROFILE IS ASSUMED.
C          THE WIND STRESS IS COMPUTED USING THE ROUGNNESSLENGTH
C
C                  Z1=Z0/SQRT(1-TAUW/TAU)
C
C           WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-
C           INDUCED STRESS AND TAU IS THE TOTAL STRESS.
C           WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH
C           TAUW/TAU < 1.
C
C     EXTERNALS.
C     ----------
C
C          NONE.
C
C     REFERENCES.
C     -----------
C
C          FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
C
C-------------------------------------------------------------------
C
      PARAMETER (ITOT=$ITOT,JTOT=15)
      PARAMETER(ITAUMAX=${ITAUMAX},JUMAX=${JUMAX})
      PARAMETER(IUSTAR=${IUSTAR},IALPHA=${IALPHA})
      PARAMETER (KL=$KL,ML=$ML,ML1=$ML1,ML54=$ML54)
      LOGICAL FRSTIME,OLD,STATIC
C
C
      COMMON/GRIDPA/SCH(ITOT)
      COMMON/MEANPA/DELTH,DELTHA,DELTH2,DFIM(ML1),DFR(100),TH(KL),
     1              EMEAN(ITOT),FMEAN(ITOT),FPEAK(ITOT),USTOKES(ITOT),
     2              FR(100),G,GZPI28,PI,PIH,ZPI,PIG18,
     3              MWD(ITOT),KLP1,IE10(ITOT) 
      COMMON/TESTS/SWIND(ITOT,ML,KL)
      COMMON/SHALLO/AKMEAN(ITOT),CGOND(ITOT,ML1),FAK(ITOT,ML54),
     1              FRK(ITOT,ML54),DPTH(ITOT)
      COMMON/SPECTR/F1(ITOT,ML,KL),DFP(ITOT,ML,KL)
      COMMON/COUPL/XIMP,XDELF,BETAMAX,ZALP,ALPHA,XKAPPA,ATMWA,XNLEV,
     1             EPS,TAU(ITOT),TAUT(0:ITAUMAX,0:JUMAX),DELTAUW,DELU,
     2             TAUWX(ITOT),TAUWY(ITOT),TAUW(ITOT),ZNUL(ITOT),
     3             Z0UNR(ITOT),Y

      COMMON/WINDPA/THW(ITOT),USTAR(ITOT),UO(ITOT),SUSTAR(ITOT),ITURN
      COMMON/OUTSO/FMEM(ITOT,ML1),SNONL(ITOT,ML1),SINP(ITOT,ML1),
     %             SDIS(ITOT,ML1),STOT(ITOT,ML1),STOT_X(ITOT,ML1),
     %             STOT_Y(ITOT,ML1),SBOT(ITOT,ML1),SBOT_X(ITOT,ML1),
     %             SBOT_Y(ITOT,ML1)
      COMMON/TIMEPA/FP,FPH,IDELT,IDLPRO,IDLSCE,IDLWND,IRDWND,NPROOS,
     1              NTMOWD,NWNDOP
C
      DIMENSION TEMP(ITOT),TEMP1(ITOT,KL),TAUHF(ITOT),UFAC1(ITOT),
     1          UFAC2(ITOT),UFAC3(ITOT),SIG(ITOT),ZNUL_OLD(ITOT),
     2          ALPHAP(ITOT)
      DIMENSION W(JTOT)
      DATA ITERTOT/0/
      DATA NTOT/0/
      SAVE ZNUL_OLD,OMS 
C
C     -----------------------------------------------------------------
C
C
C*         1.     PRELIMINARY CALCULATIONS.
C                 -------------------------
C
       CONST0 = DELTH*ZPI**4*FR(ML1)**5/G**2
       CONST1 = BETAMAX/XKAPPA**2
       OMEGAC = ZPI*FR(ML1)
       X0     = 0.05
       YMAX   = 1.0
       DO 2000 J=1,JTOT
          W(J) = 1.
 2000  CONTINUE
       W(1) = 0.5
       W(JTOT) = 0.5
CC
       DO 1100 IJ=1,ITOT
         USTAR(IJ) = MAX(0.000001,USTAR(IJ))
         SIG(IJ)   = SUSTAR(IJ)
         TEMP(IJ)  = 0.
         TAUWX(IJ) = 0.
         TAUWY(IJ) = 0.
 1100 CONTINUE
C
      ALPHAP = 0.
      DO 1200 K=1,KL
         TKD = (K-1)*DELTH
         DO 1300 IJ=1,ITOT
            COSW        = MAX(COS(TKD-THW(IJ)),0.)
            TEMP1(IJ,K) = COSW
            ALPHA_P     = CONST0*F1(IJ,ML1,K)
            ALPHAP(IJ)  = CONST0*F1(IJ,ML1,K)+ALPHAP(IJ)
            TEMP(IJ)    = TEMP(IJ)+ALPHA_P*COSW**3
 1300    CONTINUE
 1200 CONTINUE
C
C*    THE FIRSTTIME ONLY FRICTION VELOCITY IS CALCULATED WITHOUT
C     EFFECT OF WAVES BECAUSE ZNUL(ITOT) IS UNKNOWN.
C
      IF(FRSTIME) THEN
        DO 1400 IJ=1,ITOT           
          TAUW(IJ) = 0.
 1400   CONTINUE
C
C*      INITIALIZE GRAV.-CAP. PARAMETERS
C
        CALL STRESS_GC(UST,EPS,Z0,ALPHA_P,OMS,Y,XMSS,TAUUNR)
        GO TO 3999
      ENDIF  

C*         2. CALCULATE LOW FREQUENCY PART OF WAVE STRESS.
C          ----------------------------------------------
C
      DO 2100 K=1,KL
         TKD = (K-1)*DELTH
         DO 2300 M=1,ML1
            CM  = G/(ZPI*FR(M))
            FAC = ZPI*FR(M)
            XX1 = FR(M)*DFR(M)            

            DO 2500 IJ=1,ITOT
               DSI        = SWIND(IJ,M,K)
               TAUWX(IJ)  = TAUWX(IJ)+SIN(TKD)*DSI*XX1*2.*ZPI
               TAUWY(IJ)  = TAUWY(IJ)+COS(TKD)*DSI*XX1*2.*ZPI

 2500       CONTINUE
 2300    CONTINUE
 2100 CONTINUE
C
C*         3. CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS
C          --------------------------------------------------
C
      DO 3000 IJ=1,ITOT
        UST      = MAX(USTAR(IJ),0.000001)
        Z0       = ZNUL(IJ)
        OMEGACC  = MAX(OMEGAC,X0*G/UST)
        YC       = MIN(OMEGACC*SQRT(Z0/G),YMAX)
        DELZ     = (LOG(YMAX)-LOG(YC))/REAL(JTOT)

        TAUHF(IJ) = 0.
        DO J=1,JTOT
          Z        = LOG(YC)+REAL(J-1)*DELZ
          Y        = EXP(Z)
          OMEGA    = Y*SQRT(G/Z0)
          IF (OMEGA.LT.OMS) THEN
            CM       = G/OMEGA
            ZX       = UST/CM +ZALP
            ZARG     = MIN(XKAPPA/ZX,20.)
            ZMU      = MIN(G*Z0/CM**2*EXP(ZARG),1.)
C
            ZLOG         = MIN(LOG(ZMU),0.)
            ZBETA        = CONST1*ZMU*ZLOG**4
            TAUHF(IJ)    = TAUHF(IJ)+ZBETA*DELZ*W(J)
          ENDIF
        ENDDO
        TAUHF(IJ) = TEMP(IJ)*USTAR(IJ)**2*TAUHF(IJ)
 3000 CONTINUE
C
C
C
      DO 3400 IJ=1,ITOT
         TX = TAUWX(IJ)
         TY = TAUWY(IJ)
         TAUWLF = SQRT(TX**2+TY**2)

         TAUWY(IJ) = TAUWY(IJ)+TAUHF(IJ)*COS(THW(IJ))
         TAUWX(IJ) = TAUWX(IJ)+TAUHF(IJ)*SIN(THW(IJ))
         TAUW(IJ)  = SQRT(TAUWX(IJ)**2+TAUWY(IJ)**2)
 3400 CONTINUE
C
 3999 CONTINUE

 60   FORMAT(3F10.5) 
C
C***  4.DETERMINE TOTAL STRESS
C       ----------------------
C
C 
      XNUA = 0.00001363
      NITER = 10
      OLD = .FALSE.
      STATIC = .false.

      DO 4000 IJ=1,ITOT
        ZTAUW   = TAUW(IJ)
        UTOP    = UO(IJ)     
        CDRAG   = 1.285*10.**(-3)
        WCD     = SQRT(CDRAG)
        USTOLD  = UTOP*WCD
        TAUOLD  = MAX(USTOLD**2,ATMWA*ZTAUW)
C
        IF (OLD) THEN
          DO 4100 ITER=1,NITER
            X    = MIN(ATMWA*ZTAUW/TAUOLD,0.99)
            UST  = SQRT(TAUOLD)
            Z0   = XNLEV/(EXP(XKAPPA*UTOP/UST)-1.)
            ALPHA_P = ALPHAP(IJ)
            CALL STRESS_GC(UST,EPS,Z0,ALPHA_P,OMS,Y,XMSS,TAUUNR)

            ALPHA  = G*Z0/TAUOLD*SQRT(TAUUNR/TAUOLD)
            ZB     = ALPHA*TAUOLD/G
            ZNU    = 0.1*XNUA/UST
            ZNU1   = ZNU/(1.-X)
            Z0     = 0.5*ZNU1+SQRT(0.25*ZNU1**2+ZB**2/(1.-X))    
            F      = UST-XKAPPA*UTOP/(LOG(XNLEV/Z0))
            XN     = 2.*Z0**2*(1.-X)-ZNU*ZO
            ZFAC   = 2./UST*(3.*ZB**2+0.5*ZNU*Z0-Z0**2)/XN

            DELF   = 1.-XKAPPA*UTOP/(LOG(XNLEV/Z0))**2*ZFAC
      
            USTOLD = UST
            UST    = UST-F/DELF
            DEL = UST-USTOLD
            PCE = 0.01
            IF (ABS(DEL).LT.PCE*UST) GO TO 998

            TAUOLD = MAX(UST**2.,ATMWA*ZTAUW)
 4100     CONTINUE
 998      CONTINUE
        ELSE

          DO 4200 ITER=1,NITER
            UST = USTOLD              
            Z0  = XNLEV/(EXP(XKAPPA*UTOP/UST)-1.)
            ALPHA_P = ALPHAP(IJ)
            
            IF (STATIC) THEN
              ALPHA = 0.007
              ZB     = ALPHA*TAUOLD/G
              TAUUNR = (ZB/Z0)**2*TAUOLD
            ELSE
              CALL STRESS_GC(UST,EPS,Z0,ALPHA_P,OMS,Y,XMSS,TAUUNR)
              ALPHA  = G*Z0/TAUOLD*SQRT(TAUUNR/TAUOLD)
              ZB     = ALPHA*TAUOLD/G
            ENDIF

            TAUV    = 0.04*XNUA*UST/(XKAPPA*Z0)
            TAUNEW  = ATMWA*ZTAUW+TAUV+TAUUNR
            USTNEW  = SQRT(TAUNEW)

            IF (UTOP.LT.4.5) THEN
              W1 = 0.85
            ELSE
              W1 = 0.750
            ENDIF
C 
            UST = W1*USTOLD+(1.-W1)*USTNEW 
C
            DEL = UST-USTOLD
            PCE = 0.005
            IF (ABS(DEL).LT.PCE*UST) GO TO 999
            TAUOLD = UST**2
            USTOLD = UST

 4200     CONTINUE
 999      CONTINUE

        ENDIF
!
!***    2.6 UPDATE FRICTION VELOCITY.
!       -----------------------------
!
        TAU(IJ)  = TAUOLD
C     
        USTAR(IJ) = SQRT(TAU(IJ))
        CDRAG     = USTAR(IJ)**2/UO(IJ)**2
        ZNUL(IJ)  = XNLEV*EXP(-XKAPPA/SQRT(CDRAG))
        Z0UNR(IJ) = ZB

 4000 CONTINUE
C
C     END DO LOOP OVER GRID POINTS
C
C
 10   FORMAT(4F10.4)
C
C     -----------------------------------------------------------------
C
      RETURN
      END
C
C
C----------------------------------------------------------------------
C
      SUBROUTINE STRESS_GC(UST,EPS,Z0,ALPHAP,OMS,Y,XMSS,TAUW)                        
C
C----------------------------------------------------------------------
C
C***  DETERMINE MSS FOR GRAV-CAP WAVES
C
C     AUTHOR: PETER JANSSEN (JULY 1997)
C     ------
C
C     VARIABLE       TYPE         PURPOSE
C     --------       ----         -------
C
C     USTAR          REAL         FRICTION VELOCITY
C     XMSS           REAL         RMS SLOPE GRAV-CAP WAVES
C
C
C     REFERENCES:
C     ----------
C
C     VIERS PAPER EQ.(29)
C
C----------------------------------------------------------------------
C
C
      IMPLICIT NONE
      REAL, ALLOCATABLE :: OMEGA_GC(:), XK_GC(:), DELK_GC(:)
      REAL, ALLOCATABLE :: F1(:)

      REAL X,UST,EPS,Z0,ALPHAP,XMSS,TAUW,PI,KRATIO_GC,ALPHA,
     V     BB,XKS_GC,XKL_GC,XKP,EPS0,SUMT,GAM,SUM,DI,GAM_W,
     V     GAMMA_WAM,ANG,G,SURFT,C,XKAPPA,
     V     BS,BA,ALPHA3,BETA,EPSMIN,XK0,FAC,
     V     ZPI,RAD,DEG,OM,OMS,OMN,Y
   
      INTEGER NWAV_GC,NEND,N,NS,ITER
      LOGICAL FRSTIME
      DATA FRSTIME/.TRUE./
      SAVE PI,ZPI,RAD,DEG,SURFT,G,ALPHA3,ANG,KRATIO_GC,XKS_GC,XKL_GC
      SAVE NWAV_GC,XK0,NS,XK_GC,DELK_GC,OMEGA_GC,F1
C      
C
C*    1.0  DETERMINE GRAV_CAP SPECTRUM, MSS AND TAUWHF.
C          -------------------------------------------
C        
C
      IF (FRSTIME) THEN
!
!*      1.1 DETERMINE CONSTANTS.
!       ------------------------
!
        PI  = 4.*ATAN(1.)
        ZPI = 2.*PI
        RAD = PI/180.
        DEG = 180./PI

        SURFT  = 0.0000717
        G      = 9.806
        ANG    = 0.50
        XKAPPA = 0.4

        FRSTIME=.FALSE.
      ENDIF
!
      ALPHA3 = 3.*ZPI
      XK0 = SQRT(G/SURFT)    
      OMN = SQRT(G*XK0)
      Y = 1./(1.48+2.05*UST)
      XKS_GC = XK0*Y
      NS = LOG(XKS_GC/XK_GC(1))/LOG(KRATIO_GC)+1
      OMS = OMEG_GC(XKS_GC)

      BS  = 0.5*ALPHAP
      EPS0 = ALPHA3*C(XKS_GC)**4/VG_GC(XKS_GC)*BS**2   
      SUM  = 0.
      SUMT = 0.

      DO N=NS,NWAV_GC         

        X     = XK_GC(N)
!     ANALYTICAL FORM INERTIAL SUB RANGE
        BB    = SQRT(EPS0*VG_GC(X)/ALPHA3)/C_GC(X)**2
        F1(N) = X**(-4)*BB
        SUM   = SUM + X**2*F1(N)*X*DELK_GC(N)
        GAM_W = ANG*GAMMA_WAM(OMEGA_GC(N),X,UST,Z0,EPS)
        SUMT  = SUMT+OMEG_GC(X)*GAM_W*BB*DELK_GC(N)/X**3/EPS

      ENDDO
C
      XMSS = SUM
      TAUW = SUMT
C
      RETURN
      END
C
C
C---------------------------------------------------------------------
C
      FUNCTION GAMMA_WAM(OMEGA,XK,USTAR,Z0,EPS)
C
C---------------------------------------------------------------------
C
C**** *GAMMA_WAM* - COMPUTATION OF GROWTHRATE
C
C     P.A.E.M. JANSSEN
C
C
C
C     PURPOSE.
C     ---------
C
C          COMPUTE GROWTHRATE BY WIND
C
C**   INTERFACE.
C     ----------
C
C          *FUNCTION CALL*
C
C     METHOD.
C     -------
C
C          SEE REFERENCE.
C
C     EXTERNALS.
C     ----------
C
C          NONE
C
C     REFERENCES.
C     -----------
C
C          FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
C
C-------------------------------------------------------------------
C
      PARAMETER(BETAMAX=$BETAMAX)
C
C* 1. DETERMINE GROWTH ACCORDING TO JANSSEN-MILES
C     -------------------------------------------
C
      G  = 9.806   
      PI = 4.*ATAN(1.)
      ZPI = 2.*PI
C      BETAMAX = 1.3
      XLAMBDA = 1.0
      ZALP = 0.008
      XKAPPA = 0.40
      ARG_MAX = 50.
      CONST1 = BETAMAX/XKAPPA**2

      USTAR = MAX(0.000001,USTAR)
C
      FAC  = OMEGA
      CM   = FAC/XK
      ZFAK = FAC**2/(G*XK)

      ZX   = USTAR/CM 
      X    = ZX
      ZX   = ZX+ZALP
      X1   = ZX
      ZARG = MIN(XKAPPA/X1,ARG_MAX)

      ZMU    = XK*Z0*EXP(ZARG)
C
      XLOG     = ALOG(ZMU/XLAMBDA)
      ZLOG     = MIN(XLOG,0.0)
      ZBETA    = CONST1*ZMU*ZLOG**4
      CK       = EPS*ZBETA*X**2

      GAMMA_WAM = CK*ZFAK*FAC
C
C     -----------------------------------------------------------------
C
      RETURN
      END
