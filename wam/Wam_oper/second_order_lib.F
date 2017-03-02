!
!
!-----------------------------------------------------------------------
!
!***  *REAL FUNCTION* *A(XI,XJ,THI,THJ)
!
!-----------------------------------------------------------------------
      REAL FUNCTION A(XI,XJ,THI,THJ)
!
!***  *A*  DETERMINES THE MINUS INTERACTIONS.
!
!     PETER JANSSEN
!
!     PURPOSE.
!     --------
!
!              GIVES NONLINEAR TRANSFER COEFFICIENT FOR THREE
!              WAVE INTERACTIONS OF GRAVITY WAVES IN THE
!              IDEAL CASE OF NO CURRENT. (CF.ZAKHAROV)
!
!     INTERFACE.
!     ----------
!              *A(XI,XJ)*
!                      *XI*  - WAVE NUMBER
!                      *XJ*  - WAVE NUMBER
!     METHOD.
!     -------
!              NONE
!
!     EXTERNALS.
!     ----------
!              NONE.
!
!-----------------------------------------------------------------------
!
      USE YOWPCONS, ONLY : G, PI
      IMPLICIT NONE
      REAL RI,RJ,RK,XI,XJ,THI,THJ,THK,OI,OJ,OK,FI,FJ,FK,VABS,VDIR,OMEG,
     V     A1,A3
!
!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------
!

      RI = XI
      RJ = XJ
      RK  = VABS(RI,RJ,THI,THJ)
      THK = VDIR(RI,RJ,THI,THJ)
        
      OI=OMEG(RI)
      OJ=OMEG(RJ)
      OK=OMEG(RK)

      FI = SQRT(OI/(2.*G))
      FJ = SQRT(OJ/(2.*G))
      FK = SQRT(OK/(2.*G))
 
      
      A = FK/(FI*FJ)*(A1(RK,RI,RJ,THK,THI,THJ)+
     V                A3(RK,RI,RJ,THK-PI,THI,THJ))
      
      RETURN
      END
!
!***  *REAL FUNCTION* *B(XI,XJ,THI,THJ)
!
!-----------------------------------------------------------------------
      REAL FUNCTION B(XI,XJ,THI,THJ)
!
!***  *B*  DETERMINES THE PLUS INTERACTION COEFFICIENTS.
!
!     PETER JANSSEN
!
!     PURPOSE.
!     --------
!
!              GIVES NONLINEAR TRANSFER COEFFICIENT FOR THREE
!              WAVE INTERACTIONS OF GRAVITY WAVES IN THE
!              IDEAL CASE OF NO CURRENT. (CF.ZAKHAROV)
!
!     INTERFACE.
!     ----------
!              *B(XI,XJ)*
!                      *XI*  - WAVE NUMBER
!                      *XJ*  - WAVE NUMBER
!     METHOD.
!     -------
!              NONE
!
!     EXTERNALS.
!     ----------
!              NONE.
!
!-----------------------------------------------------------------------
!
      USE YOWPCONS, ONLY : G, PI
      IMPLICIT NONE
      REAL DEL,RI,RJ,RK,XI,XJ,THI,THJ,THK,OI,OJ,OK,FI,FJ,FK,VABS,VDIR,
     V     OMEG,A2
!
!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------
!
      DEL = 0.
      RI = XI
      RJ = XJ
      RK  = VABS(RJ,RI,THJ,THI-PI)
      THK = VDIR(RJ,RI,THJ,THI-PI)      

      OI=OMEG(RI)+DEL
      OJ=OMEG(RJ)+DEL
      OK=OMEG(RK)+DEL

      FI = SQRT(OI/(2.*G))
      FJ = SQRT(OJ/(2.*G))
      FK = SQRT(OK/(2.*G))

      B = 0.5*FK/(FI*FJ)*(A2(RK,RI,RJ,THK,THI,THJ)+
     V                    A2(RK,RJ,RI,THK-PI,THJ,THI))

      RETURN
      END
!
!-----------------------------------------------------------------------
!
!***  *REAL FUNCTION* *C_QL(XK0,XK1,TH0,TH1)
!
!-----------------------------------------------------------------------
      REAL FUNCTION C_QL(XK0,XK1,TH0,TH1)
!
!***  *A*  DETERMINES THE QUASI-LINEAR TERM.
!
!     PETER JANSSEN
!
!     PURPOSE.
!     --------
!
!              DETERMINE CONTRIBUTION BY QUASI-LINEAR TERMS
!
!     INTERFACE.
!     ----------
!              *C_QL(XK0,XK1)*
!                        *XK0*  - WAVE NUMBER
!                        *XK1*  - WAVE NUMBER
!     METHOD.
!     -------
!              NONE
!
!     EXTERNALS.
!     ----------
!              NONE.
!
!-----------------------------------------------------------------------
!
      USE YOWPCONS, ONLY : G, PI

      IMPLICIT NONE
      REAL XK0,XK1,TH0,TH1,OM1,F1,OMEG,B2,B3
!
!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------
!
      OM1 = OMEG(XK1)
      F1  = SQRT(OM1/(2.*G))

      C_QL = 2./F1**2*(B2(XK0,XK1,XK1,XK0,TH0,TH1,TH1,TH0)+
     V                 B3(XK0,XK0,XK1,XK1,TH0-PI,TH0,TH1,TH1))
           
      RETURN
      END
!
!-----------------------------------------------------------------------
!
!***  *REAL FUNCTION* *U(XI,XJ,XK,XL,THI,THJ,THK,THL)
!
!-----------------------------------------------------------------------
      REAL FUNCTION U(XI,XJ,XK,XL,THI,THJ,THK,THL)
!
!***  *U*  DETERMINES THE THIRD-ORDER TRANSFER COEFFICIENT FOR FOUR
!              WAVE INTERACTIONS OF GRAVITY WAVES.
!
!     PETER JANSSEN
!
!     PURPOSE.
!     --------
!
!              GIVES NONLINEAR TRANSFER COEFFICIENT FOR FOUR
!              WAVE INTERACTIONS OF GRAVITY-CAPILLARY WAVES IN THE
!              IDEAL CASE OF NO CURRENT. (CF.ZAKHAROV,AND CRAWFORD ET AL)
!
!     INTERFACE.
!     ----------
!              *U(XI,XJ,XK,XL)*
!                      *XI*  - WAVE NUMBER
!                      *XJ*  - WAVE NUMBER
!                      *XK*  - WAVE NUMBER
!                      *XL*  - WAVE NUMBER
!     METHOD.
!     -------
!              NONE
!
!     EXTERNALS.
!     ----------
!              NONE.
!
!-----------------------------------------------------------------------
!
      USE YOWPCONS, ONLY : G
      IMPLICIT NONE
      REAL XI,XJ,XK,XL,THI,THJ,THK,THL,OI,OJ,OK,OL,XIK,XJK,XIL,XJL,
     V     OIK,OJK,OIL,OJL,QI,QJ,QIK,QJK,QIL,QJL,SQIJKL,ZCONST,VABS,OMEG
!
!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------
!
      ZCONST=1./(16.)

      OI=OMEG(XI)
      OJ=OMEG(XJ)
      OK=OMEG(XK)
      OL=OMEG(XL)

      XIK = VABS(XI,XK,THI,THK)
      XJK = VABS(XJ,XK,THJ,THK)
      XIL = VABS(XI,XL,THI,THL)
      XJL = VABS(XJ,XL,THJ,THL)
      OIK=OMEG(XIK)
      OJK=OMEG(XJK)
      OIL=OMEG(XIL)
      OJL=OMEG(XJL)

      QI=OI**2/G
      QJ=OJ**2/G
      QIK=OIK**2/G
      QJK=OJK**2/G
      QIL=OIL**2/G
      QJL=OJL**2/G
      SQIJKL=SQRT(OK*OL/(OI*OJ))
      U = ZCONST*SQIJKL*( 2.*(XI**2*QJ+XJ**2*QI)-QI*QJ*(
     V                  QIK+QJK+QIL+QJL) )
      RETURN
      END
!      
!-----------------------------------------------------------------------
!
!***  *REAL FUNCTION* *W2(XI,XJ,XK,XL,THI,THJ,THK,THL)
!
!-----------------------------------------------------------------------
      REAL FUNCTION W2(XI,XJ,XK,XL,THI,THJ,THK,THL)
!
!***  *W2*  DETERMINES THE CONTRIBUTION OF THE DIRECT FOUR-WAVE
!              INTERACTIONS OF GRAVITY WAVES OF THE TYPE
!              A_2^*A_3A_4.
!
!     PETER JANSSEN
!
!     PURPOSE.
!     --------
!
!              GIVES NONLINEAR TRANSFER COEFFICIENT FOR FOUR
!              WAVE INTERACTIONS OF GRAVITY WAVES IN THE
!              IDEAL CASE OF NO CURRENT. (CF.ZAKHAROV,AND CRAWFORD ET AL)
!
!     INTERFACE.
!     ----------
!              *W(XI,XJ,XK,XL)*
!                      *XI*  - WAVE NUMBER
!                      *XJ*  - WAVE NUMBER
!                      *XK*  - WAVE NUMBER
!                      *XL*  - WAVE NUMBER
!     METHOD.
!     -------
!              NONE
!
!     EXTERNALS.
!     ----------
!              NONE.
!
!-----------------------------------------------------------------------
!
      USE YOWPCONS, ONLY : PI
      IMPLICIT NONE
      REAL XI,XJ,XK,XL,THI,THJ,THK,THL,U
!
!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------
!
      W2= U(XI,XJ,XK,XL,THI-PI,THJ-PI,THK,THL)+
     V    U(XK,XL,XI,XJ,THK,THL,THI-PI,THJ-PI)-
     V    U(XK,XJ,XI,XL,THK,THJ-PI,THI-PI,THL)-
     V    U(XI,XK,XJ,XL,THI-PI,THK,THJ-PI,THL)-
     V    U(XI,XL,XK,XJ,THI-PI,THL,THK,THJ-PI)-
     V    U(XL,XJ,XK,XI,THL,THJ-PI,THK,THI-PI) 
      RETURN
      END
!      
!-----------------------------------------------------------------------
!
!***  *REAL FUNCTION* *V2(XI,XJ,XK,XL,THI,THJ,THK,THL)
!
!-----------------------------------------------------------------------
      REAL FUNCTION V2(XI,XJ,XK,XL,THI,THJ,THK,THL)
!
!***  *V2*  DETERMINES THE CONTRIBUTION OF THE VIRTUAL 
!           FOUR-WAVE INTERACTIONS OF GRAVITY WAVES.
!
!     PETER JANSSEN
!
!     PURPOSE.
!     --------
!
!              GIVES NONLINEAR TRANSFER COEFFICIENT FOR FOUR
!              WAVE INTERACTIONS OF GRAVITY WAVES IN THE
!              IDEAL CASE OF NO CURRENT. (CF.ZAKHAROV,AND 
!              CRAWFORD ET AL)
!
!     INTERFACE.
!     ----------
!              *V2(XI,XJ,XK,XL)*
!                      *XI*  - WAVE NUMBER
!                      *XJ*  - WAVE NUMBER
!                      *XK*  - WAVE NUMBER
!                      *XL*  - WAVE NUMBER
!     METHOD.
!     -------
!              NONE
!
!
!     EXTERNALS.
!     ----------
!              NONE.
!
!-----------------------------------------------------------------------
!
      USE YOWPCONS, ONLY : PI

      IMPLICIT NONE
      REAL DEL1,XI,XJ,XK,XL,THI,THJ,THK,THL,OI,OJ,OK,OL,RI,RJ,RK,RL,
     V     RIJ,RIK,RLI,RJL,RJK,RKL,THIJ,THIK,THLI,THJL,THJK,THKL,OIJ,
     V     OIK,OJL,OJK,OLI,OKL,XNIK,XNJL,XNJK,XNIL,YNIL,YNJK,YNJL,YNIK,
     V     ZNIJ,ZNKL,ZPIJ,ZPKL,THLJ,THIL,THKJ,THKI,THJI,THLK,VABS,VDIR,
     V     VMIN,VPLUS,OMEG
!
!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------
!
      DEL1=10.**(-5)

      RI=XI+DEL1
      RJ=XJ+DEL1/2.
      RK=XK+DEL1/3.
      RL=XL+DEL1*(1.+1./2.-1./3.)

      OI=OMEG(RI)
      OJ=OMEG(RJ)
      OK=OMEG(RK)
      OL=OMEG(RL)

      RIJ  = VABS(RI,RJ,THI,THJ)
      THIJ = VDIR(RI,RJ,THI,THJ)

      RIK  = VABS(RI,RK,THI,THK-PI)
      THIK = VDIR(RI,RK,THI,THK-PI)

      RLI  = VABS(RL,RI,THL,THI-PI)
      THLI = VDIR(XL,XI,THL,THI-PI)

      RJL  = VABS(RJ,RL,THJ,THL-PI)
      THJL = VDIR(RJ,RL,THJ,THL-PI)

      RJK  = VABS(RJ,RK,THJ,THK-PI)
      THJK = VDIR(RJ,RK,THJ,THK-PI)

      RKL  = VABS(RK,RL,THK,THL)
      THKL = VDIR(RK,RL,THK,THL)
   
      OIJ=OMEG(RIJ)
      OIK=OMEG(RIK)
      OJL=OMEG(RJL)
      OJK=OMEG(RJK)
      OLI=OMEG(RLI)
      OKL=OMEG(RKL)
      
      XNIK = OK+OIK-OI
      XNJL = OJ+OJL-OL
      XNJK = OK+OJK-OJ
      XNIL = OI+OLI-OL
      
      YNIL = OL+OLI-OI
      YNJK = OJ+OJK-OK
      YNJL = OL+OJL-OJ
      YNIK = OI+OIK-OK
      
      ZNIJ = OIJ-OI-OJ
      ZNKL = OKL-OK-OL
      ZPIJ = OIJ+OI+OJ
      ZPKL = OKL+OK+OL

      THLJ = THJL-PI
      THIL = THLI-PI
      THKJ = THJK-PI
      THKI = THIK-PI
      THJI = THIJ-PI
      THLK = THKL-PI 
      
      V2= VMIN(RI,RK,RIK,THI,THK,THIK)*VMIN(RL,RJ,RJL,THL,THJ,THLJ)*
     V    (1./XNIK+1./XNJL)
     V   +VMIN(RJ,RK,RJK,THJ,THK,THJK)*VMIN(RL,RI,RLI,THL,THI,THLI)*
     V    (1./XNJK+1./XNIL)
     V   +VMIN(RI,RL,RLI,THI,THL,THIL)*VMIN(RK,RJ,RJK,THK,THJ,THKJ)*
     V    (1./YNIL+1./YNJK)
     V   +VMIN(RJ,RL,RJL,THJ,THL,THJL)*VMIN(RK,RI,RIK,THK,THI,THKI)*
     V    (1./YNJL+1./YNIK)
     V   +VMIN(RIJ,RI,RJ,THIJ,THI,THJ)*VMIN(RKL,RK,RL,THKL,THK,THL)*
     V    (1./ZNIJ+1./ZNKL)
     V   +VPLUS(RIJ,RI,RJ,THJI,THI,THJ)*VPLUS(RKL,RK,RL,THLK,THK,THL)*
     V    (1./ZPIJ+1./ZPKL)

      V2 = -V2

      RETURN
      END
!      
!-----------------------------------------------------------------------
!
!***  *REAL FUNCTION* *W1(XI,XJ,XK,XL,THI,THJ,THK,THL)
!
!-----------------------------------------------------------------------
      REAL FUNCTION W1(XI,XJ,XK,XL,THI,THJ,THK,THL)
!
!***  *W1*  DETERMINES THE NONLINEAR TRANSFER COEFFICIENT FOR FOUR
!              WAVE INTERACTIONS OF GRAVITY WAVES OF THE TYPE
!              A_2A_3A_4.
!
!     PETER JANSSEN
!
!     PURPOSE.
!     --------
!
!              GIVES NONLINEAR TRANSFER COEFFICIENT FOR FOUR
!              WAVE INTERACTIONS OF GRAVITY-CAPILLARY WAVES IN THE
!              IDEAL CASE OF NO CURRENT. (CF.ZAKHAROV,AND CRAWFORD ET AL)
!
!     INTERFACE.
!     ----------
!              *W1(XI,XJ,XK,XL)*
!                      *XI*  - WAVE NUMBER
!                      *XJ*  - WAVE NUMBER
!                      *XK*  - WAVE NUMBER
!                      *XL*  - WAVE NUMBER
!     METHOD.
!     -------
!              NONE
!
!     EXTERNALS.
!     ----------
!              NONE.
!
!-----------------------------------------------------------------------
!
      USE YOWPCONS, ONLY : PI

      IMPLICIT NONE
      REAL XI,XJ,XK,XL,THI,THJ,THK,THL,U
!
!
!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------
!
      W1= -U(XI,XJ,XK,XL,THI-PI,THJ,THK,THL)-
     V     U(XI,XK,XJ,XL,THI-PI,THK,THJ,THL)-
     V     U(XI,XL,XJ,XK,THI-PI,THL,THJ,THK)+
     V     U(XJ,XK,XI,XL,THJ,THK,THI-PI,THL)+
     V     U(XJ,XL,XI,XK,THJ,THL,THI-PI,THK)+
     V     U(XK,XL,XI,XJ,THK,THL,THI-PI,THJ)

      W1=W1/3. 
               
      RETURN
      END
!
!***  *REAL FUNCTION* *W4(XI,XJ,XK,XL,THI,THJ,THK,THL)
!
!-----------------------------------------------------------------------
      REAL FUNCTION W4(XI,XJ,XK,XL,THI,THJ,THK,THL)
!
!***  *W4*  DETERMINES THE NONLINEAR TRANSFER COEFFICIENT FOR FOUR
!              WAVE INTERACTIONS OF GRAVITY WAVES of the type
!              A_^*A_3^*A_4^*.
!
!     PETER JANSSEN
!
!     PURPOSE.
!     --------
!
!              GIVES NONLINEAR TRANSFER COEFFICIENT FOR FOUR
!              WAVE INTERACTIONS OF GRAVITY-CAPILLARY WAVES IN THE
!              IDEAL CASE OF NO CURRENT. (CF.ZAKHAROV,AND CRAWFORD ET AL)
!
!     INTERFACE.
!     ----------
!              *W4(XI,XJ,XK,XL)*
!                      *XI*  - WAVE NUMBER
!                      *XJ*  - WAVE NUMBER
!                      *XK*  - WAVE NUMBER
!                      *XL*  - WAVE NUMBER
!     METHOD.
!     -------
!              NONE
!
!     EXTERNALS.
!     ----------
!              NONE.
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
      REAL XI,XJ,XK,XL,THI,THJ,THK,THL,U
!
!
!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------
!
 
      W4= U(XI,XJ,XK,XL,THI,THJ,THK,THL)+
     V    U(XI,XK,XJ,XL,THI,THK,THJ,THL)+
     V    U(XI,XL,XJ,XK,THI,THL,THJ,THK)+
     V    U(XJ,XK,XI,XL,THJ,THK,THI,THL)+
     V    U(XJ,XL,XI,XK,THJ,THL,THI,THK)+
     V    U(XK,XL,XI,XJ,THK,THL,THI,THJ)


      W4=W4/3. 
               
      RETURN
      END
!      
!-----------------------------------------------------------------------
!
!***  *REAL FUNCTION* *B3(XI,XJ,XK,XL,THI,THJ,THK,THL)
!
!-----------------------------------------------------------------------
      REAL FUNCTION B3(XI,XJ,XK,XL,THI,THJ,THK,THL)
!
!***  *B3*  WEIGHTS OF THE A_2^*A_3^*A_4 PART OF THE
!           CANONICAL TRANSFORMATION.
!
!     PETER JANSSEN
!
!     PURPOSE.
!     --------
!
!              GIVES NONLINEAR TRANSFER COEFFICIENT FOR FOUR
!              WAVE INTERACTIONS OF GRAVITY-CAPILLARY WAVES IN THE
!              IDEAL CASE OF NO CURRENT. (CF.ZAKHAROV,AND CRAWFORD ET AL)
!
!     INTERFACE.
!     ----------
!              *B3(XI,XJ,XK,XL)*
!                      *XI*  - WAVE NUMBER
!                      *XJ*  - WAVE NUMBER
!                      *XK*  - WAVE NUMBER
!                      *XL*  - WAVE NUMBER
!     METHOD.
!     -------
!              NONE
!
!
!     EXTERNALS.
!     ----------
!              NONE.
!
!-----------------------------------------------------------------------
!
      USE YOWPCONS, ONLY : PI

      IMPLICIT NONE
      REAL DEL1,XI,XJ,XK,XL,THI,THJ,THK,THL,OI,OJ,OK,OL,RI,RJ,RK,RL,
     V     RIJ,RJI,RIK,RKI,RLJ,RJL,RJK,RKJ,RLI,RIL,RLK,RKL,THIJ,THJI,
     V     THIK,THKI,THLJ,THJL,THJK,THKJ,THLI,THIL,THLK,THKL,ZIJKL,
     V     VABS,VDIR,VMIN,VPLUS,A1,A3,W1,OMEG
!
!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------
!
      DEL1=10.**(-5)

      RI=XI
      RJ=XJ
      RK=XK
      RL=XL

      OI=OMEG(RI)+DEL1
      OJ=OMEG(RJ)+DEL1
      OK=OMEG(RK)+DEL1
      OL=OMEG(RL)+DEL1

      RIJ  = VABS(RI,RJ,THI,THJ)
      THIJ = VDIR(RI,RJ,THI,THJ)

      RJI  = VABS(RJ,RI,THJ,THI)
      THJI = VDIR(RJ,RI,THJ,THI)

      RIK  = VABS(RI,RK,THI,THK)
      THIK = VDIR(RI,RK,THI,THK)

      RKI  = VABS(RK,RI,THK,THI)
      THKI = VDIR(RK,RI,THK,THI)

      RLJ  = VABS(RL,RJ,THL,THJ-PI)
      THLJ = VDIR(RL,RJ,THL,THJ-PI)

      RJL  = VABS(RJ,RL,THJ,THL-PI)
      THJL = VDIR(RJ,RL,THJ,THL-PI)

      RJK  = VABS(RJ,RK,THJ,THK)
      THJK = VDIR(RJ,RK,THJ,THK)

      RKJ  = VABS(RK,RJ,THK,THJ)
      THKJ = VDIR(RK,RJ,THK,THJ)

      RLI  = VABS(RL,RI,THL,THI-PI)
      THLI = VDIR(RL,RI,THL,THI-PI)

      RIL  = VABS(RI,RL,THI,THL-PI)
      THIL = VDIR(RI,RL,THI,THL-PI)

      RLK  = VABS(RL,RK,THL,THK-PI)
      THLK = VDIR(RL,RK,THL,THK-PI)

      RKL  = VABS(RK,RL,THK,THL-PI)
      THKL = VDIR(RK,RL,THK,THL-PI)

      ZIJKL = OI+OJ+OK-OL

      B3= -1./ZIJKL*(2.*( 
     V    VMIN(RL,RI,RLI,THL,THI,THLI)*A1(RJK,RJ,RK,THJK,THJ,THK)
     V   -VMIN(RIJ,RI,RJ,THIJ,THI,THJ)*A1(RL,RK,RLK,THL,THK,THLK)
     V   -VMIN(RIK,RI,RK,THIK,THI,THK)*A1(RL,RJ,RLJ,THL,THJ,THLJ)
     V   -VPLUS(RJ,RI,RJI,THJ,THI,THJI-PI)*A1(RK,RL,RKL,THK,THL,THKL)
     V   -VPLUS(RK,RI,RKI,THK,THI,THKI-PI)*A1(RJ,RL,RJL,THJ,THL,THJL)
     V   +VMIN(RI,RL,RIL,THI,THL,THIL)*A3(RJ,RK,RJK,THJ,THK,THJK-PI))
     V   +3.*W1(RL,RK,RJ,RI,THL,THK,THJ,THI) )

      RETURN
      END
!      
!-----------------------------------------------------------------------
!
!***  *REAL FUNCTION* *B4(XI,XJ,XK,XL,THI,THJ,THK,THL)
!
!-----------------------------------------------------------------------
      REAL FUNCTION B4(XI,XJ,XK,XL,THI,THJ,THK,THL)
!
!***  *B4*  WEIGHTS OF THE A_2^*A_3^*A_4^* PART OF THE CANONICAL
!           TRANSFORMATION.
!
!     PETER JANSSEN
!
!     PURPOSE.
!     --------
!
!              GIVES NONLINEAR TRANSFER COEFFICIENT FOR FOUR
!              WAVE INTERACTIONS OF GRAVITY-CAPILLARY WAVES IN THE
!              IDEAL CASE OF NO CURRENT. (CF.ZAKHAROV,AND CRAWFORD ET AL)
!
!     INTERFACE.
!     ----------
!              *B4(XI,XJ,XK,XL)*
!                      *XI*  - WAVE NUMBER
!                      *XJ*  - WAVE NUMBER
!                      *XK*  - WAVE NUMBER
!                      *XL*  - WAVE NUMBER
!     METHOD.
!     -------
!              NONE
!
!
!     EXTERNALS.
!     ----------
!              NONE.
!
!-----------------------------------------------------------------------
!
      USE YOWPCONS, ONLY : PI
      IMPLICIT NONE
      REAL DEL1,XI,XJ,XK,XL,THI,THJ,THK,THL,OI,OJ,OK,OL,RI,RJ,RK,RL,
     V     RIJ,RIK,RIL,RJL,RJK,RKL,THIJ,THIK,THIL,THJL,THJK,THLK,THKL,
     V     ZIJKL,VABS,VDIR,VMIN,VPLUS,A1,A3,W4,OMEG
!
!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------
!


      RI=XI
      RJ=XJ
      RK=XK
      RL=XL

      OI=OMEG(RI)
      OJ=OMEG(RJ)
      OK=OMEG(RK)
      OL=OMEG(RL)


      RIJ  = VABS(RI,RJ,THI,THJ)
      THIJ = VDIR(RI,RJ,THI,THJ)

      RIK  = VABS(RI,RK,THI,THK)
      THIK = VDIR(RI,RK,THI,THK)

      RIL  = VABS(RI,RL,THI,THL)
      THIL = VDIR(RI,RL,THI,THL)

      RJL  = VABS(RJ,RL,THJ,THL)
      THJL = VDIR(RJ,RL,THJ,THL)

      RJK  = VABS(RJ,RK,THJ,THK)
      THJK = VDIR(RJ,RK,THJ,THK)

      RKL  = VABS(RK,RL,THK,THL)
      THKL = VDIR(RK,RL,THK,THL)      

      
      ZIJKL = OI+OJ+OK+OL

      B4= -1./ZIJKL*(2./3.*( 
     V     VPLUS(RIJ,RI,RJ,THIJ-PI,THI,THJ)*A1(RKL,RK,RL,THKL,THK,THL)
     V    +VPLUS(RIK,RI,RK,THIK-PI,THI,THK)*A1(RJL,RJ,RL,THJL,THJ,THL)
     V    +VPLUS(RIL,RI,RL,THIL-PI,THI,THL)*A1(RJK,RJ,RK,THJK,THJ,THK)
     V    +VMIN(RIK,RI,RK,THIK,THI,THK)*A3(RJL,RJ,RL,THJL-PI,THJ,THL)
     V    +VMIN(RIL,RI,RL,THIL,THI,THL)*A3(RJK,RJ,RK,THJK-PI,THJ,THK)
     V    +VMIN(RIJ,RI,RJ,THIJ,THI,THJ)*A3(RKL,RK,RL,THKL-PI,THK,THL) )
     V    +W4(RI,RJ,RK,RL,THI,THJ,THK,THL) )

      RETURN
      END
!      
!-----------------------------------------------------------------------
!
!***  *REAL FUNCTION* *B1(XI,XJ,XK,XL,THI,THJ,THK,THL)
!
!-----------------------------------------------------------------------
      REAL FUNCTION B1(XI,XJ,XK,XL,THI,THJ,THK,THL)
!
!***  *B1*  WEIGHTS OF THE A_2A_3A_4 PART OF THE CANONICAL
!           TRANSFORMATION.
!
!     PETER JANSSEN
!
!     PURPOSE.
!     --------
!
!              GIVES NONLINEAR TRANSFER COEFFICIENT FOR FOUR
!              WAVE INTERACTIONS OF GRAVITY-CAPILLARY WAVES IN THE
!              IDEAL CASE OF NO CURRENT. (CF.ZAKHAROV,AND CRAWFORD ET AL)
!
!     INTERFACE.
!     ----------
!              *B1(XI,XJ,XK,XL)*
!                      *XI*  - WAVE NUMBER
!                      *XJ*  - WAVE NUMBER
!                      *XK*  - WAVE NUMBER
!                      *XL*  - WAVE NUMBER
!     METHOD.
!     -------
!              NONE
!
!
!     EXTERNALS.
!     ----------
!              NONE.
!
!-----------------------------------------------------------------------
!
      USE YOWPCONS, ONLY : PI
      IMPLICIT NONE
      REAL DEL1,XI,XJ,XK,XL,THI,THJ,THK,THL,OI,OJ,OK,OL,RI,RJ,RK,RL,
     V     RIJ,RJI,RIK,RKI,RJL,RJK,RLI,RIL,RKL,THIJ,THJI,
     V     THIK,THKI,THJL,THJK,THLI,THIL,THKL,ZIJKL,
     V     VABS,VDIR,VMIN,A1,A3,W1,OMEG
!
!
!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------
!

      RI=XI
      RJ=XJ
      RK=XK
      RL=XL

      OI=OMEG(RI)
      OJ=OMEG(RJ)
      OK=OMEG(RK)
      OL=OMEG(RL)

      RIJ  = VABS(RI,RJ,THI,THJ-PI)
      THIJ = VDIR(RI,RJ,THI,THJ-PI)

      RJI  = VABS(RJ,RI,THJ,THI-PI)
      THJI = VDIR(RJ,RI,THJ,THI-PI)

      RIK  = VABS(RI,RK,THI,THK-PI)
      THIK = VDIR(RI,RK,THI,THK-PI)

      RKI  = VABS(RK,RI,THK,THI-PI)
      THKI = VDIR(RK,RI,THK,THI-PI)

      RIL  = VABS(RI,RL,THI,THL-PI)
      THIL = VDIR(RI,RL,THI,THL-PI)

      RLI  = VABS(RL,RI,THL,THI-PI)
      THLI = VDIR(RL,RI,THL,THI-PI)

      RJL  = VABS(RJ,RL,THJ,THL)
      THJL = VDIR(RJ,RL,THJ,THL)

      RJK  = VABS(RJ,RK,THJ,THK)
      THJK = VDIR(RJ,RK,THJ,THK)

      RKL  = VABS(RK,RL,THK,THL)
      THKL = VDIR(RK,RL,THK,THL)      
      
      ZIJKL = OI-OJ-OK-OL
      
      B1= -1./ZIJKL*(2./3.*( 
     V     VMIN(RI,RJ,RIJ,THI,THJ,THIJ)*A1(RKL,RK,RL,THKL,THK,THL)
     V    +VMIN(RI,RK,RIK,THI,THK,THIK)*A1(RJL,RJ,RL,THJL,THJ,THL)
     V    +VMIN(RI,RL,RIL,THI,THL,THIL)*A1(RJK,RJ,RK,THJK,THJ,THK)
     V    +VMIN(RK,RI,RKI,THK,THI,THKI)*A3(RJL,RJ,RL,THJL-PI,THJ,THL)
     V    +VMIN(RL,RI,RLI,THL,THI,THLI)*A3(RJK,RJ,RK,THJK-PI,THJ,THK)
     V    +VMIN(RJ,RI,RJI,THJ,THI,THJI)*A3(RKL,RK,RL,THKL-PI,THK,THL) 
     V    ) +W1(RI,RJ,RK,RL,THI,THJ,THK,THL) )    
      RETURN
      END

!      
!-----------------------------------------------------------------------
!
!***  *REAL FUNCTION* *B2(XI,XJ,XK,XL,THI,THJ,THK,THL)
!
!-----------------------------------------------------------------------
      REAL FUNCTION B2(XI,XJ,XK,XL,THI,THJ,THK,THL)
!
!
!***  *B2*  WEIGHTS OF THE A_2^*A_3A_4 PART OF THE CANONICAL
!           TRANSFORMATION.
!
!     PETER JANSSEN
!
!     PURPOSE.
!     --------
!
!              GIVES NONLINEAR TRANSFER COEFFICIENT FOR FOUR
!              WAVE INTERACTIONS OF GRAVITY-CAPILLARY WAVES IN THE
!              IDEAL CASE OF NO CURRENT. (CF.ZAKHAROV,AND CRAWFORD ET AL)
!
!     INTERFACE.
!     ----------
!              *B2(XI,XJ,XK,XL)*
!                      *XI*  - WAVE NUMBER
!                      *XJ*  - WAVE NUMBER
!                      *XK*  - WAVE NUMBER
!                      *XL*  - WAVE NUMBER
!     METHOD.
!     -------
!              NONE
!
!
!     EXTERNALS.
!     ----------
!              NONE.
!
!-----------------------------------------------------------------------
!
      USE YOWPCONS, ONLY : PI
      IMPLICIT NONE
      REAL DEL1,XI,XJ,XK,XL,THI,THJ,THK,THL,OI,OJ,OK,OL,RI,RJ,RK,RL,
     V     RIJ,RIK,RKI,RJL,RLJ,RJK,RKJ,RLI,RIL,RKL,THIJ,
     V     THIK,THKI,THJL,THLJ,THJK,THKJ,THLI,THIL,THKL,ZIJKL,
     V     VABS,VDIR,A1,A3,OMEG
!
!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------
!

      RI=XI
      RJ=XJ
      RK=XK
      RL=XL

      RIJ  = VABS(RI,RJ,THI,THJ)
      THIJ = VDIR(RI,RJ,THI,THJ)

      RIK  = VABS(RI,RK,THI,THK-PI)
      THIK = VDIR(RI,RK,THI,THK-PI)

      RKI  = VABS(RK,RI,THK,THI-PI)
      THKI = VDIR(RK,RI,THK,THI-PI)

      RIL  = VABS(RI,RL,THI,THL-PI)
      THIL = VDIR(RI,RL,THI,THL-PI)

      RLI  = VABS(RL,RI,THL,THI-PI)
      THLI = VDIR(RL,RI,THL,THI-PI)

      RJL  = VABS(RJ,RL,THJ,THL-PI)
      THJL = VDIR(RJ,RL,THJ,THL-PI)

      RLJ  = VABS(RL,RJ,THL,THJ-PI)
      THLJ = VDIR(RL,RJ,THL,THJ-PI)

      RJK  = VABS(RJ,RK,THJ,THK-PI)
      THJK = VDIR(RJ,RK,THJ,THK-PI)

      RKJ  = VABS(RK,RJ,THK,THJ-PI)
      THKJ = VDIR(RK,RJ,THK,THJ-PI)

      RKL  = VABS(RK,RL,THK,THL)
      THKL = VDIR(RK,RL,THK,THL)   

      B2=  A3(RI,RJ,RIJ,THI,THJ,THIJ-PI)*A3(RK,RL,RKL,THK,THL,THKL-PI)
     V    +A1(RJ,RK,RJK,THJ,THK,THJK)*A1(RL,RI,RLI,THL,THI,THLI)
     V    +A1(RJ,RL,RJL,THJ,THL,THJL)*A1(RK,RI,RKI,THK,THI,THKI)
     V    -A1(RIJ,RI,RJ,THIJ,THI,THJ)*A1(RKL,RK,RL,THKL,THK,THL)
     V    -A1(RI,RK,RIK,THI,THK,THIK)*A1(RL,RJ,RLJ,THL,THJ,THLJ)
     V    -A1(RI,RL,RIL,THI,THL,THIL)*A1(RK,RJ,RKJ,THK,THJ,THKJ)


      RETURN
      END
!
!-----------------------------------------------------------------------
!
!***  *REAL FUNCTION* *A1(XI,XJ,XK,THI,THJ,THK)
!
!-----------------------------------------------------------------------
      REAL FUNCTION A1(XI,XJ,XK,THI,THJ,THK)
!
!***  *A1*  AUXILIARY SECOND-ORDER COEFFICIENT.
!
!     PETER JANSSEN
!
!     PURPOSE.
!     --------
!
!              GIVES NONLINEAR TRANSFER COEFFICIENT FOR THREE
!              WAVE INTERACTIONS OF GRAVITY-CAPILLARY WAVES IN THE
!              IDEAL CASE OF NO CURRENT. (CF.ZAKHAROV)
!
!     INTERFACE.
!     ----------
!              *VMIN(XI,XJ,XK)*
!                      *XI*  - WAVE NUMBER
!                      *XJ*  - WAVE NUMBER
!                      *XK*  - WAVE NUMBER
!     METHOD.
!     -------
!              NONE
!
!     EXTERNALS.
!     ----------
!              NONE.
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE

      REAL DEL1,XI,XJ,XK,THI,THJ,THK,OI,OJ,OK,
     V     VMIN,OMEG
!
!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------
!
      DEL1 = 10.**(-8)

      OI=OMEG(XI)+DEL1
      OJ=OMEG(XJ)+DEL1
      OK=OMEG(XK)+DEL1

      A1 = -VMIN(XI,XJ,XK,THI,THJ,THK)/(OI-OJ-OK)

      RETURN
      END
!
!-----------------------------------------------------------------------
!
!***  *REAL FUNCTION* *A2(XI,XJ,XK,THI,THJ,THK)
!
!-----------------------------------------------------------------------
      REAL FUNCTION A2(XI,XJ,XK,THI,THJ,THK)
!
!***  *A2*  AUXILIARY SECOND-ORDER FUNCTION.
!
!     PETER JANSSEN
!
!     PURPOSE.
!     --------
!
!              GIVES NONLINEAR TRANSFER COEFFICIENT FOR THREE
!              WAVE INTERACTIONS OF GRAVITY-CAPILLARY WAVES IN THE
!              IDEAL CASE OF NO CURRENT. (CF.ZAKHAROV)
!
!     INTERFACE.
!     ----------
!              *VMIN(XI,XJ,XK)*
!                      *XI*  - WAVE NUMBER
!                      *XJ*  - WAVE NUMBER
!                      *XK*  - WAVE NUMBER
!     METHOD.
!     -------
!              NONE
!
!     EXTERNALS.
!     ----------
!              NONE.
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
      REAL DEL1,XI,XJ,XK,THI,THJ,THK,A1
!
!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------
!
      A2 = -2.*A1(XK,XJ,XI,THK,THJ,THI)
      RETURN
      END
!
!-----------------------------------------------------------------------
!
!***  *REAL FUNCTION* *A3(XI,XJ,XK,THI,THJ,THK)
!
!-----------------------------------------------------------------------
      REAL FUNCTION A3(XI,XJ,XK,THI,THJ,THK)
!
!***  *A3*  AUXILIARY SECOND-ORDER FUNCTION.
!
!     PETER JANSSEN
!
!     PURPOSE.
!     --------
!
!              GIVES NONLINEAR TRANSFER COEFFICIENT FOR THREE
!              WAVE INTERACTIONS OF GRAVITY-CAPILLARY WAVES IN THE
!              IDEAL CASE OF NO CURRENT. (CF.ZAKHAROV)
!
!     INTERFACE.
!     ----------
!              *VMIN(XI,XJ,XK)*
!                      *XI*  - WAVE NUMBER
!                      *XJ*  - WAVE NUMBER
!                      *XK*  - WAVE NUMBER
!     METHOD.
!     -------
!              NONE
!
!     EXTERNALS.
!     ----------
!              NONE.
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE

      REAL DEL1,OI,OJ,OK,XI,XJ,XK,THI,THJ,THK,OMEG,VPLUS
!
!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------
!    
      DEL1 = 10.**(-8)

      OI=OMEG(XI)+DEL1
      OJ=OMEG(XJ)+DEL1
      OK=OMEG(XK)+DEL1

      A3 = -VPLUS(XI,XJ,XK,THI,THJ,THK)/(OI+OJ+OK)
      RETURN
      END

!
!-----------------------------------------------------------------------
!
!
!***  *REAL FUNCTION* *OMEG(X)*
!
!-----------------------------------------------------------------------
!
      REAL FUNCTION OMEG(X)
!
!***  *OMEG*   DETERMINES THE DISPERSION RELATION FOR GRAVITY
!              WAVES.
!
!     PETER JANSSEN
!
!     PURPOSE.
!     --------
!
!              GIVES DISPERSION RELATION FOR GRAVITY-
!              WAVES IN THE IDEAL CASE OF NO CURRENT.
!
!     INTERFACE.
!     ----------
!              *OMEG(X)*
!                      *X*  - WAVE NUMBER
!
!     METHOD.
!     -------
!              NONE
!
!     EXTERNALS.
!     ----------
!              NONE.
!
!-----------------------------------------------------------------------
!
      USE YOWPCONS, ONLY : G
      USE YOWCONST_2ND, ONLY: DPTH
      IMPLICIT NONE
      REAL D,XK,X,T

      D = DPTH
      XK = ABS(X)
      T = TANH(XK*D)
      OMEG=SQRT(G*XK*T)

      RETURN
      END
!
      REAL FUNCTION VABS(XI,XJ,THI,THJ)
!
!---------------------------------------------------------------------
!
      IMPLICIT NONE
      REAL XI,XJ,THI,THJ,ARG

      ARG = XI**2+XJ**2+2.*XI*XJ*COS(THI-THJ)

      IF (ARG.LE.0.) THEN
         VABS = 0.
      ELSE
         VABS = SQRT(ARG)
      ENDIF

      RETURN
      END
!
      REAL FUNCTION VDIR(XI,XJ,THI,THJ)
!
!---------------------------------------------------------------------
!
      IMPLICIT NONE
      REAL XI,XJ,THI,THJ,EPS,Y,X

      EPS = 0.

      Y = XJ*SIN(THJ-THI)
      X = XI+XJ*COS(THJ-THI)+EPS
      VDIR = ATAN2(Y,X)+THI
      IF (X.EQ.0.) VDIR = 0.

      RETURN
      END
