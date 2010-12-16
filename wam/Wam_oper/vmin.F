!
!-----------------------------------------------------------------------
!
!***  *REAL FUNCTION* *VMIN(XI,XJ,XK,THI,THJ,THK)
!
!-----------------------------------------------------------------------
      REAL FUNCTION VMIN(XI,XJ,XK,THI,THJ,THK)
!
!***  *VMIN*  DETERMINES THE SECOND-ORDER TRANSFER COEFFICIENT FOR 
!             THREE WAVE INTERACTIONS OF GRAVITY WAVES.
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
!              *VMIN(XI,XJ,XK)*
!                      *XI*   - WAVE NUMBER
!                      *XJ*   - WAVE NUMBER
!                      *XK*   - WAVE NUMBER
!                      *THI*  - WAVE DIRECTION
!                      *THJ*  - WAVE DIRECTION
!                      *THK*  - WAVE DIRECTION
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
      REAL DEL1,RI,RJ,RK,XI,XJ,XK,THI,THJ,THK,OI,OJ,OK,QI,QJ,QK,
     V     RIJ,RIK,RJK,SQIJK,SQIKJ,SQJKI,ZCONST,OMEG
!
!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------
!
      DEL1 = 10.**(-12)
      ZCONST=1./(4*SQRT(2.))

      RI = XI
      RJ = XJ
      RK = XK

      OI=OMEG(RI)+DEL1
      OJ=OMEG(RJ)+DEL1
      OK=OMEG(RK)+DEL1

      QI=OI**2/G
      QJ=OJ**2/G
      QK=OK**2/G

      RIJ = RI*RJ*COS(THJ-THI)
      RIK = RI*RK*COS(THK-THI)
      RJK = RJ*RK*COS(THK-THJ)

      SQIJK=SQRT(G*OK/(OI*OJ))
      SQIKJ=SQRT(G*OJ/(OI*OK))
      SQJKI=SQRT(G*OI/(OJ*OK))

      VMIN=ZCONST*( (RIJ-QI*QJ)*SQIJK + (RIK-QI*QK)*SQIKJ
     V               + (RJK+QJ*QK)*SQJKI )
      RETURN
      END

