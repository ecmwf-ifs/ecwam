!
!-----------------------------------------------------------------------
!
!***  *REAL(KIND=JWRB) FUNCTION* *VPLUS(XI,XJ,XK,THI,THJ,THK)
!
!-----------------------------------------------------------------------
      REAL(KIND=JWRB) FUNCTION VPLUS(XI,XJ,XK,THI,THJ,THK)
!
!***  *VPLUS*  DETERMINES THE SECOND-ORDER TRANSFER COEFFICIENT 
!              FOR THREE WAVE INTERACTIONS OF GRAVITY WAVES.
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
!              *VPLUS(XI,XJ,XK)*
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
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS, ONLY : G

!-----------------------------------------------------------------------
      IMPLICIT NONE

      REAL(KIND=JWRB) :: DEL1, RI, RJ, RK, XI, XJ, XK,                  &
     &                   THI, THJ, THK, OI, OJ, OK, QI, QJ, QK,         &
     &                   RIJ, RIK, RJK, SQIJK, SQIKJ, SQJKI,            &
     &                   ZCONST, OMEG
!
!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------
!
      DEL1 = 10.0_JWRB**(-12)
      ZCONST=1.0_JWRB/(4.0_JWRB*SQRT(2.0_JWRB))

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

      VPLUS=ZCONST*( (RIJ+QI*QJ)*SQIJK + (RIK+QI*QK)*SQIKJ + (RJK+QJ*QK)*SQJKI )

      END FUNCTION VPLUS

