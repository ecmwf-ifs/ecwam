
!----------------------------------------------------------------

!***  *REAL FUNCTION* *VPLUS(XI,XJ,XK,THI,THJ,THK,OI,OJ,OK)

!-----------------------------------------------------------------------

      REAL FUNCTION VPLUS(XI,XJ,XK,XIJ,XIK,XJK,XOI,XOJ,XOK)

!***  *VPLUS*  DETERMINES THE NONLINEAR TRANSFER COEFFICIENT FOR THREE
!              WAVE INTERACTIONS OF GRAVITY-CAPILLARY WAVES.

!     PETER JANSSEN

!     PURPOSE.
!     --------

!              GIVES NONLINEAR TRANSFER COEFFICIENT FOR THREE
!              WAVE INTERACTIONS OF GRAVITY-CAPILLARY WAVES IN THE
!              IDEAL CASE OF NO CURRENT. (CF.ZAKHAROV)

!     INTERFACE.
!     ----------
!              *VPLUS(XI,XJ,XK)*
!                      *XI*  - WAVE NUMBER
!                      *XJ*  - WAVE NUMBER
!                      *XK*  - WAVE NUMBER
!     METHOD.
!     -------
!              NONE

!     EXTERNALS.
!     ----------
!              NONE.

!-----------------------------------------------------------------------


!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------

      DEL1=10.**(-12)
      PI=4.*ATAN(1.)
      ZCONST=1./(8.*PI*SQRT(2.))

      RI=ABS(XI)+DEL1
      RJ=ABS(XJ)+DEL1
      RK=ABS(XK)+DEL1
      OI=XOI+DEL1
      OJ=XOJ+DEL1
      OK=XOK+DEL1
      SQIJK=SQRT(OI*OJ*RK/(OK*RI*RJ))
      SQIKJ=SQRT(OI*OK*RJ/(OJ*RI*RK))
      SQJKI=SQRT(OJ*OK*RI/(OI*RJ*RK))
      VPLUS=ZCONST*( (XIJ+RI*RJ)*SQIJK + (XIK+RI*RK)*SQIKJ
     &               + (XJK+RJ*RK)*SQJKI )
      RETURN
      END FUNCTION VPLUS
