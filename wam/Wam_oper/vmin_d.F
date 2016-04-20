!-----------------------------------------------------------------------

!***  *REAL FUNCTION* *VMIN_D(XI,XJ,XK,THI,THJ,THK,OI,OJ,OK)

!-----------------------------------------------------------------------
      REAL FUNCTION VMIN_D(XI,XJ,XK,XIJ,XIK,XJK,XOI,XOJ,XOK)

!***  *VMIN_D*  DETERMINES THE NONLINEAR TRANSFER COEFFICIENT FOR THREE
!               WAVE INTERACTIONS OF DEEP WATER WAVES.

!     PETER JANSSEN

!     PURPOSE.
!     --------

!              GIVES NONLINEAR TRANSFER COEFFICIENT FOR THREE
!              WAVE INTERACTIONS OF DEEP-WATER WAVES IN THE
!              IDEAL CASE OF NO CURRENT. (CF.ZAKHAROV)

!     INTERFACE.
!     ----------
!              *VMIN_D(XI,XJ,XK)*
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

      USE YOWPCONS , ONLY : ZCONST

!-----------------------------------------------------------------------

      IMPLICIT NONE

      REAL :: XI,XJ,XK,XIJ,XIK,XJK,XOI,XOJ,XOK

      REAL, PARAMETER :: DEL1=1.0E-10

      REAL :: RI,RJ,RK,OI,OJ,OK,SQIJK,SQIKJ,SQJKI

!***  1. DETERMINE NONLINEAR TRANSFER.
!     --------------------------------

      RI=ABS(XI)+DEL1
      RJ=ABS(XJ)+DEL1
      RK=ABS(XK)+DEL1
      OI=XOI+DEL1
      OJ=XOJ+DEL1
      OK=XOK+DEL1
      SQIJK=SQRT(OI*OJ*RK/(OK*RI*RJ))
      SQIKJ=SQRT(OI*OK*RJ/(OJ*RI*RK))
      SQJKI=SQRT(OJ*OK*RI/(OI*RJ*RK))
      VMIN_D=ZCONST*( (XIJ-RI*RJ)*SQIJK + (XIK-RI*RK)*SQIKJ
     &                + (XJK+RJ*RK)*SQJKI )

      RETURN
      END FUNCTION VMIN_D
