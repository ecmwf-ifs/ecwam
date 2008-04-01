      REAL FUNCTION TRANSF2(XK,D,G)
!
!***  DETERMINE NARROW BAND LIMIT BENJAMIN-FEIR INDEX FOR
!     THE FINITE DEPTH CASE.
!
!
!     VARIABLE       TYPE         PURPOSE
!     --------       ----         -------
!
!     XK             REAL         WAVE NUMBER
!     D              REAL         DEPTH
!     G              REAL         ACCELERATION OF GRAVITY
!
!
!     BF**2 = (2 S^2)/SIG_OM^2) . TRANSF2(XK,D) 
!
       
! ----------------------------------------------------------------------
 
      IMPLICIT NONE

      REAL :: XK,D,G,X,T_0,OM,C_0,C_S,V_G,D2OM,XNL_1,XNL_2,T_NL

!----------------------------------------------------------------------
!
!*    1. DETERMINE TRANSFER FUNCTION.
!     ------------------------------
!     
      X   = MIN(XK*D,40.)
      X   = MAX(X,0.5)
      T_0 = TANH(X)
      OM  = SQRT(G*XK*T_0)
      C_0 = OM/XK
      C_S = SQRT(G*D)
      V_G = 0.5*C_0*(1.+2.*X/SINH(2.*X))
      D2OM = V_G**2-C_S**2*(1.-T_0**2)*(1.-X*T_0)
      
      XNL_1 = (9.*T_0**4-10.*T_0**2+9.)/(8.*T_0**3)
      XNL_2 = ((2.*V_G-0.5*C_0)**2/(G*D-V_G**2)+1.)/X

      T_NL = XK**3*(XNL_1-XNL_2)
      TRANSF2 = (V_G/C_0)**2*G*T_NL/XK**4/D2OM
!
      RETURN
      END
