      REAL FUNCTION TRANSF2(XK,D)
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
!
!
!     BF**2 = (2 S^2)/SIG_OM^2) . TRANSF2(XK,D) 
!
       
! ----------------------------------------------------------------------

      USE YOWPCONS , ONLY : G     ,DKMAX

!----------------------------------------------------------------------

      IMPLICIT NONE

      REAL :: XK,D,X,T_0,OM,C_0,C_S,V_G,D2OM,XNL_1,XNL_2,T_NL,ZSEC,ZEPS

!----------------------------------------------------------------------
!
!*    1. DETERMINE TRANSFER FUNCTION.
!     ------------------------------
!
!     We implicitely imply use 999 as the threshold for deep water
      IF(D.LT.999.) THEN
      !ZEPS= 1000.* TINY(ZEPS)
      ZEPS= 100.* EPSILON(ZEPS)
      X   = MIN(XK*D,DKMAX)
      X   = MAX(X,0.5)
      T_0 = SIGN(MAX(ABS(TANH(X)),ZEPS),TANH(X))
      OM  = SQRT(G*XK*T_0)
      C_0 = OM/XK
      C_S = SQRT(G*D)
      V_G = 0.5*C_0*(1.+2.*X/SINH(2.*X))
      D2OM = V_G**2-C_S**2*(1.-T_0**2)*(1.-X*T_0)
      D2OM = SIGN(MAX(ABS(D2OM),ZEPS),D2OM)
      
      XNL_1 = (9.*T_0**4-10.*T_0**2+9.)/(8.*T_0**3)
      ZSEC  = G*D-V_G**2
      ZSEC  = SIGN(MAX(ABS(ZSEC),ZEPS),ZSEC)
      XNL_2 = ((2.*V_G-0.5*C_0)**2/ZSEC+1.)/X

!      T_NL = XK**3*(XNL_1-XNL_2)
!     TRANSF2 = (V_G/C_0)**2*G*T_NL/XK**4/D2OM
!     Alternative formulation
      T_NL = XNL_1-XNL_2
      TRANSF2 =  V_G*V_G*T_NL/T_0 / D2OM

      ELSE
        TRANSF2 = 1.0
      ENDIF

!!! debile
         write(*,*) 'debile transf2 = ',transf2
!
      RETURN
      END
