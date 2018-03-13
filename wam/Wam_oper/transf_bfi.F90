      REAL FUNCTION TRANSF_BFI(XK,D,XNU,SIG_TH)
 
!***  DETERMINE NARROW BAND LIMIT BENJAMIN-FEIR INDEX FOR
!     THE FINITE DEPTH CASE             
 
 
!     VARIABLE       TYPE         PURPOSE
!     --------       ----         -------
 
!     XK             REAL         WAVE NUMBER
!     D              REAL         DEPTH
 
!     XNU            REAL         RELATIVE SPECTRAL WIDTH
!     SIG_TH         REAL         RELATIVE WIDTH in DIRECTION
 
!     BF**2 = (2 S^2)/SIG_OM^2) . (2 V_G/C_0)^2 . T_0/K_0^3 
 
!----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWPCONS , ONLY : G     ,DKMAX
      USE YOWSHAL , ONLY : BATHYMAX
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
                     
!----------------------------------------------------------------------
      IMPLICIT NONE

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: X,XK,D,T_0,T_0_SQ,OM,C_0,C_S_SQ,V_G,V_G_SQ,D2OM
      REAL(KIND=JWRB) :: XNL_1,XNL_2,XNL_3,XNL_4
      REAL(KIND=JWRB) :: XNL,ALP,ZFAC,SIG_TH,XNU,T_NL
!----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('TRANSF_BFI',0,ZHOOK_HANDLE)
#endif

!*    1. DETERMINE TRANSFER FUNCTION.
!     ------------------------------
     
      IF(D.LT.BATHYMAX .AND. D.GT.0._JWRB) THEN
        X   = XK*D
        IF ( X .GT. DKMAX) THEN
          TRANSF_BFI = 1._JWRB 
        ELSE
          X   = MAX(X,0.5_JWRB)
          T_0 = TANH(X)
          T_0_SQ = T_0**2
          OM  = SQRT(G*XK*T_0)
          C_0 = OM/XK
          C_S_SQ = G*D
          V_G = 0.5_JWRB*C_0*(1._JWRB+2._JWRB*X/SINH(2._JWRB*X))
          V_G_SQ = V_G**2
          D2OM = (T_0-X*(1._JWRB-T_0_SQ))**2+4._JWRB*X**2*T_0_SQ*(1._JWRB-T_0_SQ)
      
          XNL_1 = (9._JWRB*T_0_SQ**2-10._JWRB*T_0_SQ+9._JWRB)/(8._JWRB*T_0_SQ*T_0)
          XNL_2 = ((2._JWRB*V_G-0.5_JWRB*C_0)**2/(G*D-V_G_SQ)+1._JWRB)/X
          XNL_4 = 1./(4._JWRB*T_0)*(2._JWRB*C_0+V_G*(1._JWRB-T_0_SQ))**2/(C_S_SQ-V_G_SQ)
          ALP = (1._JWRB-V_G_SQ/C_S_SQ)*C_0**2/V_G_SQ
          ZFAC = SIG_TH**2/(SIG_TH**2+ALP*XNU**2)
          XNL_3 = ZFAC*XNL_4

!!          T_NL = XK**3*(XNL_1-XNL_2+XNL_3)
!!          TRANSF_BFI = 4._JWRB*(V_G/C_0)**2*T_NL/XK**3*T_0/D2OM
          T_NL = XNL_1-XNL_2+XNL_3
          TRANSF_BFI = 4._JWRB*(V_G/C_0)**2*T_NL*T_0/D2OM
          TRANSF_BFI = MIN(5._JWRB,TRANSF_BFI)
          TRANSF_BFI = MAX(-5._JWRB,TRANSF_BFI)
        ENDIF
      ELSE
        TRANSF_BFI = 1._JWRB
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('TRANSF_BFI',1,ZHOOK_HANDLE)
#endif
      END FUNCTION TRANSF_BFI
