      FUNCTION TRANSF_SNL(XK,D,XNU,SIG_TH)
 
!***  DETERMINE NARROW BAND LIMIT NONLINEAR TRANSFER FUNCTION            
 
!     VARIABLE       TYPE         PURPOSE
!     --------       ----         -------
 
!     XK             REAL         WAVE NUMBER
!     D              REAL         DEPTH
!     XNU            REAL         RELATIVE SPECTRAL WIDTH
!     SIG_TH         REAL         RELATIVE WIDTH in DIRECTION
 
!----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWPCONS , ONLY : G     ,DKMAX
      USE YOWSHAL , ONLY : BATHYMAX, XKDMIN
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
                     
!----------------------------------------------------------------------
      IMPLICIT NONE

      REAL(KIND=JWRB) :: TRANSF_SNL

      REAL(KIND=JWRB), PARAMETER :: EPS=0.0001_JWRB
      REAL(KIND=JWRB), PARAMETER :: TRANSF_SNL_MIN=0.1_JWRB
      REAL(KIND=JWRB), PARAMETER :: TRANSF_SNL_MAX=10._JWRB

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: X,XK,D,T_0,T_0_SQ,OM,C_0,V_G,V_G_SQ,DV_G
      REAL(KIND=JWRB) :: XNL_1,XNL_2,XNL_3,XNL_4,XNL
      REAL(KIND=JWRB) :: C_S_SQ,ALP,ZFAC,SIG_TH,XNU

!----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('TRANSF_SNL',0,ZHOOK_HANDLE)
#endif
!*    1. DETERMINE TRANSFER FUNCTION.
!     ------------------------------

      IF(D.LT.BATHYMAX .AND. D.GT.0._JWRB) THEN
        X   = XK*D
        IF ( X .GT. DKMAX) THEN
          TRANSF_SNL = 1._JWRB
        ELSE
          X   = MAX(X,XKDMIN)
          T_0 = TANH(X)
          T_0_SQ = T_0**2
          OM  = SQRT(G*XK*T_0)
          C_0 = OM/XK
          C_S_SQ = G*D
          V_G = 0.5_JWRB*C_0*(1._JWRB+2._JWRB*X/SINH(2._JWRB*X))
          V_G_SQ = V_G**2
          DV_G = (T_0-X*(1.-T_0_SQ))**2+4._JWRB*X**2*T_0_SQ*(1._JWRB-T_0_SQ)
      
          XNL_1 = (9._JWRB*T_0_SQ**2-10._JWRB*T_0_SQ+9._JWRB)/(8._JWRB*T_0_SQ*T_0)
          XNL_2 = ((2._JWRB*V_G-0.5_JWRB*C_0)**2/(G*D-V_G_SQ)+1._JWRB)/X
          XNL_4 = 1./(4._JWRB*T_0)*(2._JWRB*C_0+V_G*(1._JWRB-T_0_SQ))**2/(C_S_SQ-V_G_SQ)
          ALP = (1.-V_G_SQ/C_S_SQ)*C_0**2/V_G_SQ
          ZFAC = SIG_TH**2/(SIG_TH**2+ALP*XNU**2)
          XNL_3 = ZFAC*XNL_4 

          XNL = XNL_1-XNL_2+XNL_3
          TRANSF_SNL = XNL**2/(DV_G*T_0_SQ**4)
          TRANSF_SNL = MAX(MIN(TRANSF_SNL_MAX,TRANSF_SNL),TRANSF_SNL_MIN)
        ENDIF
      ELSE
        TRANSF_SNL = 1._JWRB
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('TRANSF_SNL',1,ZHOOK_HANDLE)
#endif
      END FUNCTION TRANSF_SNL
