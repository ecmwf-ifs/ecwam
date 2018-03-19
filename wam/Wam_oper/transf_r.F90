      FUNCTION TRANSF_R(XK0,D)
 
!***  DETERMINE DEPTH-DEPENDENT RATIO OF ANGULAR WIDTH AND FREQUENCY WIDTH
 
!     VARIABLE       TYPE         PURPOSE
!     --------       ----         -------
 
!     XK0            REAL         WAVE NUMBER
!     D              REAL         DEPTH
 
 
!     FOR R SEE EQ. (A14) OF TECH MEMO 813 
 
!----------------------------------------------------------------------
 
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWPCONS , ONLY : G     ,DKMAX
      USE YOWSHAL , ONLY : BATHYMAX, XKDMIN
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

!----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB) :: TRANSF_R
      REAL(KIND=JWRB), PARAMETER :: EPS=0.0001_JWRB
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: XK0,D
      REAL(KIND=JWRB) :: X,XK,T_0,T_0_SQ,OM,C_0,V_G,D2OM

!----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('TRANSF_R',0,ZHOOK_HANDLE)
#endif
 
!*    1. DETERMINE TRANSFER FUNCTION.
!     ------------------------------

      IF(D.LT.BATHYMAX .AND. D.GT.0._JWRB .AND. XK0.GT.0._JWRB) THEN
        X = XK0*D
        IF ( X .GT. DKMAX) THEN
          TRANSF_R = 1._JWRB 
        ELSE
          XK  = MAX(XK0,XKDMIN/D)
          X   = XK*D
          T_0 = TANH(X)
          T_0_SQ = T_0**2
          OM  = SQRT(G*XK*T_0)
          C_0 = OM/XK
          IF(X .LT. EPS) THEN
            V_G = C_0
          ELSE
            V_G = 0.5_JWRB*C_0*(1._JWRB+2._JWRB*X/SINH(2._JWRB*X))
          ENDIF
          D2OM = (T_0-X*(1._JWRB-T_0_SQ))**2+4._JWRB*X**2*T_0_SQ*(1._JWRB-T_0_SQ)
          TRANSF_R = 4._JWRB*(V_G/C_0)**3*T_0_SQ/D2OM
        ENDIF
      ELSE
        TRANSF_R = 1._JWRB
      ENDIF
 
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('TRANSF_R',1,ZHOOK_HANDLE)
#endif
      END FUNCTION TRANSF_R
