      SUBROUTINE STAT_NL(IJS, IJL,                                      &
     &                   XM0, OM0, XK0, BF2, XNU, SIG_TH, DPTH,         &
     &                   C3, C4, ETA_M, C4_B, C4_DYN)
 
!***  DETERMINE SKEWNESS AND ENVELOPE KURTOSIS FOR A NARROW-BAND WAVE 
!     TRAIN IN ARBITRARY DEPTH.
 
!     PURPOSE:
!     -------
     
!     VARIABLE       TYPE         PURPOSE
!     --------       ----         -------
 
!     XM0            REAL         WAVE VARIANCE
!     XK0            REAL         PEAK WAVE NUMBER
!     DEPTH          REAL         DEPTH
 
!     OUTPUT:
!     ------
 
!     C3             REAL         SKEWNESS
!     C4             REAL         ENVELOPE KURTOSIS
!     ETA_M          REAL         WAVE-INDUCED MEAN SURFACE ELEVATION
!     XNU            REAL         RELATIVE SPECTRAL WIDTH
!     SIG_TH         REAL         RELATIVE WIDTH in DIRECTION

 
!     AUTHOR:
!     ------
 
!     P.A.E.M. JANSSEN, NOVEMBER 2017
 
!----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWPCONS , ONLY : G        ,PI       ,C4MAX    ,C4MIN    ,DKMAX
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

!----------------------------------------------------------------------
      IMPLICIT NONE 

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(IN) :: XM0, OM0, XK0, BF2, XNU, SIG_TH, DPTH
      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(OUT) :: C3, C4, ETA_M, C4_B, C4_DYN

      REAL(KIND=JWRB), PARAMETER :: EPS = 0.0001_JWRB  
      REAL(KIND=JWRB), PARAMETER :: CONST_C3 = 1.12_JWRB*2._JWRB
      REAL(KIND=JWRB), PARAMETER :: CONST_C4 = 0.91_JWRB*8._JWRB  
      REAL(KIND=JWRB), PARAMETER :: CONST_THR = 1.1_JWRB
      REAL(KIND=JWRB), PARAMETER :: SQRT3 = SQRT(3._JWRB)

      INTEGER(KIND=JWIM) :: IJ

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: X,XK,D,T0,T0_SQ,OM,ALPH,GAM,DELTA,DELTA_1D
      REAL(KIND=JWRB) :: DELTA_2D,C_0,C_S_SQ,V_G,V_G_SQ,ZFAC,ZFAC1,ZFAC2
      REAL(KIND=JWRB) :: XKAPPA1,ALPHA,R_0,BETA_R0,XJ
      REAL(KIND=JWRB) :: ZEPSILON, ZSQREPSILON
      REAL(KIND=JWRB) :: TRANSF_R
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: TRANSF, R 

!-----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('STAT_NL',0,ZHOOK_HANDLE)
#endif

!     CONSTANTS.

      ZEPSILON=10._JWRB*EPSILON(ZEPSILON)
      ZSQREPSILON=SQRT(ZEPSILON)
 
      R_0      = 7.44_JWRB*SQRT3/(4._JWRB*PI**3)
      BETA_R0  = R_0*PI/(3._JWRB*SQRT3)
 
!     RESULTS FOR A NARROW BAND WAVE TRAIN
  
      DO IJ = IJS,IJL
        TRANSF(IJ) = TRANSF_R(XK0(IJ),DPTH(IJ))
      ENDDO

      DO IJ = IJS,IJL
        OM  = OM0(IJ)
        XK  = XK0(IJ)
        IF (XM0(IJ).GT.ZEPSILON .AND. XK.GT.0._JWRB) THEN
          D   = DPTH(IJ)
          X   = XK*D
          T0  = TANH(X)
          T0_SQ = T0**2
          ALPH = XK/(4._JWRB*T0_SQ*T0)*(3._JWRB-T0_SQ)
          GAM = -0.5_JWRB*ALPH**2
          C_0      = OM/XK
          C_S_SQ   = G*D
          IF(X .GT. DKMAX) THEN
            V_G = 0.5_JWRB*C_0
          ELSE IF(X .LT. EPS) THEN
            V_G = C_0
          ELSE
            V_G = 0.5_JWRB*C_0*(1._JWRB+2._JWRB*X/SINH(2._JWRB*X))
          ENDIF
          V_G_SQ=V_G**2

          ZFAC     = -0.25_JWRB*XK*C_S_SQ/(C_S_SQ-V_G_SQ)
          DELTA_1D = ZFAC*(2._JWRB*(1._JWRB-T0_SQ)/T0+1._JWRB/X)

          ZFAC1    = 0.5_JWRB*C_0*C_S_SQ*V_G/T0

          XKAPPA1  = ZFAC1*(2._JWRB*C_0+V_G*(1._JWRB-T0_SQ))/(C_S_SQ-V_G_SQ)
          ALPHA    = (1._JWRB-V_G_SQ/C_S_SQ)*C_0**2/V_G_SQ
          ZFAC2    = SIG_TH(IJ)**2/(SIG_TH(IJ)**2+ALPHA*XNU(IJ)**2)
          DELTA_2D = 0.5_JWRB*XK**2*XKAPPA1/(OM*C_S_SQ)*ZFAC2
  
          DELTA    = DELTA_1D+DELTA_2D
 
!         WAVE-INDUCED MEAN SEA LEVEL
 
          ETA_M(IJ) = 2._JWRB*XM0(IJ)*DELTA
 
!         SKEWNESS
 
          C3(IJ) = CONST_C3*SQRT(XM0(IJ))*(ALPH+DELTA)
 
!         BOUND KURTOSIS
 
          C4_B(IJ) = CONST_C4*XM0(IJ)*(GAM+ALPH**2+(ALPH+DELTA)**2)
 
!         DYNAMIC KURTOSIS
 
          R(IJ) = TRANSF(IJ)*(SIG_TH(IJ)/XNU(IJ))**2
          IF (R(IJ).GT.1._JWRB) THEN
            XJ = BETA_R0*(1._JWRB-CONST_THR*R(IJ))/(R(IJ)*(1._JWRB+R_0*R(IJ)))
          ELSE
            XJ = BETA_R0*(CONST_THR-R(IJ))/(R(IJ)+R_0)
          ENDIF

          C4_DYN(IJ)  = XJ*BF2(IJ)
 
!         DETERMINE TOTAL KURTOSIS AND PROTECT
 
          C4(IJ)  = C4_DYN(IJ)+C4_B(IJ)
          C4(IJ)  = MAX(MIN(C4MAX,C4(IJ)),C4MIN)
        ELSE
          ETA_M(IJ)  = 0._JWRB
          C3(IJ)     = 0._JWRB
          C4(IJ)     = 0._JWRB
          C4_DYN(IJ) = 0._JWRB
          C4_B(IJ)   = 0._JWRB
        ENDIF
      ENDDO
 
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('STAT_NL',1,ZHOOK_HANDLE)
#endif
      END SUBROUTINE STAT_NL
