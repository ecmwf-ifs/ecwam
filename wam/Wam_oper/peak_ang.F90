      SUBROUTINE PEAK_ANG(F1, IJS, IJL, XNU, SIG_TH)

!***  *PEAK_ANG*   DETERMINES ANGULAR WIDTH NEAR PEAK OF SPECTRUM

!     PETER JANSSEN

!     PURPOSE.
!     --------

!              DETERMINATION OF PEAK PARAMETERS

!     INTERFACE.
!     ----------
!              *CALL*  *PEAK_ANG(F1,IJS,IJL,XNU,SIG_TH)*

!               INPUT:
!                  *F1*   - FREQUENCY SPECTRUM
!                  *IJS*   - FIRST GRIDPOINT              
!                  *IJL*   - LAST GRIDPOINT              
!               OUTPUT:
!                  *XNU*    - RELATIVE SPECTRAL WIDTH
!                  *SIG_TH* - RELATIVE WIDTH IN DIRECTION

!     METHOD.
!     -------
!              NONE

!     EXTERNALS.
!     ----------
!              NONE

!-----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMFR  ,DFIMOFR ,      &
     &               DELTH, TH       ,SINTH    ,COSTH   ,WETAIL  ,      &
     &               WP1TAIL ,FRTAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
 
! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F1
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: XNU, SIG_TH

      INTEGER(KIND=JWIM), PARAMETER :: NSH = 5 
      INTEGER(KIND=JWIM) :: IJ, M, K
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL)::MMAX, MSTART, MSTOP

      REAL(KIND=JWRB), PARAMETER :: CONST_SIG = 0.86_JWRB
      REAL(KIND=JWRB) :: R1
      REAL(KIND=JWRB) :: DELT25, COEF_FR, DELT2
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: ZEPSILON
      REAL(KIND=JWRB),DIMENSION(IJS:IJL) :: SUM0, SUM1, SUM2, SUM4, SUM6, XMAX
      REAL,DIMENSION(IJS:IJL) :: TEMP
      REAL(KIND=JWRB),DIMENSION(IJS:IJL,NFRE) :: THMEAN, SUM_S, SUM_C

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('PEAK_ANG',0,ZHOOK_HANDLE)
#endif

!***  1. DETERMINE L-H SPECTRAL WIDTH OF THE 2-D SPECTRUM.
!     ---------------------------------------------------

      ZEPSILON=10._JWRB*EPSILON(ZEPSILON)

      DO IJ=IJS,IJL
        SUM0(IJ)= ZEPSILON
        SUM1(IJ)= 0._JWRB
        SUM6(IJ)= 0._JWRB          
      ENDDO

      DO M=1,NFRE
        K=1
        DO IJ=IJS,IJL
          TEMP(IJ) = F1(IJ,K,M)
        ENDDO
        DO K=2,NANG
          DO IJ=IJS,IJL
            TEMP(IJ) = TEMP(IJ)+F1(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          SUM0(IJ) = SUM0(IJ)+TEMP(IJ)*DFIM(M)
          SUM1(IJ) = SUM1(IJ)+TEMP(IJ)*DFIMFR(M)
          SUM6(IJ) = SUM6(IJ)+TEMP(IJ)*DFIMOFR(M)
        ENDDO
      ENDDO

!     ADD TAIL CORRECTIONS
      DELT25 = WETAIL*FR(NFRE)*DELTH
      COEF_FR = WP1TAIL*DELTH*FR(NFRE)**2
      DELT2 = FRTAIL*DELTH
      DO IJ=IJS,IJL
        SUM0(IJ) = SUM0(IJ)+DELT25*TEMP(IJ)
        SUM1(IJ) = SUM1(IJ)+COEF_FR*TEMP(IJ)
        SUM6(IJ) = SUM6(IJ)+DELT2*TEMP(IJ)
      ENDDO

      DO IJ=IJS,IJL
        IF (SUM0(IJ).GT.ZEPSILON) THEN
          XNU(IJ) = SQRT(MAX(0._JWRB,SUM1(IJ)*SUM6(IJ)/SUM0(IJ)**2-1._JWRB))
        ELSE
          XNU(IJ) = ZEPSILON 
        ENDIF 
      ENDDO

!***  2. DETERMINE ANGULAR WIDTH OF THE 2-D SPECTRUM.
!     ----------------------------------------------
 
!     MAX OF 2-D SPECTRUM
 
      DO IJ=IJS,IJL
        XMAX(IJ) = 0._JWRB
        MMAX(IJ) = 2
      ENDDO

      DO M=2,NFRE-1
        DO K=1,NANG
          DO IJ=IJS,IJL
            IF (F1(IJ,K,M).GT.XMAX(IJ)) THEN
              MMAX(IJ) = M
              XMAX(IJ) = F1(IJ,K,M)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        SUM1(IJ) = ZEPSILON
        SUM2(IJ) = 0._JWRB
      ENDDO
 
      DO M=1,NFRE
        DO IJ=IJS,IJL
          SUM_S(IJ,M) = 0._JWRB
          SUM_C(IJ,M) = ZEPSILON 
        ENDDO
      ENDDO

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            SUM_S(IJ,M) = SUM_S(IJ,M) +SINTH(K)*F1(IJ,K,M)
            SUM_C(IJ,M) = SUM_C(IJ,M) +COSTH(K)*F1(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

      DO M=1,NFRE
        DO IJ=IJS,IJL
          THMEAN(IJ,M) = ATAN2(SUM_S(IJ,M),SUM_C(IJ,M))
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        MSTART(IJ) = MAX(1,MMAX(IJ)-NSH) 
        MSTOP(IJ)  = MIN(NFRE,MMAX(IJ)+NSH) 
      ENDDO

      DO IJ=IJS,IJL
        DO M=MSTART(IJ),MSTOP(IJ)
          DO K=1,NANG
            SUM1(IJ) = SUM1(IJ) +F1(IJ,K,M)*DFIM(M)
            SUM2(IJ) = SUM2(IJ) +COS(TH(K)-THMEAN(IJ,M))*F1(IJ,K,M)*DFIM(M)
          ENDDO
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        IF (SUM1(IJ).GT.ZEPSILON) THEN
          R1 = SUM2(IJ)/SUM1(IJ)
          SIG_TH(IJ) = CONST_SIG*SQRT(2._JWRB*(1._JWRB-R1))
        ELSE
          SIG_TH(IJ) = 0._JWRB
        ENDIF
      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('PEAK_ANG',1,ZHOOK_HANDLE)
#endif
      END SUBROUTINE PEAK_ANG
