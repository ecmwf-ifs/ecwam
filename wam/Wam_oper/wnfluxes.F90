      SUBROUTINE WNFLUXES (IJS, IJL,                                    &
     &                     MIJ, RHOWGDFTH,                              &
     &                     SSURF, CICVR,                                &
     &                     PHIWA,                                       &
     &                     EM, F1, U10, THW,                            &
     &                     USNEW, ROAIRN, LNUPD)

! ----------------------------------------------------------------------

!**** *WNFLUXES* - WAVE FLUXES CALCULATION

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *WNFLUXES* (IJS, IJL,
!    &                     MIJ, RHOWGDFTH,
!    &                     SSURF, CICVR,
!    &                     PHIWA,
!    &                     EM, F1, U10, THW,
!    &                     USNEW, ROAIRN, LNUPD)
!          *IJS*    - INDEX OF FIRST GRIDPOINT.
!          *IJL*    - INDEX OF LAST GRIDPOINT.
!         *MIJ*         - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!         *RHOWGDFTH    - WATER DENSITY * G * DF * DTHETA
!                         FOR TRAPEZOIDAL INTEGRATION BETWEEN FR(1) and FR(MIJ) 
!                         !!!!!!!!  RHOWGDFTH=0 FOR FR > FR(MIJ)
!          *SSURF*  - CONTRIBUTION OF ALL SOURCE TERMS ACTING ON 
!                     THE SURFACE MOMENTUM AND ENERGY FLUXES.
!          *CICVR*  - SEA ICE COVER.
!          *PHIWA*  - ENERGY FLUX FROM WIND INTO WAVES INTEGRATED
!                     OVER THE FULL FREQUENCY RANGE.
!          *EM*     - MEAN WAVE ENERGY.
!          *F1*     - MEAN WAVE FREQUENCY BASED ON f*F INTEGRATION.
!          *U10*    - NEW WIND SPEED IN M/S.
!          *THW*    - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC CONVENTION
!          *USNEW*  - NEW FRICTION VELOCITY IN M/S.
!          *ROAIRN* - AIR DENSITY IN KG/M3.
!          *LNUPD*  - UPDATE NEMO FIELDS.


!     METHOD.
!     -------

!     EXTERNALS.
!     ---------

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWCOUP  , ONLY : LWNEMOCOU, LWNEMOTAUOC, NEMOTAUX, NEMOTAUY, &
     &                      NEMONEW10, NEMOPHIF
      USE YOWFRED  , ONLY : COSTH    ,SINTH
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE  ,CITHRSH
      USE YOWMEAN  , ONLY : PHIEPS   ,PHIAW    ,TAUOC    ,              &
     &                      TAUXD    ,TAUYD    ,                        &
     &                      TAUOCXD  ,TAUOCYD  ,PHIOCD   ,              &
     &                      NPHIEPS  ,NTAUOC   ,NSWH     ,NMWP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : PHIEPSMIN,PHIEPSMAX
      USE YOWSHAL  , ONLY : CINV     ,INDEP
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS,IJL
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL), INTENT(IN) :: MIJ

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE), INTENT(IN) :: RHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: SSURF
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: CICVR 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: PHIWA
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: EM, F1, U10, THW
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: USNEW, ROAIRN

      LOGICAL, INTENT(IN) :: LNUPD

      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JWRB), PARAMETER :: XMMIN=0.5_JWRB
      REAL(KIND=JWRB), PARAMETER :: XMMAX=3.0_JWRB
      REAL(KIND=JWRB) :: TAU, XN, TAUO
      REAL(KIND=JWRB) :: CNST
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: XSTRESS, YSTRESS
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: PHILF
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: REDCI, XM
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: CMRHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: SUMT, SUMX, SUMY

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WNFLUXES',0,ZHOOK_HANDLE)


!*    DETERMINE NORMALIZED FLUXES FROM AIR TO WAVE AND FROM WAVE TO OCEAN.
!     -------------------------------------------------------------------

!     ENERGY FLUX from SSURF
!     MOMENTUM FLUX FROM SSURF
      DO IJ=IJS,IJL
        PHILF(IJ) = 0.0_JWRB
        XSTRESS(IJ) = 0.0_JWRB
        YSTRESS(IJ) = 0.0_JWRB
      ENDDO

!     THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
      DO M=1,MAXVAL(MIJ(:))
        K=1
        DO IJ=IJS,IJL
          SUMT(IJ) = SSURF(IJ,K,M)
          SUMX(IJ) = SINTH(K)*SSURF(IJ,K,M)
          SUMY(IJ) = COSTH(K)*SSURF(IJ,K,M)
        ENDDO
        DO K=2,NANG
          DO IJ=IJS,IJL
            SUMT(IJ) = SUMT(IJ) + SSURF(IJ,K,M)
            SUMX(IJ) = SUMX(IJ) + SINTH(K)*SSURF(IJ,K,M)
            SUMY(IJ) = SUMY(IJ) + COSTH(K)*SSURF(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          PHILF(IJ)   = PHILF(IJ)   + SUMT(IJ)*RHOWGDFTH(IJ,M)
          CMRHOWGDFTH(IJ) = CINV(INDEP(IJ),M)*RHOWGDFTH(IJ,M)
          XSTRESS(IJ) = XSTRESS(IJ) + SUMX(IJ)*CMRHOWGDFTH(IJ)
          YSTRESS(IJ) = YSTRESS(IJ) + SUMY(IJ)*CMRHOWGDFTH(IJ)
        ENDDO
      ENDDO

      IF (LICERUN .AND. LMASKICE) THEN
        DO IJ=IJS,IJL
          IF(CICVR(IJ).GT.CITHRSH) THEN
!           PERCENTAGE OF THE WAVE FLUXES THAT ARE ACCOUNTED FOR UNDER SEA ICE
            REDCI(IJ)=MAX(1.0_JWRB-CICVR(IJ),0.0_JWRB)
!           FICTITIOUS VALUE OF THE NORMALISED WAVE ENERGY FLUX UNDER THE SEA ICE 
!           (negative because it is defined as leaving the waves)
            XM(IJ)=-XMMIN-(XMMAX-XMMIN)*REDCI(IJ)
          ELSE
            REDCI(IJ)=1.0_JWRB
            XM(IJ)=0.0_JWRB
          ENDIF
        ENDDO
      ELSE
        REDCI(:)=1.0_JWRB
        XM(:)=0.0_JWRB
      ENDIF

      DO IJ=IJS,IJL

        TAU        = ROAIRN(IJ)*USNEW(IJ)**2
        TAUXD(IJ)  = TAU*SIN(THW(IJ))
        TAUYD(IJ)  = TAU*COS(THW(IJ))

        TAUOCXD(IJ)= TAUXD(IJ)-REDCI(IJ)*XSTRESS(IJ)
        TAUOCYD(IJ)= TAUYD(IJ)-REDCI(IJ)*YSTRESS(IJ)
        TAUO       = SQRT(TAUOCXD(IJ)**2+TAUOCYD(IJ)**2)
        TAUOC(IJ)  = TAUO/TAU

        XN        = ROAIRN(IJ)*USNEW(IJ)**3
        PHIOCD(IJ)= REDCI(IJ)*(PHILF(IJ)-PHIWA(IJ))+XM(IJ)*XN

        PHIEPS(IJ)= PHIOCD(IJ)/XN 
        PHIEPS(IJ)= MIN(MAX(PHIEPS(IJ),PHIEPSMIN),PHIEPSMAX)

        PHIAW(IJ) = PHIWA(IJ)/XN

        IF (LWNEMOCOU.AND.LNUPD) THEN
          NPHIEPS(IJ) = PHIEPS(IJ)
          NTAUOC(IJ)  = TAUOC(IJ)
          IF (EM(IJ)/=0.0_JWRB) THEN
             NSWH(IJ) = 4.0_JWRB*SQRT(EM(IJ))
          ELSE
             NSWH(IJ) = 0.0_JWRB
          ENDIF
          IF (F1(IJ)/=0.0_JWRB) THEN
             NMWP(IJ) = 1.0_JWRB/F1(IJ)
          ELSE
             NMWP(IJ) = 0.0_JWRB
          ENDIF

          IF (LWNEMOTAUOC) THEN
            NEMOTAUX(IJ) = NEMOTAUX(IJ) + TAUOCXD(IJ)
            NEMOTAUY(IJ) = NEMOTAUY(IJ) + TAUOCYD(IJ)
          ELSE
            NEMOTAUX(IJ) = NEMOTAUX(IJ) + TAUXD(IJ)
            NEMOTAUY(IJ) = NEMOTAUY(IJ) + TAUYD(IJ)
          ENDIF 
          NEMONEW10(IJ) = NEMONEW10(IJ) + U10(IJ)
          NEMOPHIF(IJ)  = NEMOPHIF(IJ) + PHIOCD(IJ) 
        ENDIF
      ENDDO

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WNFLUXES',1,ZHOOK_HANDLE)

      END SUBROUTINE WNFLUXES
