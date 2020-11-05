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
      USE YOWALTAS , ONLY : EGRCRV   ,AFCRV       ,BFCRV
      USE YOWCOUP  , ONLY : LWNEMOCOU, LWNEMOTAUOC, NEMOTAUX, NEMOTAUY, &
     &                      NEMONEW10, NEMOPHIF   , LWFLUX,             &
     &                      NPHIEPS  ,NTAUOC      ,NSWH     ,NMWP
      USE YOWFRED  , ONLY : FR       ,COSTH       ,SINTH    ,FRIC
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE  , LWAMRSETCI, CITHRSH
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,PHIEPS   ,PHIAW    ,    &
     &                      TAUOC    ,TAUXD    ,TAUYD    ,              &
     &                      TAUOCXD  ,TAUOCYD  ,PHIOCD
      USE PARKIND1 , ONLY : JPRO
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : TAUOCMIN ,TAUOCMAX ,PHIEPSMIN,PHIEPSMAX,    &
     &               EPSUS ,EPSU10   ,G        ,ZPI
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

!     FICTITIOUS VALUE OF THE NORMALISED WAVE ENERGY FLUX UNDER THE SEA ICE 
!     (negative because it is defined as leaving the waves)
      REAL(KIND=JWRB), PARAMETER :: PHIOC_ICE=-3.75_JWRB

!     USE HERSBACH 2011 FOR CD(U10) (SEE ALSO EDSON et al. 2013)
      REAL(KIND=JWRB), PARAMETER :: C1 = 1.03E-3_JWRB
      REAL(KIND=JWRB), PARAMETER :: C2 = 0.04E-3_JWRB
      REAL(KIND=JWRB), PARAMETER :: P1 = 1.48_JWRB
      REAL(KIND=JWRB), PARAMETER :: P2 = -0.21_JWRB
      REAL(KIND=JWRB), PARAMETER :: CDMAX = 0.003_JWRB

      REAL(KIND=JWRB), PARAMETER :: EFD_MIN = 0.0625_JWRB  ! corresponds to min Hs=1m under sea ice
      REAL(KIND=JWRB), PARAMETER :: EFD_MAX = 6.25_JWRB    ! corresponds to max Hs=10m under sea ice

      REAL(KIND=JWRB) :: TAU, XN, TAUO
      REAL(KIND=JWRB) :: U10P, CD_BULK, CD_WAVE, CD_ICE 
      REAL(KIND=JWRB) :: CNST
      REAL(KIND=JWRB) :: EPSUS3 
      REAL(KIND=JWRB) :: CITHRSH_INV
      REAL(KIND=JWRB) :: EFD, FFD, EFD_FAC, FFD_FAC
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: XSTRESS, YSTRESS
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: USTAR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: PHILF
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: OOVAL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: EM_OC, F1_OC
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: CMRHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: SUMT, SUMX, SUMY

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WNFLUXES',0,ZHOOK_HANDLE)

      EPSUS3 =  EPSUS*SQRT(EPSUS)

      CITHRSH_INV=1._JWRB/MAX(CITHRSH,0.01_JWRB)

      EFD_FAC = 4.0_JWRB*EGRCRV/G**2 
      FFD_FAC = (EGRCRV/AFCRV)**(1.0_JWRB/BFCRV) * G

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

      IF (LICERUN .AND. LMASKICE .AND. LWAMRSETCI) THEN
        DO IJ=IJS,IJL
          IF(CICVR(IJ) .GT. 0.0_JWRB) THEN
            OOVAL(IJ)=EXP(-MIN((CICVR(IJ)*CITHRSH_INV)**4,10._JWRB))
!           ADJUST USTAR FOR THE PRESENCE OF SEA ICE
            U10P = MAX(U10(IJ),EPSU10)
            CD_BULK = MIN((C1 + C2*U10P**P1)*U10P**P2, CDMAX)
            CD_WAVE = (USNEW(IJ)/U10P)**2
            CD_ICE = OOVAL(IJ)*CD_WAVE + (1.0_JWRB-OOVAL(IJ))*CD_BULK
            USTAR(IJ) = MAX(SQRT(CD_ICE)*U10P,EPSUS)

            ! EM_OC and F1_OC with fully developed model ENERGY 
            ! The significant wave height derived from EM_OC will be used
            ! by NEMO as a scaling factor as if it was open ocean
            EFD = MIN(EFD_FAC*USTAR(IJ)**4, EFD_MAX)
            EM_OC(IJ) = MAX(OOVAL(IJ)*EM(IJ)+(1.0_JWRB-OOVAL(IJ))*EFD, EFD_MIN)
            FFD = FFD_FAC/USTAR(IJ)
            F1_OC(IJ) = OOVAL(IJ)*F1(IJ) + (1.0_JWRB-OOVAL(IJ))*FFD
            F1_OC(IJ) = MIN(MAX(F1_OC(IJ), FR(2)),FR(NFRE))
          ELSE
            OOVAL(IJ) = 1.0_JWRB
            USTAR(IJ) = USNEW(IJ)
            EM_OC(IJ) = EM(IJ)
            F1_OC(IJ) = F1(IJ)
          ENDIF
        ENDDO
      ELSE
        OOVAL(:) = 1.0_JWRB
        USTAR(:) = USNEW(:)
        EM_OC(:) = EM(:)
        F1_OC(:) = F1(:)
      ENDIF

      IF(LWFLUX) THEN
        DO IJ=IJS,IJL
          EMEAN(IJ)  = EM_OC(IJ) 
          FMEAN(IJ)  = F1_OC(IJ) 
        ENDDO
      ENDIF

      DO IJ=IJS,IJL

        TAU        = ROAIRN(IJ)*MAX(USTAR(IJ)**2,EPSUS)
        TAUXD(IJ)  = TAU*SIN(THW(IJ))
        TAUYD(IJ)  = TAU*COS(THW(IJ))

        TAUOCXD(IJ)= TAUXD(IJ)-OOVAL(IJ)*XSTRESS(IJ)
        TAUOCYD(IJ)= TAUYD(IJ)-OOVAL(IJ)*YSTRESS(IJ)
        TAUO       = SQRT(TAUOCXD(IJ)**2+TAUOCYD(IJ)**2)
        TAUOC(IJ)  = MIN(MAX(TAUO/TAU,TAUOCMIN),TAUOCMAX)

        XN        = ROAIRN(IJ)*MAX(USTAR(IJ)**3,EPSUS3)
        PHIOCD(IJ)= OOVAL(IJ)*(PHILF(IJ)-PHIWA(IJ))+(1.0_JWRB-OOVAL(IJ))*PHIOC_ICE*XN

        PHIEPS(IJ)= PHIOCD(IJ)/XN 
        PHIEPS(IJ)= MIN(MAX(PHIEPS(IJ),PHIEPSMIN),PHIEPSMAX)

        PHIOCD(IJ)= PHIEPS(IJ)*XN

        PHIAW(IJ) = PHIWA(IJ)/XN

        IF (LWNEMOCOU.AND.LNUPD) THEN
          NPHIEPS(IJ) = PHIEPS(IJ)
          NTAUOC(IJ)  = TAUOC(IJ)
          IF (EM_OC(IJ)/=0.0_JWRB) THEN
             NSWH(IJ) = 4.0_JPRO*SQRT(EM_OC(IJ))
          ELSE
             NSWH(IJ) = 0.0_JPRO
          ENDIF
          IF (F1_OC(IJ)/=0.0_JWRB) THEN
             NMWP(IJ) = 1.0_JPRO/F1_OC(IJ)
          ELSE
             NMWP(IJ) = 0.0_JPRO
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
