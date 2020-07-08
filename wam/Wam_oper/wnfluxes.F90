      SUBROUTINE WNFLUXES (IJS, IJL,                                    &
     &                     MIJ, RHOWGDFTH,                           &
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
!          *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE for flux calculation
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
     &                      NEMONEW10, NEMOPHIF   , LWFLUX,             &
     &                      NPHIEPS  ,NTAUOC      ,NSWH     ,NMWP
      USE YOWFRED  , ONLY : COSTH    ,SINTH
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,PHIEPS   ,PHIAW    ,    &
     &                      TAUOC    ,TAUXD    ,TAUYD    ,              &
     &                      TAUOCXD  ,TAUOCYD  ,PHIOCD
      USE YOWNEMOP , ONLY : NEMODP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : TAUOCMIN ,TAUOCMAX ,PHIEPSMIN,PHIEPSMAX, EPSUS
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

      REAL(KIND=JWRB) :: TAU, XN, TAUO
      REAL(KIND=JWRB) :: CNST
      REAL(KIND=JWRB) :: EPSUS3 
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: XSTRESS, YSTRESS
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: PHILF
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: CMRHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: SUMT, SUMX, SUMY

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WNFLUXES',0,ZHOOK_HANDLE)


      EPSUS3 =  EPSUS*SQRT(EPSUS)

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

      IF(LWFLUX) THEN
        DO IJ=IJS,IJL
          EMEAN(IJ)  = EM(IJ) 
          FMEAN(IJ)  = F1(IJ) 
        ENDDO
      ENDIF

      DO IJ=IJS,IJL

        TAU        = ROAIRN(IJ)*MAX(USNEW(IJ)**2,EPSUS)
        TAUXD(IJ)  = TAU*SIN(THW(IJ))
        TAUYD(IJ)  = TAU*COS(THW(IJ))

        TAUOCXD(IJ)= TAUXD(IJ)-XSTRESS(IJ)
        TAUOCYD(IJ)= TAUYD(IJ)-YSTRESS(IJ)
        TAUO       = SQRT(TAUOCXD(IJ)**2+TAUOCYD(IJ)**2)
        TAUOC(IJ)  = MIN(MAX(TAUO/TAU,TAUOCMIN),TAUOCMAX)

        XN        = ROAIRN(IJ)*MAX(USNEW(IJ)**3,EPSUS3)
        PHIOCD(IJ)= PHILF(IJ)-PHIWA(IJ)

        PHIEPS(IJ)= PHIOCD(IJ)/XN 
        PHIEPS(IJ)= MIN(MAX(PHIEPS(IJ),PHIEPSMIN),PHIEPSMAX)

        PHIOCD(IJ)= PHIEPS(IJ)*XN

        PHIAW(IJ) = PHIWA(IJ)/XN

        IF (LWNEMOCOU.AND.LNUPD) THEN
          NPHIEPS(IJ) = PHIEPS(IJ)
          NTAUOC(IJ)  = TAUOC(IJ)
          IF (EM(IJ)/=0.0_JWRB) THEN
             NSWH(IJ) = 4.0_NEMODP*SQRT(EM(IJ))
          ELSE
             NSWH(IJ) = 0.0_NEMODP
          ENDIF
          IF (F1(IJ)/=0.0_JWRB) THEN
             NMWP(IJ) = 1.0_NEMODP/F1(IJ)
          ELSE
             NMWP(IJ) = 0.0_NEMODP
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
