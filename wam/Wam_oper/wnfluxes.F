      SUBROUTINE WNFLUXES (IJS, IJL,
     &                     MIJ, RHOWGDFTH,
     &                     SSURF, CICVR,
     &                     PHIWA,
     &                     EM, F1, U10, THW,
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

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWCOUP  , ONLY : LWNEMOCOU, LWNEMOTAUOC, NEMOTAUX, NEMOTAUY,
     &                      NEMONEW10, NEMOPHIF
      USE YOWFRED  , ONLY : COSTH    ,SINTH
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE  ,CITHRSH
      USE YOWMEAN  , ONLY : PHIEPS   ,PHIAW    ,TAUOC    , 
     &                      NPHIEPS  ,NTAUOC   ,NSWH     ,NMWP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : PHIEPSMIN,PHIEPSMAX
      USE YOWSHAL  , ONLY : CINV     ,INDEP
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS,IJL
      INTEGER, DIMENSION(IJS:IJL), INTENT(IN) :: MIJ

      REAL,DIMENSION(IJS:IJL,NFRE), INTENT(IN) :: RHOWGDFTH
      REAL, DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: SSURF
      REAL, DIMENSION(IJS:IJL), INTENT(IN) :: CICVR 
      REAL, DIMENSION(IJS:IJL), INTENT(IN) :: PHIWA
      REAL, DIMENSION(IJS:IJL), INTENT(IN) :: EM, F1, U10, THW
      REAL, DIMENSION(IJS:IJL), INTENT(IN) :: USNEW, ROAIRN

      LOGICAL, INTENT(IN) :: LNUPD

      INTEGER :: IJ, K, M

      REAL, PARAMETER :: XMMIN=0.5
      REAL, PARAMETER :: XMMAX=3.0
      REAL :: TAU, XN, TAUO, PHIOCD
      REAL :: TAUX, TAUY, TAUOX, TAUOY, CNST
      REAL :: ZHOOK_HANDLE

      REAL, DIMENSION(IJS:IJL) :: XSTRESS, YSTRESS
      REAL, DIMENSION(IJS:IJL) :: PHILF
      REAL, DIMENSION(IJS:IJL) :: REDCI, XM 
      REAL, DIMENSION(IJS:IJL) :: CMRHOWGDFTH
      REAL, DIMENSION(IJS:IJL) :: SUMT, SUMX, SUMY

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WNFLUXES',0,ZHOOK_HANDLE)
#endif

!*    DETERMINE NORMALIZED FLUXES FROM AIR TO WAVE AND FROM WAVE TO OCEAN.
!     -------------------------------------------------------------------

!     ENERGY FLUX from SSURF
!     MOMENTUM FLUX FROM SSURF
      DO IJ=IJS,IJL
        PHILF(IJ) = 0.
        XSTRESS(IJ) = 0.
        YSTRESS(IJ) = 0.
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
            REDCI(IJ)=MAX(1.-CICVR(IJ),0.)
!           FICTITIOUS VALUE OF THE WAVE ENERGY FLUX UNDER THE SEA ICE 
            XM(IJ)=XMMIN+(XMMAX-XMMIN)*REDCI(IJ)
          ELSE
            REDCI(IJ)=1.
            XM(IJ)=0.
          ENDIF
        ENDDO
      ELSE
        REDCI(:)=1.
        XM(:)=0.
      ENDIF

      DO IJ=IJS,IJL

        TAU        = ROAIRN(IJ)*USNEW(IJ)**2
        TAUX       = TAU*SIN(THW(IJ))
        TAUY       = TAU*COS(THW(IJ))

        TAUOX      = TAUX-REDCI(IJ)*XSTRESS(IJ)
        TAUOY      = TAUY-REDCI(IJ)*YSTRESS(IJ)
        TAUO       = SQRT(TAUOX**2+TAUOY**2)
        TAUOC(IJ)  = TAUO/TAU

        XN        = ROAIRN(IJ)*USNEW(IJ)**3
        PHIOCD    = REDCI(IJ)*(PHILF(IJ)-PHIWA(IJ))+XM(IJ)*XN

        PHIEPS(IJ)= PHIOCD/XN 
        PHIEPS(IJ)= MIN(MAX(PHIEPS(IJ),PHIEPSMIN),PHIEPSMAX)

        PHIAW(IJ) = PHIWA(IJ)/XN

        IF (LWNEMOCOU.AND.LNUPD) THEN
          NPHIEPS(IJ) = PHIEPS(IJ)
          NTAUOC(IJ)  = TAUOC(IJ)
          IF (EM(IJ)/=0.) THEN
             NSWH(IJ) = 4.0*SQRT(EM(IJ))
          ELSE
             NSWH(IJ) = 0.0
          ENDIF
          IF (F1(IJ)/=0.) THEN
             NMWP(IJ) = 1.0/F1(IJ)
          ELSE
             NMWP(IJ) = 0.0
          ENDIF

          IF (LWNEMOTAUOC) THEN
            NEMOTAUX(IJ)  = NEMOTAUX(IJ) + TAUOX
            NEMOTAUY(IJ)  = NEMOTAUY(IJ) + TAUOY
          ELSE
            NEMOTAUX(IJ)  = NEMOTAUX(IJ) + TAUX
            NEMOTAUY(IJ)  = NEMOTAUY(IJ) + TAUY
          ENDIF 
          NEMONEW10(IJ) = NEMONEW10(IJ) + U10(IJ)
          NEMOPHIF(IJ)  = NEMOPHIF(IJ) + PHIOCD 
        ENDIF
      ENDDO

! ----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WNFLUXES',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE WNFLUXES
