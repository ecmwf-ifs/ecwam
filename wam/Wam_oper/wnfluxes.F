      SUBROUTINE WNFLUXES (IJS, IJL,
     &                     TAUWLF, PHIWA, PHIAWDIAG, PHIAWUNR,
     &                     PHIOC, TAUWD,
     &                     EM, FM, U10, THW,
     &                     USNEW, ROAIRN, LNUPD)

! ----------------------------------------------------------------------

!**** *WNFLUXES* - FINAL STAGE OF THE WAVE FLUXES CALCULATION AND
!                  NORMALISATION.

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *WNFLUXES* (IJS, IJL,
!    &                     PHIWA, PHIAWDIAG, PHIAWUNR,
!    &                     PHIOC, TAUWD,
!    &                     EM, FM, U10, THW,
!    &                     USNEW, ROAIRN)

!          *IJS*    - INDEX OF FIRST GRIDPOINT.
!          *IJL*    - INDEX OF LAST GRIDPOINT.
!          *TAUWLF* - LOW-FREQUENCY PART OF WAVE STRESS
!                    (i.e. INTEGRATED OVER  THE PROGNOSTIC RANGE)
!          *PHIWA*  - ENERGY FLUX FROM WIND INTO WAVES INTEGRATED
!                     OVER THE FULL DISCRETISED FREQUENCY RANGE.
!          *PHIAWDIAG* - DIAGNOSTIC PART OF ENERGY FLUX INTO WAVES
!                        (i.e. INTEGRATED FROM THE END OF THE PROGNOSTIC RANGE
!                         TO THE LAST DISCRETISED FREQUENCY) 
!          *PHIAWUNR*- UNRESOLVED PART OF ENERGY FLUX
!                      (i.e. INTEGRATED FROM THE LAST DISCRETISED FREQUENCY
!                       TO INFINITY).
!          *PHIOC*  - ENERGY FLUX FROM WAVES TO OCEAN INTEGRATED OVER 
!                     THE PROGNOSTIC RANGE.
!          *TAUWD*  - DISSIPATION STRESS INTEGRATED OVER
!                     THE PROGNOSTIC RANGE.
!          *EM*     - MEAN WAVE ENERGY.
!          *FM*     - MEAN WAVE FREQUENCY BASED ON (1/f)*F INTEGRATION.
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
      USE YOWMEAN  , ONLY : PHIEPS   ,PHIAW    ,TAUOC    ,
     &                      NPHIEPS  ,NTAUOC   ,NSWH     ,NMWP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJ,IJS,IJL

      REAL :: TAU, XN, PHIDIAG, TAUO
      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(IJS:IJL) :: TAUWLF, PHIWA, PHIAWDIAG, PHIAWUNR
      REAL, DIMENSION(IJS:IJL) :: PHIOC, TAUWD
      REAL, DIMENSION(IJS:IJL) :: EM, FM, U10, THW
      REAL, DIMENSION(IJS:IJL) :: USNEW, ROAIRN
      LOGICAL :: LNUPD

      IF (LHOOK) CALL DR_HOOK('WNFLUXES',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    DETERMINE NORMALIZED FLUXES FROM AIR TO WAVE AND FROM WAVE TO OCEAN.
!     -------------------------------------------------------------------

      DO IJ=IJS,IJL
        TAU       = ROAIRN(IJ)*USNEW(IJ)**2
        XN        = ROAIRN(IJ)*USNEW(IJ)**3
        PHIDIAG   = PHIAWDIAG(IJ)+PHIAWUNR(IJ)
        PHIEPS(IJ)= (PHIOC(IJ)-PHIDIAG)/XN 
        PHIAW(IJ) = (PHIWA(IJ)+PHIAWUNR(IJ))/XN
        TAUO       = TAU-TAUWLF(IJ)-TAUWD(IJ)
        TAUOC(IJ)  = TAUO/TAU
        IF (LWNEMOCOU.AND.LNUPD) THEN
          NPHIEPS(IJ) = PHIEPS(IJ)
          NTAUOC(IJ)  = TAUOC(IJ)
          IF (EM(IJ)/=0.) THEN
             NSWH(IJ) = 4.0*SQRT(EM(IJ))
          ELSE
             NSWH(IJ) = 0.0
          ENDIF
          IF (FM(IJ)/=0.) THEN
             NMWP(IJ) = 1.0/FM(IJ)
          ELSE
             NMWP(IJ) = 0.0
          ENDIF

          IF (LWNEMOTAUOC) THEN
            NEMOTAUX(IJ)  = NEMOTAUX(IJ) + TAUO*SIN(THW(IJ))
            NEMOTAUY(IJ)  = NEMOTAUY(IJ) + TAUO*COS(THW(IJ))
          ELSE
            NEMOTAUX(IJ)  = NEMOTAUX(IJ) + TAU*SIN(THW(IJ))
            NEMOTAUY(IJ)  = NEMOTAUY(IJ) + TAU*COS(THW(IJ))
          ENDIF 
          NEMONEW10(IJ) = NEMONEW10(IJ) + U10(IJ)
          NEMOPHIF(IJ)  = NEMOPHIF(IJ) + (PHIOC(IJ)-PHIDIAG)
        ENDIF
      ENDDO

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WNFLUXES',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE WNFLUXES
