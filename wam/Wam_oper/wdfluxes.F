      SUBROUTINE WDFLUXES (IJS, IJL, IG,
     &                     FL3, FL, SL,
     &                     CICVR,
     &                     THWNEW, USNEW, Z0NEW, ROAIRN, ZIDLNEW)

! ----------------------------------------------------------------------

!**** *WDFLUXES* - IF NEEDED, IT EVALUATES SINPUT AND SDISSIP
!****              FOR THE CALCULATION OF OUTPUT WAVE DEPENDANT FLUXES.
!****              THIS IS NORMALLY DONE IN *IMPLSCH* BUT AT INITILISATION,
!****              THE SOURCES TERMS ARE NOT YET EVALUATED.

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *WDFLUXES (FL3, FL, SL, IJS, IJL, IG,
!    &                    CICVR,
!    &                    THWNEW,USNEW,Z0NEW,ROAIRN,ZIDLNEW)
!          *FL3*    - FREQUENCY SPECTRUM(INPUT).
!          *FL*     - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE (dummy)
!          *SL*     - TOTAL SOURCE FUNCTION ARRAY (dummy).
!          *IJS*    - INDEX OF FIRST GRIDPOINT.
!          *IJL*    - INDEX OF LAST GRIDPOINT.
!          *IG*     - BLOCK NUMBER.
!          *THWNEW* - WIND DIRECTION IN RADIANS.
!          *USNEW*  - NEW FRICTION VELOCITY IN M/S.
!          *Z0NEW*  - ROUGHNESS LENGTH IN M.
!          *ROAIRN* - AIR DENSITY IN KG/M3.
!          *ZIDLNEW*- Zi/L (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).

!     METHOD.
!     -------

!     EXTERNALS.
!     ---------

!       *FEMEAN*    - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT.
!       *SDISSIP*   - COMPUTATION OF DISSIPATION SOURCE FUNCTION
!                     AND LINEAR CONTRIBUTION OF DISSIPATION TO
!                     FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
!       *SEMEAN*    - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.
!       *SINPUT*    - COMPUTATION OF INPUT SOURCE FUNCTION, AND
!                     LINEAR CONTRIBUTION OF INPUT SOURCE FUNCTION
!                     TO FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
!       *STRESSO*   - COMPUTATION NORMALISED WAVE STRESS.
!           !!!!!!! MAKE SURE THAT SINPUT IS CALLED FIRST, STRESSO
!           !!!!!!! NEXT, AND THEN THE REST OF THE SOURCE FUNCTIONS.
!       *FRCUTINDEX*

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWFLUX 
      USE YOWCOUT  , ONLY : LWFLUXOUT 
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,PHIEPS   ,PHIAW    ,
     &                      TAUOC    
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      IMPLICIT NONE

      REAL,DIMENSION(IJS:IJL,NANG,NFRE) :: FL,FL3,SL
      REAL,DIMENSION(IJS:IJL) :: CICVR
      REAL,DIMENSION(IJS:IJL) :: THWNEW,USNEW,Z0NEW,
     &                           ROAIRN,ZIDLNEW

! ----------------------------------------------------------------------

      INTEGER :: IJ,IJS,IJL,IG
      INTEGER :: MIJ(IJS:IJL)

      REAL :: TAU, XN, PHIDIAG
      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(IJS:IJL) :: TAUW_LOC  ! TAUW should not be updated do use a local array
      REAL, DIMENSION(IJS:IJL) :: EMEANWS, FMEANWS
      REAL, DIMENSION(IJS:IJL) :: F1MEAN, AKMEAN, XKMEAN 
      REAL, DIMENSION(IJS:IJL) :: TAUWLF,TAUWD,PHIAWDIAG,PHIAWUNR,PHIOC,
     &                            PHIWA 
      REAL, DIMENSION(IJS:IJL,NANG,NFRE) :: XLLWS

      LOGICAL :: LCFLX

      IF (LHOOK) CALL DR_HOOK('WDFLUXES',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. INITIALISATION.
!        ---------------

      LCFLX=LWFLUX.OR.LWFLUXOUT
! ----------------------------------------------------------------------

      IF(LCFLX) THEN

!*      1.2 COMPUTE MEAN PARAMETERS.
!          ------------------------

        CALL FKMEAN(FL3, IJS, IJL, EMEAN(IJS), FMEAN(IJS),
     &              F1MEAN, AKMEAN, XKMEAN)

! ----------------------------------------------------------------------

!*      1.2 COMPUTATION OF RELEVANT SOURCE FUNCTIONS.
!           -----------------------------------------

        CALL SINPUT (FL3, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,
     &               ROAIRN, ZIDLNEW, SL, XLLWS)
        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '   SUB. WDFLUXES: SINPUT CALLED'
          CALL FLUSH (IU06)
        ENDIF

!       MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
        CALL FEMEANWS(FL3,IJS,IJL,EMEANWS,FMEANWS,XLLWS)

!       COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
        CALL FRCUTINDEX(IJS, IJL, FMEAN(IJS), FMEANWS, CICVR, MIJ)

        CALL STRESSO (FL3, IJS, IJL, THWNEW, USNEW, Z0NEW,
     &                ROAIRN, TAUW_LOC, TAUWLF, PHIWA, 
     &                PHIAWDIAG, PHIAWUNR, SL, 
     &                MIJ, LCFLX)
        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '   SUB. WDFLUXES: STRESSO CALLED'
          CALL FLUSH (IU06)
        ENDIF

        CALL SDISSIP (FL3 ,FL, IJS, IJL, IG, SL, F1MEAN, XKMEAN,
     &                PHIOC, TAUWD, MIJ, LCFLX)
        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '   SUB. WDFLUXES: SDISSIP CALLED'
          CALL FLUSH (IU06)
        ENDIF

!*      1.3 DETERMINE FLUXES FROM AIR TO WAVE AND FROM WAVE TO OCEAN.
!           -------------------------------------------------------

        DO IJ=IJS,IJL
          TAU       = ROAIRN(IJ)*USNEW(IJ)**2
          XN        = ROAIRN(IJ)*USNEW(IJ)**3
          PHIDIAG    = PHIAWDIAG(IJ)+PHIAWUNR(IJ)
          PHIEPS(IJ) = (PHIOC(IJ)-PHIDIAG)/XN 
          PHIAW(IJ)  = (PHIWA(IJ)+PHIAWUNR(IJ))/XN
          TAUOC(IJ)  = (TAU-TAUWLF(IJ)-TAUWD(IJ))/TAU
        ENDDO

      ENDIF

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WDFLUXES',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE WDFLUXES
