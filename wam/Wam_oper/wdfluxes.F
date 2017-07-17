      SUBROUTINE WDFLUXES (IJS, IJL, IG,
     &                     MIJ,
     &                     FL3, XLLWS,
     &                     CICVR,
     &                     U10NEW, THWNEW, USNEW,
     &                     Z0NEW, ROAIRN, WSTAR)

! ----------------------------------------------------------------------

!**** *WDFLUXES* - IF NEEDED, IT EVALUATES SINPUT AND SDISSIP
!****              FOR THE CALCULATION OF OUTPUT WAVE DEPENDANT FLUXES.
!****              THIS IS NORMALLY DONE IN *IMPLSCH* BUT AT INITILISATION,
!****              THE SOURCES TERMS ARE NOT YET EVALUATED.

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *WDFLUXES (FL3, IJS, IJL, IG,
!    &                    CICVR,
!    &                    THWNEW,USNEW,Z0NEW,ROAIRN,WSTAR)
!          *FL3*    - FREQUENCY SPECTRUM(INPUT).
!          *IJS*    - INDEX OF FIRST GRIDPOINT.
!          *IJL*    - INDEX OF LAST GRIDPOINT.
!          *IG*     - BLOCK NUMBER.
!          *U10NEW* - WIND SPEED.
!          *THWNEW* - WIND DIRECTION IN RADIANS.
!          *USNEW*  - NEW FRICTION VELOCITY IN M/S.
!          *Z0NEW*  - ROUGHNESS LENGTH IN M.
!          *ROAIRN* - AIR DENSITY IN KG/M3.
!          *WSTAR*  - FREE CONVECTION VELOCITY SCALE (M/S).


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

      USE YOWCOUP  , ONLY : LWFLUX   ,LWVFLX_SNL
      USE YOWCOUT  , ONLY : LWFLUXOUT 
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWMEAN  , ONLY : PHIEPS   ,PHIAW    ,TAUOC
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS,IJL,IG
      INTEGER, INTENT(OUT) :: MIJ(IJS:IJL)

      REAL,DIMENSION(IJS:IJL),INTENT(IN) :: CICVR
      REAL,DIMENSION(IJS:IJL),INTENT(IN) :: U10NEW, THWNEW, ROAIRN,WSTAR
      REAL,DIMENSION(IJS:IJL),INTENT(INOUT) :: USNEW, Z0NEW
      REAL,DIMENSION(IJS:IJL,NANG,NFRE),INTENT(IN) :: FL3
      REAL,DIMENSION(IJS:IJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      INTEGER :: IJ, K, M

      REAL :: TAU, XN, PHIDIAG, TAUO
      REAL :: ZHOOK_HANDLE

      REAL, DIMENSION(IJS:IJL) :: TAUW_LOC  ! TAUW should not be updated do use a local array
      REAL, DIMENSION(IJS:IJL) :: EMEANALL, FMEANALL
      REAL, DIMENSION(IJS:IJL) :: EMEANWS, FMEANWS
      REAL, DIMENSION(IJS:IJL) :: F1MEAN, AKMEAN, XKMEAN
      REAL, DIMENSION(IJS:IJL) :: PHIWA
      REAL, DIMENSION(IJS:IJL,NFRE) :: RHOWGDFTH
      REAL, DIMENSION(IJS:IJL,NANG,NFRE) :: FL, SL

!????
      REAL, DIMENSION(IJS:IJL,NANG,NFRE) :: SSOURCE 

      LOGICAL :: LCFLX

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WDFLUXES',0,ZHOOK_HANDLE)
#endif
! ----------------------------------------------------------------------

!*    1. INITIALISATION.
!        ---------------

      LCFLX=LWFLUX.OR.LWFLUXOUT
! ----------------------------------------------------------------------

      IF(LCFLX) THEN

!*      1.2 COMPUTE MEAN PARAMETERS.
!          ------------------------

        CALL FKMEAN(FL3, IJS, IJL, EMEANALL, FMEANALL,
     &              F1MEAN, AKMEAN, XKMEAN)

! ----------------------------------------------------------------------

!*      1.2 COMPUTATION OF RELEVANT SOURCE FUNCTIONS.
!           -----------------------------------------

        CALL SINPUT (FL3, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,
     &               ROAIRN, WSTAR, SL, XLLWS)
        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '   SUB. WDFLUXES: SINPUT CALLED'
          CALL FLUSH (IU06)
        ENDIF
  
       IF(.NOT. LWVFLX_SNL) THEN
!      in order to reproduce the wrong way to compute the wave fluxes, we need the wind input source term
         DO M=1,NFRE
           DO K=1,NANG
             DO IJ=IJS,IJL
               SSOURCE(IJ,K,M) = SL(IJ,K,M)
             ENDDO
           ENDDO
         ENDDO
       ENDIF


!       MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
        CALL FEMEANWS(FL3,IJS,IJL,EMEANWS,FMEANWS,XLLWS)

!       COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
        CALL FRCUTINDEX(IJS, IJL, FMEANALL, FMEANWS,USNEW, CICVR,
     &                  MIJ, RHOWGDFTH)

        CALL STRESSO (FL3, SL, IJS, IJL,
     &                MIJ, RHOWGDFTH,
     &                THWNEW, USNEW, Z0NEW, ROAIRN,
     &                TAUW_LOC, PHIWA)

        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '   SUB. WDFLUXES: STRESSO CALLED'
          CALL FLUSH (IU06)
        ENDIF

        CALL SDISSIP (FL3 ,FL, SL, IJS, IJL,
     &                EMEANALL, F1MEAN, XKMEAN,
     &                USNEW, THWNEW, ROAIRN)
        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '   SUB. WDFLUXES: SDISSIP CALLED'
          CALL FLUSH (IU06)
        ENDIF

        IF(.NOT. LWVFLX_SNL) THEN
          CALL WRONG_WNFLUXES (IJS, IJL,
     &                         MIJ, RHOWGDFTH,
     &                         SSOURCE, SL,
     &                         PHIWA,
     &                         EMEANALL, F1MEAN, U10NEW, THWNEW,
     &                         USNEW, ROAIRN, .FALSE.)
        ENDIF

        CALL SNONLIN (FL3, FL, IJS, IJL, IG, SL, AKMEAN)
        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '   SUB. WDFLUXES: SNONLIN CALLED'
          CALL FLUSH (IU06)
        ENDIF

        IF(LWVFLX_SNL) THEN
          CALL WNFLUXES (IJS, IJL,
     &                   MIJ, RHOWGDFTH,
     &                   SL,
     &                   PHIWA,
     &                   EMEANALL, F1MEAN, U10NEW, THWNEW,
     &                   USNEW, ROAIRN, .FALSE.)
        ENDIF

      ENDIF
! ----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WDFLUXES',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE WDFLUXES
