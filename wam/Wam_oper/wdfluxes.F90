      SUBROUTINE WDFLUXES (KIJS, KIJL,                        &
     &                     MIJ,                               &
     &                     FL1, XLLWS,                        &
     &                     WAVNUM, CINV, CGROUP, STOKFAC,     &
     &                     DEPTH,                             &
     &                     CICVR,                             &
     &                     U10NEW, THWNEW, USNEW,             &
     &                     Z0NEW, Z0B, ROAIRN, WSTAR,         &
     &                     WSEMEAN, WSFMEAN,                  &
     &                     USTOKES, VSTOKES, STRNMS )

! ----------------------------------------------------------------------

!**** *WDFLUXES* - IF NEEDED, IT EVALUATES SINPUT AND SDISSIP
!****              FOR THE CALCULATION OF OUTPUT WAVE DEPENDANT FLUXES.
!****              THIS IS NORMALLY DONE IN *IMPLSCH* BUT AT INITILISATION,
!****              THE SOURCES TERMS ARE NOT YET EVALUATED.
!****              ALSO DETERMINE THE SURFACE STOKES DRIFT ANd STRIN IN THE ICE.

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *WDFLUXES (KIJS, KIJL,
!    &                    MIJ,
!    &                    FL1, XLLWS,
!    &                    WAVNUM, CINV, CGROUP, STOKFAC,
!    &                    CICVR,
!    &                    THWNEW,USNEW,Z0NEW,Z0B,ROAIRN,WSTAR,
!    &                    WSEMEAN, WSFMEAN,
!    &                    USTOKES, VSTOKES, STRNMS)

!          *KIJS*   - INDEX OF FIRST GRIDPOINT.
!          *KIJL*   - INDEX OF LAST GRIDPOINT.
!          *FL1*    - SPECTRUM(INPUT).
!          *XLLWS* - TOTAL WINDSEA MASK FROM INPUT SOURCE TERM
!          *WAVNUM* - WAVE NUMBER.
!          *CINV*   - INVERSE PHASE VELOCITY.
!          *CGROUP* - GROUP SPPED.
!          *STOKFAC* - STOKES DRIFT FACTOR.
!          *DEPTH*  - WATER DEPTH.
!          *U10NEW* - WIND SPEED.
!          *THWNEW* - WIND DIRECTION IN RADIANS.
!          *USNEW*  - NEW FRICTION VELOCITY IN M/S.
!          *Z0NEW*  - ROUGHNESS LENGTH IN M.
!          *Z0B*    - BACKGROUND ROUGHNESS LENGTH.
!          *ROAIRN* - AIR DENSITY IN KG/M3.
!          *WSTAR*  - FREE CONVECTION VELOCITY SCALE (M/S).
!          *WSEMEAN*  WINDSEA VARIANCE.
!          *WSFMEAN*  WINDSEA MEAN FREQUENCY.
!          *USTOKES*  U-COMP SURFACE STOKES DRIFT.
!          *VSTOKES*  V-COMP SURFACE STOKES DRIFT.
!          *STRNMS*   MEAN SQUARE STRAIN INTO THE SEA ICE (only if LWNEMOCOUSTRN).


!     METHOD.
!     -------

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWFLUX   ,LWVFLX_SNL, LWNEMOCOUSTRN
      USE YOWCOUT  , ONLY : LWFLUXOUT 
      USE YOWFRED  , ONLY : FR       ,TH

      USE YOWMEAN  , ONLY : PHIEPS   ,PHIAW    ,TAUOC
!                        !!!!!!!!!!!!!!!!!!!!!!!!!

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : WSEMEAN_MIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "cimsstrn.intfb.h"
#include "femeanws.intfb.h"
#include "fkmean.intfb.h"
#include "sdissip.intfb.h"
#include "sinflx.intfb.h"
#include "snonlin.intfb.h"
#include "stokesdrift.intfb.h"
#include "wnfluxes.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM),  DIMENSION(KIJS:KIJL), INTENT(OUT) :: MIJ

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM, CINV, CGROUP, STOKFAC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: CICVR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: THWNEW, ROAIRN, WSTAR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: U10NEW, USNEW, Z0NEW, Z0B
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: WSEMEAN, WSFMEAN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: USTOKES, VSTOKES, STRNMS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      INTEGER(KIND=JWIM) :: IJ, K, M
      INTEGER(KIND=JWIM) :: ICALL, NCALL

      REAL(KIND=JWRB) :: TAU, XN, PHIDIAG, TAUO
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TAUW_LOC  ! TAUW should not be updated do use a local array
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TAUWDIR_LOC  ! TAUW should not be updated do use a local array
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EMEAN, FMEAN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EMEANWS, FMEANWS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: F1MEAN, AKMEAN, XKMEAN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: PHIWA
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: FLM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: RHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: FLD, SL, SPOS

      LOGICAL :: LCFLX
      LOGICAL :: LUPDTUS

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WDFLUXES',0,ZHOOK_HANDLE)

!*    1. INITIALISATION.
!        ---------------

      LCFLX=LWFLUX.OR.LWFLUXOUT
! ----------------------------------------------------------------------

!*    1.2 COMPUTATION OF RELEVANT SOURCE FUNCTIONS.
!         -----------------------------------------

      CALL FKMEAN(KIJS, KIJL, FL1, WAVNUM,                         &
     &            EMEAN, FMEAN, F1MEAN, AKMEAN, XKMEAN)

      TAUW_LOC(:) = 0.0_JWRB
      TAUWDIR_LOC(:) = THWNEW(:)

      LUPDTUS = .FALSE.
      FMEANWS(:) = FMEAN(:)
      FLM(:,:) = 0.0_JWRB
      NCALL = 1
      ICALL = 1
      CALL SINFLX (ICALL, NCALL, KIJS, KIJL,                        &
     &             LUPDTUS,                                         &
     &             FL1,                                             &
     &             WAVNUM, CINV, CGROUP,                            &
     &             U10NEW, THWNEW, ROAIRN, WSTAR, CICVR,            &
     &             FMEAN, FMEANWS,                                  &
     &             FLM,                                             &
     &             USNEW, TAUW_LOC, TAUWDIR_LOC, Z0NEW, Z0B, PHIWA, &
     &             FLD, SL, SPOS,                                   &
     &             MIJ, RHOWGDFTH, XLLWS)

      IF (LCFLX) THEN

        CALL SDISSIP (KIJS, KIJL, FL1 ,FLD, SL,    &
     &                WAVNUM, CGROUP,              &
     &                EMEAN, F1MEAN, XKMEAN,       &
     &                USNEW, THWNEW, ROAIRN)

        IF (.NOT. LWVFLX_SNL) THEN
          CALL WNFLUXES (KIJS, KIJL,                        &
     &                   MIJ, RHOWGDFTH,                    &
     &                   CINV,                              &
     &                   SL, CICVR,                         &
     &                   PHIWA,                             &
     &                   EMEAN, F1MEAN, U10NEW, THWNEW,     &
     &                   USNEW, ROAIRN, .FALSE.)
        ENDIF

        CALL SNONLIN (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN)

        IF (LWVFLX_SNL) THEN
          CALL WNFLUXES (KIJS, KIJL,                       &
     &                   MIJ, RHOWGDFTH,                   &
     &                   SL, CICVR,                        &
     &                   PHIWA,                            &
     &                   EMEAN, F1MEAN, U10NEW, THWNEW,    &
     &                   USNEW, ROAIRN, .FALSE.)
        ENDIF

        IF (LWFLUX) THEN
         CALL FEMEANWS(KIJS, KIJL, FL1, XLLWS, EMEANWS, FMEANWS)

          DO IJ=KIJS,KIJL
            IF (EMEANWS(IJ) < WSEMEAN_MIN) THEN
              WSEMEAN(IJ) = WSEMEAN_MIN 
              WSFMEAN(IJ) = 2._JWRB*FR(NFRE)
            ELSE
              WSEMEAN(IJ) = EMEANWS(IJ)
              WSFMEAN(IJ) = FMEANWS(IJ) 
            ENDIF
          ENDDO
        ENDIF

        CALL STOKESDRIFT(KIJS, KIJL, FL1, STOKFAC, U10NEW, THWNEW, CICVR, USTOKES, VSTOKES)

        IF (LWNEMOCOUSTRN) CALL CIMSSTRN(KIJS, KIJL, FL1, WAVNUM, DEPTH, STRNMS)

      ENDIF
! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WDFLUXES',1,ZHOOK_HANDLE)

      END SUBROUTINE WDFLUXES

