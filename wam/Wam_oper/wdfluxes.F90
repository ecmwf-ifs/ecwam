      SUBROUTINE WDFLUXES (IJS, IJL, IG,                                &
     &                     MIJ,                                         &
     &                     FL1, XLLWS,                                  &
     &                     CICVR,                                       &
     &                     U10NEW, THWNEW, USNEW,                       &
     &                     Z0NEW, ROAIRN, WSTAR,                        &
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

!       *CALL* *WDFLUXES (FL1, IJS, IJL, IG,
!    &                    CICVR,
!    &                    THWNEW,USNEW,Z0NEW,ROAIRN,WSTAR,
!    &                    USTOKES, VSTOKES, STRNMS)
!          *FL1*    - FREQUENCY SPECTRUM(INPUT).
!          *IJS*    - INDEX OF FIRST GRIDPOINT.
!          *IJL*    - INDEX OF LAST GRIDPOINT.
!          *IG*     - BLOCK NUMBER.
!          *U10NEW* - WIND SPEED.
!          *THWNEW* - WIND DIRECTION IN RADIANS.
!          *USNEW*  - NEW FRICTION VELOCITY IN M/S.
!          *Z0NEW*  - ROUGHNESS LENGTH IN M.
!          *ROAIRN* - AIR DENSITY IN KG/M3.
!          *WSTAR*  - FREE CONVECTION VELOCITY SCALE (M/S).
!          *USTOKES*   U-COMP SURFACE STOKES DRIFT.
!          *VSTOKES*   V-COMP SURFACE STOKES DRIFT.
!          *STRNMS*    MEAN SQUARE STRAIN INTO THE SEA ICE (only if LWNEMOCOUSTRN).



!     METHOD.
!     -------

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWFLUX   ,LWVFLX_SNL, LWNEMOCOUSTRN
      USE YOWCOUT  , ONLY : LWFLUXOUT 
      USE YOWFRED  , ONLY : TH
      USE YOWICE   , ONLY : FLMIN
      USE YOWMEAN  , ONLY : PHIEPS   ,PHIAW    ,TAUOC
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "cimsstrn.intfb.h"
#include "fkmean.intfb.h"
#include "sdissip.intfb.h"
#include "sinflx.intfb.h"
#include "snonlin.intfb.h"
#include "stokesdrift.intfb.h"
#include "wnfluxes.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS,IJL,IG
      INTEGER(KIND=JWIM), INTENT(OUT) :: MIJ(IJS:IJL)

      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: CICVR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: THWNEW, ROAIRN, WSTAR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(INOUT) :: U10NEW, USNEW, Z0NEW
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: USTOKES, VSTOKES, STRNMS
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL1
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      INTEGER(KIND=JWIM) :: IJ, K, M
      INTEGER(KIND=JWIM) :: ICALL, NCALL
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL) :: MIJFLX

      REAL(KIND=JWRB) :: TAU, XN, PHIDIAG, TAUO
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: TAUW_LOC  ! TAUW should not be updated do use a local array
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: EMEANALL, FMEANALL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: FMEANWS
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: F1MEAN, AKMEAN, XKMEAN
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: PHIWA
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG) :: FLM
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE) :: RHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: FL, SL, SPOS

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

      CALL FKMEAN(FL1, IJS, IJL, EMEANALL, FMEANALL,                    &
     &            F1MEAN, AKMEAN, XKMEAN)

      DO K=1,NANG
        DO IJ=IJS,IJL
          FLM(IJ,K)=FLMIN*MAX(0.0_JWRB,COS(TH(K)-THWNEW(IJ)))**2
        ENDDO
      ENDDO

      TAUW_LOC(:) = 0.0_JWRB

      LUPDTUS = .FALSE.
      FMEANWS(:) = FMEANALL(:)
      NCALL = 1
      ICALL = 1
      CALL SINFLX (ICALL, NCALL, IJS, IJL, &
     &             LUPDTUS, &
     &             U10NEW, THWNEW, ROAIRN, WSTAR, &
     &             CICVR, &
     &             FMEANALL, FLM, &
     &             FMEANWS, FL1, &
     &             USNEW, TAUW_LOC, Z0NEW, PHIWA, &
     &             FL, SL, SPOS, &
     &             MIJ, MIJFLX, RHOWGDFTH, XLLWS)

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. WDFLUXES: SINFLX CALLED ', ICALL
        CALL FLUSH (IU06)
      ENDIF

      IF(LCFLX) THEN

        CALL SDISSIP (FL1 ,FL, SL, IJS, IJL,                            &
     &                EMEANALL, F1MEAN, XKMEAN,                         &
     &                USNEW, THWNEW, ROAIRN)
        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '   SUB. WDFLUXES: SDISSIP CALLED'
          CALL FLUSH (IU06)
        ENDIF

        IF(.NOT. LWVFLX_SNL) THEN
          CALL WNFLUXES (IJS, IJL,                                      &
     &                   MIJFLX, RHOWGDFTH,                             &
     &                   SL, CICVR,                                     &
     &                   PHIWA,                                         &
     &                   EMEANALL, F1MEAN, U10NEW, THWNEW,              &
     &                   USNEW, ROAIRN, .FALSE.)
        ENDIF

        CALL SNONLIN (FL1, FL, IJS, IJL, IG, SL, AKMEAN)
        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '   SUB. WDFLUXES: SNONLIN CALLED'
          CALL FLUSH (IU06)
        ENDIF

        IF(LWVFLX_SNL) THEN
          CALL WNFLUXES (IJS, IJL,                                      &
     &                   MIJFLX, RHOWGDFTH,                             &
     &                   SL, CICVR,                                     &
     &                   PHIWA,                                         &
     &                   EMEANALL, F1MEAN, U10NEW, THWNEW,              &
     &                   USNEW, ROAIRN, .FALSE.)
        ENDIF

        CALL STOKESDRIFT(FL1, IJS, IJL, USTOKES, VSTOKES)

        IF(LWNEMOCOUSTRN) CALL CIMSSTRN(FL1, IJS, IJL, STRNMS)

      ENDIF
! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WDFLUXES',1,ZHOOK_HANDLE)

      END SUBROUTINE WDFLUXES
