      SUBROUTINE WDFLUXES (IJS, IJL, IG,                                &
     &                     MIJ,                                         &
     &                     FL3, XLLWS,                                  &
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

!       *CALL* *WDFLUXES (FL3, IJS, IJL, IG,
!    &                    CICVR,
!    &                    THWNEW,USNEW,Z0NEW,ROAIRN,WSTAR,
!    &                    USTOKES, VSTOKES, STRNMS)
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
!          *USTOKES*   U-COMP SURFACE STOKES DRIFT.
!          *VSTOKES*   V-COMP SURFACE STOKES DRIFT.
!          *STRNMS*    MEAN SQUARE STRAIN INTO THE SEA ICE (only if LWNEMOCOUSTRN).



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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWFLUX   ,LWVFLX_SNL, LWNEMOCOUSTRN
      USE YOWCOUT  , ONLY : LWFLUXOUT 
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWMEAN  , ONLY : PHIEPS   ,PHIAW    ,TAUOC
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "cimsstrn.intfb.h"
#include "femeanws.intfb.h"
#include "fkmean.intfb.h"
#include "frcutindex.intfb.h"
#include "sdissip.intfb.h"
#include "sinput.intfb.h"
#include "snonlin.intfb.h"
#include "stokesdrift.intfb.h"
#include "stresso.intfb.h"
#include "wnfluxes.intfb.h"
#include "wrong_wnfluxes.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS,IJL,IG
      INTEGER(KIND=JWIM), INTENT(OUT) :: MIJ(IJS:IJL)

      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: CICVR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: U10NEW, THWNEW, ROAIRN, WSTAR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: USNEW, Z0NEW
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: USTOKES, VSTOKES, STRNMS
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: FL3
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JWRB) :: TAU, XN, PHIDIAG, TAUO
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: TAUW_LOC  ! TAUW should not be updated do use a local array
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: EMEANALL, FMEANALL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: EMEANWS, FMEANWS
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: F1MEAN, AKMEAN, XKMEAN
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: PHIWA
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE) :: RHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: FL, SL, SPOS

!????
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: SSOURCE 

      LOGICAL :: LCFLX

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WDFLUXES',0,ZHOOK_HANDLE)

!*    1. INITIALISATION.
!        ---------------

      LCFLX=LWFLUX.OR.LWFLUXOUT
! ----------------------------------------------------------------------

!*    1.2 COMPUTATION OF RELEVANT SOURCE FUNCTIONS.
!         -----------------------------------------

      CALL SINPUT (FL3, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,             &
     &             ROAIRN, WSTAR, SL, SPOS, XLLWS)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. WDFLUXES: SINPUT CALLED'
        CALL FLUSH (IU06)
      ENDIF

!     MEAN FREQUENCY OF THE TOTAL SEA
      CALL FKMEAN(FL3, IJS, IJL, EMEANALL, FMEANALL,                    &
     &            F1MEAN, AKMEAN, XKMEAN)

!     MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
      CALL FEMEANWS(FL3,IJS,IJL,EMEANWS,FMEANWS,XLLWS)

!     COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
      CALL FRCUTINDEX(IJS, IJL, FMEANALL, FMEANWS,USNEW, CICVR,         &
     &                MIJ, RHOWGDFTH)


      IF(LCFLX) THEN
  
       DO M=1,NFRE
         DO K=1,NANG
           DO IJ=IJS,IJL
             SSOURCE(IJ,K,M) = SL(IJ,K,M)
           ENDDO
         ENDDO
       ENDDO

        CALL STRESSO (FL3, SPOS, IJS, IJL,                              &
     &                MIJ, RHOWGDFTH,                                   &
     &                THWNEW, USNEW, Z0NEW, ROAIRN,                     &
     &                TAUW_LOC, PHIWA)

        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '   SUB. WDFLUXES: STRESSO CALLED'
          CALL FLUSH (IU06)
        ENDIF

        CALL SDISSIP (FL3 ,FL, SL, IJS, IJL,                            &
     &                EMEANALL, F1MEAN, XKMEAN,                         &
     &                USNEW, THWNEW, ROAIRN)
        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '   SUB. WDFLUXES: SDISSIP CALLED'
          CALL FLUSH (IU06)
        ENDIF

        IF(.NOT. LWVFLX_SNL) THEN
          CALL WRONG_WNFLUXES (IJS, IJL,                                &
     &                         MIJ, RHOWGDFTH,                          &
     &                         SSOURCE, SL,                             &
     &                         PHIWA,                                   &
     &                         EMEANALL, F1MEAN, U10NEW, THWNEW,        &
     &                         USNEW, ROAIRN, .FALSE.)
        ENDIF

        CALL SNONLIN (FL3, FL, IJS, IJL, IG, SL, AKMEAN)
        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '   SUB. WDFLUXES: SNONLIN CALLED'
          CALL FLUSH (IU06)
        ENDIF

        IF(LWVFLX_SNL) THEN
          DO M=1,NFRE
            DO K=1,NANG
              DO IJ=IJS,IJL
!!!           SSOURCE should only contain the positive contribution from sinput
              SSOURCE(IJ,K,M) = SL(IJ,K,M)-SSOURCE(IJ,K,M)+SPOS(IJ,K,M)
              ENDDO
            ENDDO
          ENDDO

          CALL WNFLUXES (IJS, IJL,                                      &
     &                   MIJ, RHOWGDFTH,                                &
     &                   SSOURCE, CICVR,                                &
     &                   PHIWA,                                         &
     &                   EMEANALL, F1MEAN, U10NEW, THWNEW,              &
     &                   USNEW, ROAIRN, .FALSE.)
        ENDIF

      CALL STOKESDRIFT(FL3, IJS, IJL, USTOKES, VSTOKES)

      IF(LWNEMOCOUSTRN) CALL CIMSSTRN(FL3, IJS, IJL, STRNMS)


      ENDIF
! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WDFLUXES',1,ZHOOK_HANDLE)

      END SUBROUTINE WDFLUXES
