      SUBROUTINE MEANSQS(XKMSS, IJS, IJL, KIJS, KIJL, F, U10, USTAR, UDIR, XMSS)

! ----------------------------------------------------------------------

!**** *MEANSQS* - COMPUTATION OF MEAN SQUARE SLOPE UP TO WAVE NUMBER XKMSS.

!*    PURPOSE.
!     --------

!       COMPUTE MEAN SQUARE SLOPE AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *MEANSQS (XKMSS, IJS, IJL, KIJS, KIJL, F, U10, USTAR, UDIR, XMSS)*
!              *XKMSS*   - WAVE NUMBER CUT OFF
!              *IJS:IJL* - 1st DIMENSION of F
!              *KIJS*    - INDEX OF FIRST GRIDPOINT
!              *KIJL*    - INDEX OF LAST GRIDPOINT
!              *F*       - SPECTRUM.
!              *U10*     - 10m wind speed
!              *USTAR*   - NEW FRICTION VELOCITY IN M/S (INPUT).
!              *UDIR*    - WIND SPEED DIRECTION
!              *XMSS*    - MEAN SQUARE SLOPE (OUTPUT).

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : G        ,GM1      ,ZPI
      USE YOWFRED  , ONLY : FR       ,ZPIFR    ,FRATIO
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "halphap.intfb.h"
#include "meansqs_gc.intfb.h"
#include "meansqs_lf.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F

      REAL(KIND=JWRB), INTENT(IN) :: XKMSS 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: U10
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: USTAR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: UDIR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: XMSS 

      INTEGER(KIND=JWIM) :: IJ, NFRE_MSS, NFRE_EFF

      REAL(KIND=JWRB) :: XLOGFS, FCUT
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE) :: FD
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TEMP1, TEMP2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: XMSSLF, XMSS_TAIL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: HALP, FRGC

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('MEANSQS',0,ZHOOK_HANDLE)

!*    1. COMPUTE THE GRAVITY-CAPILLARY CONTRIBUTION TO MSS 
!        -------------------------------------------------

!     COMPUTE THE PHILLIPS PARAMETER
      CALL HALPHAP(IJS, IJL, KIJS, KIJL, USTAR, UDIR, F, HALP)

!     GRAVITY-CAPILLARY CONTRIBUTION TO MSS
      CALL MEANSQS_GC(XKMSS, KIJS, KIJL, HALP, U10, USTAR, XMSS, FRGC)

!     MODEL SPECTRUM PART
      FCUT = SQRT(G*XKMSS)/ZPI
      NFRE_MSS = INT(LOG(FCUT/FR(1))/LOG(FRATIO))+1 
      NFRE_EFF = MIN(NFRE,NFRE_MSS)

      CALL MEANSQS_LF (NFRE_EFF, IJS, IJL, KIJS, KIJL, F, XMSSLF)
      XMSS(:) = XMSS(:) + XMSSLF(:)

!     ADD TAIL CORRECTION TO MEAN SQUARE SLOPE (between FR(NFRE_EFF) and FRGC).

      XLOGFS = LOG(FR(NFRE_EFF))
      DO IJ=KIJS,KIJL
        XMSS_TAIL(IJ) = 2.0_JWRB*HALP(IJ)*MAX((LOG(MIN(FRGC(IJ),FCUT)) - XLOGFS),0.0_JWRB)
        XMSS(IJ) = XMSS(IJ) + XMSS_TAIL(IJ)
      ENDDO

      IF (LHOOK) CALL DR_HOOK('MEANSQS',1,ZHOOK_HANDLE)

      END SUBROUTINE MEANSQS
