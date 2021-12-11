      SUBROUTINE SDISSIP (KIJS, KIJL, FL1, FLD, SL,  &
     &                    INDEP, WAVNUM, XK2CG,      &
     &                    EMEAN, F1MEAN, XKMEAN,     &
     &                    UFRIC, WDWAVE, AIRD)
! ----------------------------------------------------------------------

!**** *SDISSIP* - COMPUTATION OF DEEP WATER DISSIPATION SOURCE FUNCTION.


!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP (KIJS, KIJL, FL1, FLD, SL, *
!                        INDEP, WAVNUM, XK2CG,  
!                        EMEAN, F1MEAN, XKMEAN,*
!                        UFRIC, WDWAVE, AIRD)*
!         *KIJS* - INDEX OF FIRST GRIDPOINT
!         *KIJL* - INDEX OF LAST GRIDPOINT
!         *FL1*  - SPECTRUM.
!         *FLD*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY
!       *INDEP*  - DEPTH INDEX
!       *WAVNUM* - WAVE NUMBER
!       *XK2CG*  - (WAVNUM)**2 * GROUP SPEED
!        *EMEAN* - MEAN ENERGY DENSITY 
!       *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
!       *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.
!       *UFRIC*  - FRICTION VELOCITY IN M/S.
!       *AIRD*   - AIR DENSITY IN KG/M3
!       *WDWAVE* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC.


! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : IPHYS

      USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "sdissip_ard.intfb.h"
#include "sdissip_jan.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: INDEP
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM, XK2CG
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: EMEAN, F1MEAN, XKMEAN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: UFRIC, WDWAVE, AIRD

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SDISSIP',0,ZHOOK_HANDLE)

      SELECT CASE (IPHYS)
      CASE(0)
         CALL SDISSIP_JAN (KIJS, KIJL, FL1 ,FLD, SL,  &
     &                     WAVNUM,                    &
     &                     EMEAN, F1MEAN, XKMEAN)

      CASE(1) 
         CALL SDISSIP_ARD (KIJS, KIJL, FL1 ,FLD, SL,   &
     &                     INDEP, WAVNUM, XK2CG,       &
     &                     UFRIC, WDWAVE, AIRD)
      END SELECT 

      IF (LHOOK) CALL DR_HOOK('SDISSIP',1,ZHOOK_HANDLE)

      END SUBROUTINE SDISSIP
