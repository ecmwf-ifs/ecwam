      SUBROUTINE SDISSIP (GFL, FLD, SL, IJS, IJL, KIJS, KIJL,   &
     &                    EMEAN, F1MEAN, XKMEAN,                &
     &                    USNEW, THWNEW, ROAIRN)
! ----------------------------------------------------------------------

!**** *SDISSIP* - COMPUTATION OF DEEP WATER DISSIPATION SOURCE FUNCTION.


!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP (GFL, FLD, SL, IJS, IJL, KIJS, KIJL, *
!                        EMEAN, F1MEAN, XKMEAN,*
!                        USNEW, THWNEW, ROAIRN)*
!        *GFL*   - SPECTRUM.
!         *FLD*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY
!        *IJS:IJL- 1st DIMENSION OF GFL
!         *KIJS* - INDEX OF FIRST GRIDPOINT
!         *KIJL* - INDEX OF LAST GRIDPOINT
!        *EMEAN* - MEAN ENERGY DENSITY 
!       *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
!       *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.
!       *USNEW*  - NEW FRICTION VELOCITY IN M/S.
!       *ROAIRN* - AIR DENSITY IN KG/M3
!       *THWNEW* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC.

!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       ARDHUIN et AL. JPO DOI:10.1175/20110JPO4324.1


! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : IPHYS

      USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "sdissip_ard.intfb.h"
#include "sdissip_jan.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: GFL 

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: EMEAN, F1MEAN, XKMEAN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: USNEW, THWNEW, ROAIRN 

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SDISSIP',0,ZHOOK_HANDLE)


      SELECT CASE (IPHYS)
      CASE(0)
         CALL SDISSIP_JAN (GFL ,FLD, SL, IJS, IJL, KIJS, KIJL,  &
     &                     EMEAN, F1MEAN, XKMEAN)

      CASE(1) 
         CALL SDISSIP_ARD (GFL ,FLD, SL, IJS, IJL, KIJS, KIJL,  &
     &                     USNEW, THWNEW, ROAIRN)
      END SELECT 

      IF (LHOOK) CALL DR_HOOK('SDISSIP',1,ZHOOK_HANDLE)

      END SUBROUTINE SDISSIP
