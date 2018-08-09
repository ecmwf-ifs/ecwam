      SUBROUTINE SINPUT (F, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,&
     &                   ROAIRN, WSTAR, SL, XLLWS)
! ----------------------------------------------------------------------

!**** *SINPUT* - COMPUTATION OF INPUT SOURCE FUNCTION.


!**   INTERFACE.
!     ----------

!     *CALL* *SINPUT (F, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,
!    &                   ROAIRN,WSTAR, SL, XLLWS)
!            *F* - SPECTRUM.
!           *FL* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!          *IJS* - INDEX OF FIRST GRIDPOINT.
!          *IJL* - INDEX OF LAST GRIDPOINT.
!       *THWNEW* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!        *USNEW* - NEW FRICTION VELOCITY IN M/S.
!        *Z0NEW* - ROUGHNESS LENGTH IN M.
!       *ROAIRN* - AIR DENSITY IN KG/M3
!        *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
!           *SL* - TOTAL SOURCE FUNCTION ARRAY.
!         *XLLWS*- 1 WHERE SINPUT IS POSITIVE

!     METHOD.
!     -------

!     DEPENDING ON THE VALUE OF IPHYS, DIFFERENT INPUTE SOURCE TERm WILL BE CALLED

!     EXTERNALS.
!     ----------

!     MODIFICATIONS
!     -------------

!     REFERENCE.
!     ----------


! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : IPHYS
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "sinput_ard.intfb.h"
#include "sinput_jan.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS,IJL

      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: THWNEW, USNEW, Z0NEW, ROAIRN, WSTAR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL,SL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(OUT) :: XLLWS

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SINPUT',0,ZHOOK_HANDLE)

      SELECT CASE (IPHYS)
      CASE(0)
        CALL SINPUT_JAN (F, FL, IJS, IJL, THWNEW, USNEW, Z0NEW, &
     &                   ROAIRN, WSTAR, SL, XLLWS)
      CASE(1) 
        CALL SINPUT_ARD (F, FL, IJS, IJL, THWNEW, USNEW, Z0NEW, &
     &                   ROAIRN, WSTAR, SL, XLLWS)
      END SELECT 

      IF (LHOOK) CALL DR_HOOK('SINPUT',1,ZHOOK_HANDLE)

      END SUBROUTINE SINPUT
