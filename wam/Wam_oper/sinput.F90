      SUBROUTINE SINPUT (NGST, KIJS, KIJL, FL1,           & 
     &                   WAVNUM, CINV, CGROUP,            &
     &                   THWNEW, U10NEW, USNEW, Z0NEW,    &
     &                   ROAIRN, WSTAR, RNFAC,            &
     &                   FLD, SL, SPOS, XLLWS)
! ----------------------------------------------------------------------

!**** *SINPUT* - COMPUTATION OF INPUT SOURCE FUNCTION.


!**   INTERFACE.
!     ----------

!     *CALL* *SINPUT (NGST, KIJS, KIJL, FL1,
!    &                WAVNUM, CINV, CGROUP,
!    &                THWNEW, USNEW, Z0NEW,
!    &                ROAIRN, WSTAR, FLD, SL, SPOS, XLLWS)
!         *NGST* - IF = 1 THEN NO GUSTINESS PARAMETERISATION
!                - IF = 2 THEN GUSTINESS PARAMETERISATION
!         *KIJS* - INDEX OF FIRST GRIDPOINT.
!         *KIJL* - INDEX OF LAST GRIDPOINT.
!          *FL1* - SPECTRUM.
!       *WAVNUM* - WAVE NUMBER.
!         *CINV* - INVERSE PHASE VELOCITY.
!       *CGROUP* - GROUP SPPED.
!       *THWNEW* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!        *USNEW* - NEW FRICTION VELOCITY IN M/S.
!        *Z0NEW* - ROUGHNESS LENGTH IN M.
!       *ROAIRN* - AIR DENSITY IN KG/M3
!        *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
!        *RNFAC* - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
!          *FLD* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!           *SL* - TOTAL SOURCE FUNCTION ARRAY.
!         *SPOS* - POSITIVE SOURCE FUNCTION ARRAY.
!        *XLLWS* - = 1 WHERE SINPUT IS POSITIVE

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

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "sinput_ard.intfb.h"
#include "sinput_jan.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: NGST
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM, CINV, CGROUP

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: THWNEW, U10NEW, USNEW, Z0NEW
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: ROAIRN, WSTAR, RNFAC

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: FLD, SL, SPOS

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: XLLWS

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SINPUT',0,ZHOOK_HANDLE)

      SELECT CASE (IPHYS)
      CASE(0)
        CALL SINPUT_JAN (NGST, KIJS, KIJL, FL1,           &
     &                   WAVNUM, CINV, CGROUP,            &
     &                   THWNEW, U10NEW, USNEW, Z0NEW,    &
     &                   ROAIRN, WSTAR, RNFAC,            &
     &                   FLD, SL, SPOS, XLLWS)
      CASE(1) 
        CALL SINPUT_ARD (NGST, KIJS, KIJL, FL1,           &
     &                   WAVNUM, CINV, CGROUP,            &
     &                   THWNEW, U10NEW, USNEW, Z0NEW,    &
     &                   ROAIRN, WSTAR, RNFAC,            &
     &                   FLD, SL, SPOS, XLLWS)
      END SELECT 

      IF (LHOOK) CALL DR_HOOK('SINPUT',1,ZHOOK_HANDLE)

      END SUBROUTINE SINPUT
