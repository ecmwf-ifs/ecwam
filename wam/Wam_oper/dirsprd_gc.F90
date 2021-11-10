REAL(KIND=JWRB) FUNCTION DIRSPRD_GC (USTAR)

! ----------------------------------------------------------------------

!**** *DIRSPRD_GC* - FUNCTION TO SPECIFY THE ANGULAR SPREADING FACTOR IN
!                    THE GRAVITY-CAPILLARY MODEL

!**   INTERFACE.
!     ----------

!       *FUNCTION* *DIRSPRD_GC (USTAR)*
!          *USTAR*  - FRICTION VELOCITY. 

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWPHYS  , ONLY : ANG_GC_A, ANG_GC_B, ANG_GC_C
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN) :: USTAR

      REAL(KIND=JWRB) :: CC 
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('DIRSPRD_GC',0,ZHOOK_HANDLE)

DIRSPRD_GC = ANG_GC_A + ANG_GC_B * TANH(ANG_GC_C * USTAR**2)

IF (LHOOK) CALL DR_HOOK('DIRSPRD_GC',1,ZHOOK_HANDLE)

END FUNCTION DIRSPRD_GC
