FUNCTION DIRSPRD_GC (U10)

! ----------------------------------------------------------------------

!**** *DIRSPRD_GC* - FUNCTION TO SPECIFY THE ANGULAR SPREADING FACTOR IN
!                    THE GRAVITY-CAPILLARY MODEL

!**   INTERFACE.
!     ----------

!       *FUNCTION* *DIRSPRD_GC (U10)*
!          *U10*     - 10m WIND SPEED

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWPHYS  , ONLY : ANG_GC_A, ANG_GC_B, ANG_GC_C, ANG_GC_D, ANG_GC_E, ANG_GC_N
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB) :: DIRSPRD_GC

      REAL(KIND=JWRB), INTENT(IN) :: U10

      REAL(KIND=JWRB) :: ANG_GC_10, C10 
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('DIRSPRD_GC',0,ZHOOK_HANDLE)

IF(U10 > 10.0_JWRB) THEN
  DIRSPRD_GC= ANG_GC_A+ANG_GC_B*TANH(ANG_GC_C*(U10-ANG_GC_D))
ELSE
  ANG_GC_10= ANG_GC_A+ANG_GC_B*TANH(ANG_GC_C*(10.0_JWRB-ANG_GC_D))
  C10=ANG_GC_10*0.1_JWRB**ANG_GC_N
  DIRSPRD_GC= MAX(C10*U10**ANG_GC_N,ANG_GC_E)
ENDIF

IF (LHOOK) CALL DR_HOOK('DIRSPRD_GC',1,ZHOOK_HANDLE)

END FUNCTION DIRSPRD_GC
