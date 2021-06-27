FUNCTION DIRSPRD_GC (USTAR)

! ----------------------------------------------------------------------

!**** *DIRSPRD_GC* - FUNCTION TO SPECIFY THE ANGULAR SPREADING FACTOR IN
!                    THE GRAVITY-CAPILLARY MODEL

!**   INTERFACE.
!     ----------

!       *FUNCTION* *DIRSPRD_GC (USTAR)*
!          *USTAR*  - FRICTION VELOCITY. 

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWPHYS  , ONLY : ANG_GC_A, ANG_GC_B, ANG_GC_C, ANG_GC_D, ANG_GC_E, ANG_GC_N, ANG_GC_U
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB) :: DIRSPRD_GC

      REAL(KIND=JWRB), INTENT(IN) :: USTAR

      REAL(KIND=JWRB) :: CC 
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

!     INLINE FUNCTION.
!     ----------------

      REAL(KIND=JWRB) :: ANG_GC, X 

      ANG_GC(X)= ANG_GC_A+ANG_GC_B*TANH(ANG_GC_C*(X-ANG_GC_D))

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('DIRSPRD_GC',0,ZHOOK_HANDLE)

IF(USTAR > ANG_GC_U) THEN
  DIRSPRD_GC= ANG_GC(USTAR)
ELSE
  CC = ANG_GC(ANG_GC_U)/ANG_GC_U**ANG_GC_N
  DIRSPRD_GC= MAX(CC*USTAR**ANG_GC_N,ANG_GC_E)
ENDIF

!!!1debile
!!1DIRSPRD_GC = 0.1_JWRB + 1.15_JWRB * TANH(3.7_JWRB * USTAR**2)
DIRSPRD_GC = 0.3_JWRB + 0.5_JWRB * TANH(4.0_JWRB * USTAR**4)


IF (LHOOK) CALL DR_HOOK('DIRSPRD_GC',1,ZHOOK_HANDLE)

END FUNCTION DIRSPRD_GC
