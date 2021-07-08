      FUNCTION CHNKMIN (U10)

! ----------------------------------------------------------------------

!**** *CHNKMIN* - FUNCTION TO COMPUTE THE MINMUM CHARNOCK

!*    PURPOSE.
!     -------


!**   INTERFACE.
!     ----------

!       *FUNCTION* *CHNKMIN (U10)*

!     METHOD.
!     -------

!     CHNKMIN = ALPHAMIN + (ALPHA-ALPHAMIN)*0.5_JWRB*(1.0_JWRB-TANH(U10-A)) 

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPHYS  , ONLY : ALPHA, ALPHAMIN, CHNKMIN_U 
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB) :: CHNKMIN

      REAL(KIND=JWRB), INTENT(IN) :: U10
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CHNKMIN',0,ZHOOK_HANDLE)

      CHNKMIN = ALPHAMIN + (ALPHA-ALPHAMIN)*0.5_JWRB*(1.0_JWRB-TANH(U10-CHNKMIN_U))

      IF (LHOOK) CALL DR_HOOK('CHNKMIN',1,ZHOOK_HANDLE)

      END FUNCTION CHNKMIN
