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

      USE YOWCOUP  , ONLY : ALPHA 
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB) :: CHNKMIN

      REAL(KIND=JWRB), PARAMETER :: A = 33.0_JWRB
      REAL(KIND=JWRB), PARAMETER :: ALPHAMIN=0.0001_JWRB
      ! Parameter for the linear correction of Charnock values for winds below 5m/s
      ! i.e Minimum Charnock: ALPHA0 for U10=0, ramping down to ALPHA for U10=5, then ALPHA for U10>5m/s
      REAL(KIND=JWRB), PARAMETER :: ALPHA0=0.015_JWRB


      REAL(KIND=JWRB), INTENT(IN) :: U10
      REAL(KIND=JWRB) :: ALPHAT 
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CHNKMIN',0,ZHOOK_HANDLE)

      ! ALPHA is adjusted for winds below 5m/s  (0.2=1/5)
      ALPHAT = MAX(ALPHA, ALPHA0 - (ALPHA0-ALPHA)*U10*0.2_JWRB )
      CHNKMIN = ALPHAMIN + (ALPHAT-ALPHAMIN)*0.5_JWRB*(1.0_JWRB-TANH(U10-A))

      IF (LHOOK) CALL DR_HOOK('CHNKMIN',1,ZHOOK_HANDLE)

      END FUNCTION CHNKMIN
