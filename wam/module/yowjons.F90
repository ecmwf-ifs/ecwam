      MODULE YOWJONS

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    **  *JONS* - JONSWAP PARAMETERS.

      REAL(KIND=JWRB), PARAMETER :: AJONS = 2.84_JWRB
      REAL(KIND=JWRB), PARAMETER :: BJONS = 0.033_JWRB
      REAL(KIND=JWRB), PARAMETER :: DJONS = -3.0_JWRB/10.0_JWRB
      REAL(KIND=JWRB), PARAMETER :: EJONS = 2.0_JWRB/3.0_JWRB
      REAL(KIND=JWRB), ALLOCATABLE :: FP(:) 
      REAL(KIND=JWRB), ALLOCATABLE :: ALPHJ(:) 
      REAL(KIND=JWRB), ALLOCATABLE :: THES(:) 
      REAL(KIND=JWRB)              :: FM 
      REAL(KIND=JWRB)              :: ALFA 
      REAL(KIND=JWRB)              :: GAMMA 
      REAL(KIND=JWRB)              :: SA 
      REAL(KIND=JWRB)              :: SB 
      REAL(KIND=JWRB)              :: THETAQ 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *AJONS*     REAL      CONSTANT USED IN FETCH LIMITED GROWTH
!                            SEE ROUTINE PEAK
!      *BJONS*     REAL      CONSTANT USED IN FETCH LIMITED GROWTH
!                            SEE ROUTINE PEAK
!      *DJONS*     REAL      CONSTANT USED IN FETCH LIMITED GROWTH
!                            SEE ROUTINE PEAK
!      *EJONS*     REAL      CONSTANT USED IN FETCH LIMITED GROWTH
!                            SEE ROUTINE PEAK
!      *FP*        REAL      PEAK FREQUENCY OF SPECTRA IN A BLOCK (HZ).
!      *ALPHJ*     REAL      ALPHA PARAMETER OF SPECTRA IN A BLOCK.
!      *THES*      REAL      MEAN DIRECTION OF SPECTRA IN A BLOCK (RAD).
!      *FM*        REAL      PEAK FREQUENCY AS DEFINED BY INPUT (HZ).
!      *ALFA*      REAL      ALPHA PARAMETER AS DEFINED BY INPUT.
!      *GAMMA*     REAL      OVERSHOOT FACTOR.
!      *SA*        REAL      LEFT PEAK WIDTH.
!      *SB*        REAL      RIGHT PEAK WIDTH.
!      *THETAQ*    REAL      MEAN DIRECTION AS DEFINED BY INPUT (RAD).

! ----------------------------------------------------------------------
      END MODULE YOWJONS
