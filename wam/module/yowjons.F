      MODULE YOWJONS

      IMPLICIT NONE

!*    **  *JONS* - JONSWAP PARAMETERS.

      REAL, ALLOCATABLE :: FP(:) 
      REAL, ALLOCATABLE :: ALPHJ(:) 
      REAL, ALLOCATABLE :: THES(:) 
      REAL FM 
      REAL ALFA 
      REAL GAMMA 
      REAL SA 
      REAL SB 
      REAL THETAQ 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
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
