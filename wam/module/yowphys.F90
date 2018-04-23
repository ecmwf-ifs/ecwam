      MODULE YOWPHYS

!*    ** *YOWPHYS* - PARAMETERS FOR WAVE PHYSICS PARAMETERISATION 

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!     MAXIMUM PHILLIPS PARAMETER USED TO CONTROL MAXIMUM STEEPNESS
      REAL(KIND=JWRB) :: ALPHAPMAX

!     ARDHUIN et al. 2010:
!     TEST 473:
!     Br:
      REAL(KIND=JWRB), PARAMETER :: SDSBR = 9.0E-4_JWRB

!     Saturation dissipation coefficient
      INTEGER(KIND=JWIM), PARAMETER :: ISDSDTH = 80_JWIM 
      INTEGER(KIND=JWIM), PARAMETER :: ISB=2_JWIM
      INTEGER(KIND=JWIM), PARAMETER :: IPSAT=2_JWIM
      REAL(KIND=JWRB), PARAMETER :: SSDSC2 = -2.2E-5_JWRB 
      REAL(KIND=JWRB), PARAMETER :: SSDSC4 = 1.0_JWRB
      REAL(KIND=JWRB), PARAMETER :: SSDSC6 = 0.3_JWRB 
      REAL(KIND=JWRB), PARAMETER :: MICHE = 1.0_JWRB 


!     Cumulative dissipation coefficient
!!!      REAL(KIND=JWRB), PARAMETER :: SSDSC3 = -0.40344_JWRB 
!!!debile
      REAL(KIND=JWRB), PARAMETER :: SSDSC3 = 0.0_JWRB 
      REAL(KIND=JWRB), PARAMETER :: SSDSBRF1   = 0.5_JWRB
!     28.16 = 22.0 * 1.6² * 1/2 with  
!     22.0 (Banner & al. 2000, figure 6) 
!     1.6  the coefficient that transforms  SQRT(B) to Banner et al. (2000)'s epsilon
!     1/2  factor to correct overestimation of Banner et al. (2000)'s breaking probability due to zero-crossing analysis
      REAL(KIND=JWRB), PARAMETER :: BRKPBCOEF=28.16_JWRB  

!     Wave-turbulence interaction coefficient 
      REAL(KIND=JWRB), PARAMETER :: SSDSC5  = 0.0_JWRB 

!     NSDSNTH is the number of directions on both used to compute the spectral saturation  
      INTEGER(KIND=JWIM) :: NSDSNTH
!     NDIKCUMUL is the  integer difference in frequency bands
      INTEGER(KIND=JWIM) :: NDIKCUMUL

      INTEGER(KIND=JWIM), ALLOCATABLE :: INDICESSAT(:,:)
      REAL(KIND=JWRB), ALLOCATABLE :: SATWEIGHTS(:,:)
      REAL(KIND=JWRB), ALLOCATABLE :: CUMULW(:,:,:,:)
! ----------------------------------------------------------------------
      END MODULE YOWPHYS
