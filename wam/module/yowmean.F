      MODULE YOWMEAN

      IMPLICIT NONE

!*    **  *MEANPA* - INTEGRATED PARAMETERS OF A BLOCK.

      REAL,ALLOCATABLE :: EMEAN(:) 
      REAL,ALLOCATABLE :: FMEAN(:) 
      REAL,ALLOCATABLE :: THQ(:) 
      REAL,ALLOCATABLE :: AKMEAN(:) 
      REAL,ALLOCATABLE :: SMEAN(:) 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *EMEAN*     REAL      TOTAL ENERGY.
!      *FMEAN*     REAL      MEAN FREQUENCY (1./MEAN PERIOD).
!      *THQ*       REAL      MEAN WAVE DIRECTION IN RADIANS.
!      *AKMEAN*    REAL      MEAN WAVE NUMBER.
!      *SMEAN*     REAL      MEAN SQUARE SLOPE.

! ----------------------------------------------------------------------
      END MODULE YOWMEAN
