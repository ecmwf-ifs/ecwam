      MODULE YOWMEAN

      IMPLICIT NONE

!*    **  *MEANPA* - INTEGRATED PARAMETERS OF A BLOCK.

      REAL,ALLOCATABLE :: EMEAN(:) 
      REAL,ALLOCATABLE :: FMEAN(:) 
      REAL,ALLOCATABLE :: THQ(:) 
      REAL,ALLOCATABLE :: AKMEAN(:) 
      REAL,ALLOCATABLE :: SMEAN(:) 
      REAL,ALLOCATABLE :: P1MEAN(:) 
      REAL,ALLOCATABLE :: P2MEAN(:) 
      REAL,ALLOCATABLE :: SPRDMEAN(:) 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *EMEAN*     REAL      TOTAL ENERGY.
!      *FMEAN*     REAL      MEAN FREQUENCY (1./MEAN PERIOD).
!      *THQ*       REAL      MEAN WAVE DIRECTION IN RADIANS.
!      *AKMEAN*    REAL      MEAN WAVE NUMBER.
!      *SMEAN*     REAL      MEAN SQUARE SLOPE.
!      *P1MEAN*    REAL      MEAN PERIOD (1).
!      *P2MEAN*    REAL      MEAN PERIOD (2).
!      *SPRDMEAN*  REAL      MEAN DIRECTIONAl SPREAD.

! ----------------------------------------------------------------------
      END MODULE YOWMEAN
