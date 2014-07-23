      MODULE YOWTRAINS

      IMPLICIT NONE

!*    SWELL TRAINS PARAMETERS.

      REAL, ALLOCATABLE :: EMTRAIN(:,:)
      REAL, ALLOCATABLE :: THTRAIN(:,:)
      REAL, ALLOCATABLE :: PMTRAIN(:,:)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *ETRAIN*   REAL      TOTAL ENERGY.
!      *THTRAIN*  REAL      MEAN DIRECTION.
!      *P1TRAIN*  REAL      MEAN PERIOD (-1).
! ----------------------------------------------------------------------
      END MODULE YOWTRAINS
