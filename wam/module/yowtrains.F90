      MODULE YOWTRAINS

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    SWELL TRAINS PARAMETERS.

      REAL(KIND=JWRB), ALLOCATABLE :: EMTRAIN(:,:)
      REAL(KIND=JWRB), ALLOCATABLE :: THTRAIN(:,:)
      REAL(KIND=JWRB), ALLOCATABLE :: PMTRAIN(:,:)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *ETRAIN*   REAL      TOTAL ENERGY.
!      *THTRAIN*  REAL      MEAN DIRECTION.
!      *P1TRAIN*  REAL      MEAN PERIOD (-1).
! ----------------------------------------------------------------------
      END MODULE YOWTRAINS
