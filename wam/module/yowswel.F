      MODULE YOWSWEL

      IMPLICIT NONE

!*    ** *SWELL* - SEA AND SWELL PARAMETERS.

      REAL, ALLOCATABLE :: ESWELL(:) 
      REAL, ALLOCATABLE :: FSWELL(:) 
      REAL, ALLOCATABLE :: THSWELL(:) 
      REAL, ALLOCATABLE :: ESEA(:) 
      REAL, ALLOCATABLE :: FSEA(:) 
      REAL, ALLOCATABLE :: THWISEA(:) 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *ESWELL*    REAL      TOTAL SWELL ENERGY
!      *FSWELL*    REAL      MEAN SWELL FREQUENCY
!      *THSWELL*   REAL      MEAN SWELL DIRECTION
!      *ESEA*      REAL      TOTAL SEA ENERGY
!      *FSEA*      REAL      MEAN SEA FREQUENCY
!      *THWISEA*   REAL      MEAN WIND INTERACTING SEA DIRECTION

! ----------------------------------------------------------------------
      END MODULE YOWSWEL
