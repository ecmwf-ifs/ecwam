      MODULE YOWREFD

      IMPLICIT NONE

!*    ** *REFDOT* - DEPTH AND CURRENT PART OF THETA DOT.

      REAL,ALLOCATABLE :: THDD(:,:) 
      REAL,ALLOCATABLE :: THDC(:,:) 
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SDOT

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *THDD*      REAL      DEPTH GRADIENT PART OF THETA DOT.
!      *THDC*      REAL      CURRENT GRADIENT PART OF THETA DOT.
!      *SDOT*      SIGMA DOT ARRAY

! ----------------------------------------------------------------------
      END MODULE YOWREFD
