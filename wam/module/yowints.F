      MODULE YOWINTS

      IMPLICIT NONE

!*    ** *INTSPAR*  GRIDDED INTEGRATED SWELL PARAMETERS.

      REAL, ALLOCATABLE :: WHGTSG(:,:) 
      REAL, ALLOCATABLE :: WDIRSG(:,:) 
      REAL, ALLOCATABLE :: WMNFSG(:,:) 
      REAL, ALLOCATABLE :: WHGTWG(:,:) 
      REAL, ALLOCATABLE :: WDIRWG(:,:) 
      REAL, ALLOCATABLE :: WMNFWG(:,:) 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *WHGTSG*    REAL      GRIDDED MEAN SWELL ENERGY (M*M) CHANGED TO
!                            WAVE HEIGHT (M) IN SUB OUTINT.
!      *WDIRSG*    REAL      GRIDDED MEAN SWELL DIRECTION (RAD) CHANGED
!                            TO DEGREE IN SUB OUTINT.
!      *WMNFSG*    REAL      GRIDDED MEAN SWELL FREQUENCY (HZ).
!      *WHGTWG*    REAL      GRIDDED WIND SEA ENERGY (M*M) CHANGED TO
!                            WAVE HEIGHT (M) IN SUB OUTINT. 
!      *WDIRWG*    REAL      GRIDDED MEAN SEA   DIRECTION (RAD) CHANGED
!                            TO DEGREE IN SUB OUTINT.
!      *WMNFWG*    REAL      GRIDDED MEAN WIND SEA FREQUENCY (HZ).

! ----------------------------------------------------------------------
      END MODULE YOWINTS
