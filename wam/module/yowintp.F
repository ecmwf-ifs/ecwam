      MODULE YOWINTP

      IMPLICIT NONE

!*    ** *INTPAR*  GRIDDED INTEGRATED PARAMETER TOTAL SPECTRUM.

      REAL, ALLOCATABLE :: WHGTTG(:,:) 
      REAL, ALLOCATABLE :: WDIRTG(:,:) 
      REAL, ALLOCATABLE :: WPKFTG(:,:) 
      REAL, ALLOCATABLE :: WMNFTG(:,:) 
      REAL, ALLOCATABLE :: USTARG(:,:) 
      REAL, ALLOCATABLE :: UDIRG(:,:) 
      REAL, ALLOCATABLE :: TAUWG(:,:) 
      REAL, ALLOCATABLE :: CDG(:,:) 
      REAL, ALLOCATABLE :: SMEANG(:,:) 
      REAL, ALLOCATABLE :: U10G(:,:) 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *WHGTTG*    REAL      GRIDDED MEAN ENERGY (M*M) CHANGED TO
!                            WAVE HEIGHT (M) IN SUB OUTINT.
!      *WDIRTG*    REAL      GRIDDED MEAN DIRECTION (RAD) CHANGED TO
!                            DEGREE IN SUB OUTINT.
!      *WPKFTG*    REAL      GRIDDED PEAK FREQUENCY (HZ) until gribbing
!                            where it is transformed into PERIOD.
!      *WMNFTG*    REAL      GRIDDED MEAN FREQUENCY (HZ) until gribbing
!                            where it is transformed into PERIOD.
!      *USTARG*    REAL      GRIDDED FRICTION VELOCITY (M/S).
!      *UDIRG*     REAL      GRIDDED WIND DIRECTION (RAD) CHANGED TO
!                            DEGREE IN SUB OUTINT.
!      *TAUWG*     REAL      GRIDDED NORMALISED WAVE STRESS.
!      *CDG*       REAL      GRIDDED DRAG COEFFICENT.
!      *SMEANG*    REAL      GRIDDED MEAN SQUARE SLOPE.
!      *U10G*      REAL      GRIDDED WIND SPEED (M/S).

! ----------------------------------------------------------------------
      END MODULE YOWINTP
