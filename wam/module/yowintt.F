      MODULE YOWINTT

      IMPLICIT NONE

!*    ** *INTSPAR*  GRIDDED INTEGRATED TRAIN PARAMETERS.

      REAL, ALLOCATABLE :: WHGTTRG(:,:,:) 
      REAL, ALLOCATABLE :: WDIRTRG(:,:,:) 
      REAL, ALLOCATABLE :: MWPMTRG(:,:,:) 
!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *WHGTTRG*    REAL     GRIDDED MEAN TRAIN ENERGY (M*M) CHANGED TO
!                            WAVE HEIGHT (M) IN SUB MAKEGRID.
!      *WDIRTRG*    REAL     GRIDDED MEAN TRAIN DIRECTION (RAD) CHANGED
!                            TO DEGREE AND METEOROLOGICAL CONVENTION
!                            IN SUB OUTINT.
!      *MWPMTRG*    REAL     GRIDDED MEAN TRAIN WAVE PERIOD -1 (S.). 

! ----------------------------------------------------------------------
      END MODULE YOWINTT
