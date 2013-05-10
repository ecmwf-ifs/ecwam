      MODULE YOWINTT

      IMPLICIT NONE

!*    ** *INTSPAR*  GRIDDED INTEGRATED TRAIN PARAMETERS.

      REAL, ALLOCATABLE :: WHGTT1G(:,:) 
      REAL, ALLOCATABLE :: WDIRT1G(:,:) 
      REAL, ALLOCATABLE :: MWP1T1G(:,:) 
      REAL, ALLOCATABLE :: WHGTT2G(:,:) 
      REAL, ALLOCATABLE :: WDIRT2G(:,:) 
      REAL, ALLOCATABLE :: MWP1T2G(:,:) 
      REAL, ALLOCATABLE :: WHGTT3G(:,:) 
      REAL, ALLOCATABLE :: WDIRT3G(:,:) 
      REAL, ALLOCATABLE :: MWP1T3G(:,:) 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *WHGTT1G*    REAL     GRIDDED MEAN TRAIN 1 ENERGY (M*M) CHANGED TO
!                            WAVE HEIGHT (M) IN SUB MAKEGRID.
!      *WDIRT1G*    REAL     GRIDDED MEAN TRAIN 1 DIRECTION (RAD) CHANGED
!                            TO DEGREE AND METEOROLOGICAL CONVENTION
!                            IN SUB OUTINT.
!      *MWP1T1G*    REAL     GRIDDED MEAN TRAIN 1 WAVE PERIOD 1 (S.). 
!      *WHGTT2G*    REAL     GRIDDED TRAIN 2 ENERGY (M*M) CHANGED TO
!                            WAVE HEIGHT (M) IN SUB MAKEGRID. 
!      *WDIRT2G*    REAL     GRIDDED TRAIN 2 DIRECTION (RAD) CHANGED
!                            TO DEGREE AND METEOROLOGICAL CONVENTION
!                            IN SUB OUTINT.
!      *MWP1T2G*    REAL     GRIDDED MEAN TRAIN 2 WAVE PERIOD 1 (S.). 
!      *WHGTT3G*    REAL     GRIDDED TRAIN 3 ENERGY (M*M) CHANGED TO
!                            WAVE HEIGHT (M) IN SUB MAKEGRID. 
!      *WDIRT3G*    REAL     GRIDDED TRAIN 3 DIRECTION (RAD) CHANGED
!                            TO DEGREE AND METEOROLOGICAL CONVENTION
!                            IN SUB OUTINT.
!      *MWP1T3G*    REAL     GRIDDED MEAN TRAIN 3 WAVE PERIOD 1 (S.). 

! ----------------------------------------------------------------------
      END MODULE YOWINTT
