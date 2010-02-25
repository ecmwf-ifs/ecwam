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
      REAL, ALLOCATABLE :: WHGTAG(:,:) 
      REAL, ALLOCATABLE :: CWHTAG(:,:) 
      REAL, ALLOCATABLE :: RANGAG(:,:) 
      REAL, ALLOCATABLE :: MWP1G(:,:) 
      REAL, ALLOCATABLE :: MWP2G(:,:) 
      REAL, ALLOCATABLE :: WSPRDG(:,:) 
      REAL, ALLOCATABLE :: C4G(:,:) 
      REAL, ALLOCATABLE :: BFG(:,:) 
      REAL, ALLOCATABLE :: QPG(:,:) 
      REAL, ALLOCATABLE :: DEPTHG(:,:) 
      REAL, ALLOCATABLE :: HMAXG(:,:) 
      REAL, ALLOCATABLE :: TMAXG(:,:) 
      REAL, ALLOCATABLE :: USTOKESG(:,:) 
      REAL, ALLOCATABLE :: VSTOKESG(:,:) 
      REAL, ALLOCATABLE :: UCURG(:,:) 
      REAL, ALLOCATABLE :: VCURG(:,:) 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *WHGTTG*    REAL      GRIDDED MEAN ENERGY (M*M) CHANGED TO
!                            WAVE HEIGHT (M) IN SUB MAKEGRID.
!      *WDIRTG*    REAL      GRIDDED MEAN DIRECTION (RAD) CHANGED TO
!                            DEGREE AND METEOROLOGICAL DEFINITION
!                            IN SUB OUTINT.
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
!      *WHGTAG*    REAL      GRIDDED ALTIMETER WAVE HEIGHT (M).
!      *CWHTAG*    REAL      GRIDDED CORRECTED ALTIMETER WAVE HEIGHT (M) 
!      *RANGAG*    REAL      GRIDDED RELATIVE CORRECTION TO THE
!                            ALTIMETER RANGE.
!      *MWP1G*     REAL      GRIDDED MEAN WAVE PERIOD 1 (S.).
!      *MWP2G*     REAL      GRIDDED MEAN WAVE PERIOD 2 (S.).
!      *WSPRDG*    REAL      GRIDDED DIRECTIONAL SPREAD (-). 
!      *C4G*       REAL      GRIDDED SPECTRAL KURTOSIS (-). 
!      *BFG*       REAL      GRIDDED BENJAMIN-FEIR INDEX (-). 
!      *QPG*       REAL      GRIDDED SPECTRAL PEAKEDNESS (-). 
!      *DEPTHG*    REAL      GRIDDED MODEL DEPTH.
!      *HMAXG*     REAL      GRIDDED MAXIMUM WAVE HEIGHT FIELD (M).
!      *TMAXG*     REAL      GRIDDED MAXIMUM WAVE PERIOD FIELD (S).
!      *USTOKESG*  REAL      GRIDDED U-COMP STOKES DRIFT (M/S).
!      *VSTOKESG*  REAL      GRIDDED V-COMP STOKESDRIFT (M/S).
!      *UCURG*     REAL      GRIDDED U-COMP SURFACE CURRENT (M/S).
!      *VCURG*     REAL      GRIDDED V-COMP SURFACE CURRENT (M/S).


! ----------------------------------------------------------------------
      END MODULE YOWINTP
