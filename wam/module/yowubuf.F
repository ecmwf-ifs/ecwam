      MODULE YOWUBUF

      IMPLICIT NONE

!*    ** *UBUF*  GRID POINT DEPENDENT CONSTANTS

      INTEGER, ALLOCATABLE :: KLAT(:,:) 
      INTEGER, ALLOCATABLE :: KLON(:,:) 
      INTEGER, ALLOCATABLE :: KLAT2(:,:) 
      REAL, ALLOCATABLE :: WADV(:,:) 
      REAL, ALLOCATABLE :: OBSLAT(:,:,:) 
      REAL, ALLOCATABLE :: OBSLON(:,:,:) 
      LOGICAL, ALLOCATABLE :: LSAMEDEPTH(:)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *KLAT*      INTEGER   INDEX OF GRIDPOINT SOUTH AND NORTH
!                            LANDPOINTS ARE MARKED BY ZERO OR NINF-1.
!      *KLON*      INTEGER   INDEX OF GRIDPOINT WEST AND EAST
!                            LANDPOINTS ARE MARKED BY ZERO OR NINF-1.
!      *KLAT2*     INTEGER   INDEX OF 2nd closest GRIDPOINT SOUTH AND
!                            NORTH LANDPOINTS ARE MARKED BY ZERO
!                            OR NINF-1.
!      *WADV*      REAL      WEIGHT IN ADVECTION SCHEME FOR CLOSEST
!                            AND THE SECOND CLOSEST GRIDPOINTS.
!      *OBSLAT*    REAL      REDUCTION FACTOR DUE TO OBSTRUCTIONS IN
!                            THE MERIDIONAL DIRECTION. 
!      *OBSLON*    REAL      REDUCTION FACTOR DUE TO OBSTRUCTIONS IN
!                            THE ZONAL DIRECTION. 
!      *LSAMEDEPTH* LOGICAL  IF TRUE THEN ALL GRID POINTS THAT CAN BE 
!                            USED BY THE ADVECTION SCHEME FOR A GIVEN
!                            GRID POINT HAVE THE SAME DEPTH INDEX (AND
!                            THUS THE SAME GROUP SPEED).

! ----------------------------------------------------------------------
      END MODULE YOWUBUF
