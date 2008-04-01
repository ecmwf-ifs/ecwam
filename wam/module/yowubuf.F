      MODULE YOWUBUF

      IMPLICIT NONE

!*    ** *UBUF*  GRID POINT DEPENDENT CONSTANTS

      INTEGER, ALLOCATABLE :: KLAT(:,:) 
      INTEGER, ALLOCATABLE :: KLON(:,:) 
      INTEGER, ALLOCATABLE :: KLAT2(:,:) 
      INTEGER, ALLOCATABLE :: KCOR(:,:,:) 
      INTEGER, ALLOCATABLE :: KLATN(:,:,:) 
      INTEGER, ALLOCATABLE :: KLONN(:,:) 
      INTEGER, ALLOCATABLE :: KRLAT(:,:) 
      INTEGER, ALLOCATABLE :: KRLON(:,:) 
      INTEGER, ALLOCATABLE :: KRLAT2(:,:) 
      INTEGER, ALLOCATABLE :: KRLON2(:,:) 
      REAL, ALLOCATABLE :: WLAT(:,:) 
      REAL, ALLOCATABLE :: SUMWN(:,:,:) 
      REAL, ALLOCATABLE :: WLATN(:,:,:,:) 
      REAL, ALLOCATABLE :: WLONN(:,:,:) 
      REAL, ALLOCATABLE :: WKPMN(:,:,:,:) 
      REAL, ALLOCATABLE :: WRLAT(:,:) 
      REAL, ALLOCATABLE :: WRLON(:,:) 
      REAL, ALLOCATABLE :: OBSLAT(:,:,:) 
      REAL, ALLOCATABLE :: OBSLON(:,:,:) 
      REAL, ALLOCATABLE :: OBSCOR(:,:,:) 
      REAL, ALLOCATABLE :: OBSRLAT(:,:,:) 
      REAL, ALLOCATABLE :: OBSRLON(:,:,:) 
      LOGICAL, ALLOCATABLE :: LSAMEDEPTH(:)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *KLAT*      INTEGER   INDEX OF GRIDPOINT SOUTH AND NORTH
!                            LANDPOINTS ARE MARKED BY ZERO OR NINF-1.
!      *KLON*      INTEGER   INDEX OF GRIDPOINT WEST AND EAST
!                            LANDPOINTS ARE MARKED BY ZERO OR NINF-1.
!      *KLATN*     INTEGER   INDEX OF GRIDPOINT SOUTH AND NORTH
!                            FOR A GIVEN DIRECTION
!                            LANDPOINTS ARE MARKED BY ZERO OR NINF-1.
!      *KLONN*     INTEGER   INDEX OF GRIDPOINT WEST AND EAST
!                            FOR A GIVEN DIRECTION
!                            LANDPOINTS ARE MARKED BY ZERO OR NINF-1.
!      *KLAT2*     INTEGER   INDEX OF 2nd closest GRIDPOINT SOUTH AND
!                            NORTH LANDPOINTS ARE MARKED BY ZERO
!                            OR NINF-1.
!      *KCOR*      INTEGER   INDEX OF THE GRID CORNER POINTS IN DIAGONAL
!                            DIRECTIONS. IT IS USED IN MPDECOMP TO DETERMINE
!                            THE MAXIMUM SIZE OF THE HALO WHEN IPROPAGS=2.
!                            BECAUSE OF THE CTUW SCHEME EACH CORNER HAS TO BE
!                            SPECIFIED BY ITS CLOSEST GRID POINTS AND THE POINT
!                            ON EACH SIDE OF IT 
!                            THE SECOND INDEX IS 1 FOR NORTH EAST CORNER
!                                                2 FOR SOUTH EAST CORNER
!                                                3 FOR SOUTH WEST CORNER
!                                                4 FOR NORTH WEST CORNER
!                            THE THIRD INDEX IS 0 FOR THE CLOSEST GRID POINT
!                                               1 FOR THE SECOND CLOSEST GRID POINT
!                                              -1 FOR THE THIRD CLOSEST GRID POINT
!      *KRLAT*     INTEGER   INDEX OF GRIDPOINT SOUTH-EAST AND NORTH-WEST
!                            LANDPOINTS ARE MARKED BY ZERO OR NINF-1.
!      *KRLON*     INTEGER   INDEX OF GRIDPOINT SOUTH-WEST AND NORTH-EAST
!                            LANDPOINTS ARE MARKED BY ZERO OR NINF-1.
!      *KRLAT2*    INTEGER   INDEX OF 2nd closest GRIDPOINT SOUTH-EAST
!                            NORTH-WEST. LANDPOINTS ARE MARKED BY ZERO
!                            OR NINF-1.
!      *KRLON2*    INTEGER   INDEX OF 2nd closest GRIDPOINT SOUTH-WEST
!                            AND NORTH-EAST. LANDPOINTS ARE MARKED BY ZERO
!                            OR NINF-1.
!      *WLAT*      REAL      WEIGHT IN ADVECTION SCHEME FOR CLOSEST
!                            AND THE SECOND CLOSEST GRIDPOINTS IN THE
!                            NORTH-SOUTH DIRECTION.
!      *SUMWN*     REAL      SUM OF WEIGHTS IN ADVECTION CTU SCHEME
!      *WLATN*     REAL      WEIGHTS IN ADVECTION CTU SCHEME
!                            FOR NORTH-SOUTH DIRECTION.
!      *WLONN*     REAL      WEIGHTS IN ADVECTION CTU SCHEME
!                            FOR EAST-WEST DIRECTION.
!      *WKPMN*     REAL      WEIGHTS WHEN USING ADVECTION CTU SCHEME
!                            FOR REFRACTION TERM K+-1.
!      *WRLAT*     REAL      WEIGHT IN ADVECTION SCHEME FOR CLOSEST
!                            AND THE SECOND CLOSEST GRIDPOINTS IN THE
!                            SOUTH-EAST AND NORTH-WEST DIRECTION.
!      *WRLON*     REAL      WEIGHT IN ADVECTION SCHEME FOR CLOSEST
!                            AND THE SECOND CLOSEST GRIDPOINTS IN THE
!                            SOUTH-WEST AND NORTH-EAST DIRECTION.
!      *OBSLAT*    REAL      TRANSMISSION COEFFICIENT DUE TO OBSTRUCTIONS
!                            IN THE MERIDIONAL DIRECTION. 
!      *OBSLON*    REAL      TRANSMISSION COEFFICIENT DUE TO OBSTRUCTIONS
!                            IN THE ZONAL DIRECTION. 
!      *OBSCOR*    REAL      TRANSMISSION COEFFICIENT DUE TO OBSTRUCTIONS
!                            IN THE GRID CORNER DIRECTIONS. 
!      *OBSRLAT*   REAL      TRANSMISSION COEFFICIENT DUE TO OBSTRUCTIONS
!                            IN THE SOUTH-EAST AND NORTH-WEST DIRECTION. 
!      *OBSRLON*   REAL      TRANSMISSION COEFFICIENT DUE TO OBSTRUCTIONS
!                            IN THE SOUTH-WEST AND NORTH-EAST DIRECTION. 
!      *LSAMEDEPTH* LOGICAL  IF TRUE THEN ALL GRID POINTS THAT CAN BE 
!                            USED BY THE ADVECTION SCHEME FOR A GIVEN
!                            GRID POINT HAVE THE SAME DEPTH INDEX (AND
!                            THUS THE SAME GROUP SPEED).
! ----------------------------------------------------------------------
      END MODULE YOWUBUF
