      MODULE YOWUBUF

      IMPLICIT NONE

!*    ** *UBUF*  GRID POINT DEPENDENT CONSTANTS

      INTEGER, ALLOCATABLE :: KLAT(:,:) 
      INTEGER, ALLOCATABLE :: KLON(:,:) 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *KLAT*      INTEGER   INDEX OF GRIDPOINT SOUTH AND NORTH
!                            LANDPOINTS ARE MARKED BY ZERO OR NINF-1.
!      *KLON*      INTEGER   INDEX OF GRIDPOINT WEST AND EAST
!                            LANDPOINTS ARE MARKED BY ZERO OR NINF-1.

! ----------------------------------------------------------------------
      END MODULE YOWUBUF
