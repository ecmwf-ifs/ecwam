      MODULE YOWCBOU

      IMPLICIT NONE

!*    ** *CBOUND* COURSE GRID BOUNDARY OUTPUT.

      INTEGER IBOUNC 
      INTEGER NBOUNC 
      INTEGER, ALLOCATABLE :: IJARC(:) 
      INTEGER, ALLOCATABLE :: IGARC(:) 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -----     --------
!      *IBOUNC*    REAL      FLAG FOR THE COURSE GRID
!                            = 1; THE RUN INCLUDEDS BOUNDARY POINTS.
!                            ELSE; NO BOUNDARY POINTS.
!      *NBOUNC*    INTEGER   NUMBER OF BOUNDARY POINTS.
!      *IGARC*     INTEGER   INDEX OF BLOCK FOR A BOUNDARY POINT.
!      *IJARC*     INTEGER   INDEX IN A BLOCK FOR A BOUNDARY POINT.
! ----------------------------------------------------------------------
      END MODULE YOWCBOU
