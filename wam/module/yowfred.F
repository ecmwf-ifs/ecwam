      MODULE YOWFRED

      IMPLICIT NONE

!*    ** *FREDIR* - FREQUENCY AND DIRECTION GRID.

      REAL, ALLOCATABLE :: FR(:)
      REAL, ALLOCATABLE :: DFIM(:) 
      REAL, ALLOCATABLE :: GOM(:) 
      REAL, ALLOCATABLE :: C(:) 
      REAL DELTH 
      REAL DELTR 
      REAL, ALLOCATABLE :: TH(:) 
      REAL, ALLOCATABLE :: COSTH(:) 
      REAL, ALLOCATABLE :: SINTH(:) 
      REAL, PARAMETER :: FRATIO = 1.1
      REAL, PARAMETER :: FRTAIL = 0.2

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *FR*        REAL      FREQUENCIES IN HERTZ.
!      *DFIM*      REAL      FREQUENCY INTERVAL*DIRECTION INTERVAL.
!      *GOM*       REAL      DEEP WATER GROUP VELOCITIES (M/S).
!      *C*         REAL      DEEP WATER PHASE VELOCITIES (M/S).
!      *DELTH*     REAL      ANGULAR INCREMENT OF SPECTRUM (RADIANS).
!      *DELTR*     REAL      DELTH TIMES RADIUS OF EARTH (METRES).
!      *TH*        REAL      DIRECTIONS IN RADIANS.
!      *COSTH*     REAL      COS OF DIRECTION.
!      *SINTH*     REAL      SIN OF DIRECTION.
!      *FRATIO*    REAL      FREQUENCY RATIO.
!      *FRTAIL*    REAL      FREQUENCY TAIL CONSTANT FACTOR.

! ----------------------------------------------------------------------
      END MODULE YOWFRED
