      MODULE YOWFRED

      IMPLICIT NONE

!*    ** *FREDIR* - FREQUENCY AND DIRECTION GRID.

      REAL, ALLOCATABLE :: FR(:)
      REAL, ALLOCATABLE :: DFIM(:) 
      REAL, ALLOCATABLE :: DFIMOFR(:) 
      REAL, ALLOCATABLE :: GOM(:) 
      REAL, ALLOCATABLE :: C(:) 
      REAL :: DELTH 
      REAL :: DELTR 
      REAL, ALLOCATABLE :: TH(:) 
      REAL, ALLOCATABLE :: COSTH(:) 
      REAL, ALLOCATABLE :: SINTH(:) 
      REAL, ALLOCATABLE :: FR5(:)
      REAL, ALLOCATABLE :: FRM5(:)
      REAL, ALLOCATABLE :: COFRM4(:)

      REAL, PARAMETER :: FRATIO = 1.1
      REAL, PARAMETER :: WETAIL = 0.25
      REAL, PARAMETER :: FRTAIL = 0.2
      REAL, PARAMETER :: WP1TAIL = 1./3. 
      REAL, PARAMETER :: WP2TAIL = 0.5
      REAL, PARAMETER :: COEF4 = 5.E-07 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *FR*        REAL      FREQUENCIES IN HERTZ.
!      *DFIM*      REAL      FREQUENCY INTERVAL*DIRECTION INTERVAL.
!      *DFIMOFR*   REAL      DFIM/FR
!      *GOM*       REAL      DEEP WATER GROUP VELOCITIES (M/S).
!      *C*         REAL      DEEP WATER PHASE VELOCITIES (M/S).
!      *DELTH*     REAL      ANGULAR INCREMENT OF SPECTRUM (RADIANS).
!      *DELTR*     REAL      DELTH TIMES RADIUS OF EARTH (METRES).
!      *TH*        REAL      DIRECTIONS IN RADIANS.
!      *COSTH*     REAL      COS OF DIRECTION.
!      *SINTH*     REAL      SIN OF DIRECTION.
!      *FR5*       REAL      FR(M)**5 
!      *FRM5*      REAL      (1./FR(M))**5 
!      *COFRM4*    REAL      COEF4*G*FR(M)**(-4.) 
!      *FRATIO*    REAL      FREQUENCY RATIO.
!      *WETAIL*    REAL      WAVE ENERGY TAIL CONSTANT FACTOR.
!      *FRTAIL*    REAL      FREQUENCY TAIL CONSTANT FACTOR.
!      *WP1TAIL*   REAL      PERIOD 1 TAIL CONSTANT FACTOR.
!      *WP2TAIL*   REAL      PERIOD 2 TAIL CONSTANT FACTOR.
!      *COEF4*     REAL      COEFFICIENT USED TO COMPUTE THE SPECTRAL
!                            LIMITER IN IMPLSCH.

! ----------------------------------------------------------------------
      END MODULE YOWFRED
