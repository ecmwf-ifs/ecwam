      MODULE YOWPCONS

      IMPLICIT NONE

!*    *PARAMETER* OF GLOBAL CONSTANTS.
!     THE NON PARAMETRIC CONSTANTS WILL BE RESET BY CALLING INIWCST !!!!

      REAL, PARAMETER :: G = 9.806
      REAL, PARAMETER :: OLDPI = 3.1415927 
      REAL :: PI = OLDPI 
      REAL, PARAMETER :: CIRC = 40000000.
      REAL :: ZPI = 6.2831854
      REAL :: ZCONST = 0.0281349 
      REAL :: RAD = 0.017453293 
      REAL :: DEG = 57.295778667 
      REAL :: R = 6366198. 
      REAL, PARAMETER :: EPSMIN = 0.1E-32
      REAL, PARAMETER :: ZMISS = -999.0
      REAL, PARAMETER :: ROAIR = 1.225
      REAL, PARAMETER :: ROWATER = 1000.
      REAL, PARAMETER :: YEPS = ROAIR/ROWATER
      REAL, PARAMETER :: YINVEPS = 1./YEPS

!*    VARIABLE.   TYPE.     PURPOSE.
!     ---------   -------   --------
!     *G*         REAL      ACCELLERATION OF GRAVITY.
!     *OLDPI*     REAL      OLD VALUE USED FOR PI (PRIOR TO CY21R3). 
!     *PI*        REAL      PI (See SUB. INIWCST).
!     *CIRC*      REAL      EARTH CIRCUMFERENCE (METRES).
!     *RAD*       REAL      PI / 180.
!     *DEG*       REAL      180. / PI.
!     *ZPI*       REAL      2. * PI.
!     *ZCONST*    REAL      1./(8.*PI*SQRT(2.))
!     *R*         REAL      EARTH RADIUS        (METRES).
!     *EPSMIN*    REAL      SMALL NUMBER
!     *ZMISS*     REAL      GRIB DATA MISSING INDICATOR.
!     *ROAIR*     REAL      AIR DENSITY.
!     *ROWATER*   REAL      WATER DENSITY.
!     *YEPS*      REAL      ROAIR/ROWATER.
!     *YINVEPS    REAL      1./YEPS.
!
!----------------------------------------------------------------------

      END MODULE YOWPCONS

