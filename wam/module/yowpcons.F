      MODULE YOWPCONS

      IMPLICIT NONE

!*    *PARAMETER* OF GLOBAL CONSTANTS.

      REAL, PARAMETER :: G = 9.806
      REAL, PARAMETER :: PI = 3.1415927
      REAL, PARAMETER :: CIRC = 40000000.
      REAL, PARAMETER :: ZPI = 2.*PI
      REAL, PARAMETER :: RAD = PI/180.
      REAL, PARAMETER :: DEG = 180./PI
      REAL, PARAMETER :: R = CIRC/ZPI
      REAL, PARAMETER :: EPSMIN=0.1E-32
      REAL, PARAMETER :: ZMISS=-999.0
      REAL, PARAMETER :: ROAIR = 1.225
      REAL, PARAMETER :: ROWATER = 1000.
      REAL, PARAMETER :: YEPS = ROAIR/ROWATER
      REAL, PARAMETER :: YINVEPS = 1./YEPS


!*    VARIABLE.   TYPE.     PURPOSE.
!     ---------   -------   --------
!     *G*         REAL      ACCELLERATION OF GRAVITY.
!     *PI*        REAL      PI.
!     *CIRC*      REAL      EARTH CIRCUMFERENCE (METRES).
!     *RAD*       REAL      PI / 180.
!     *DEG*       REAL      180. / PI.
!     *ZPI*       REAL      2. * PI.
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

