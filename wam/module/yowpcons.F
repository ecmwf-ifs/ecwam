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
      REAL, PARAMETER :: EGRCRV = 1736.2
      REAL, PARAMETER :: AGRCRV = 1.915E+6
      REAL, PARAMETER :: BGRCRV = 0.907
      REAL, PARAMETER :: AFCRV = 7.2431E-4
      REAL, PARAMETER :: BFCRV = -2.8659
      REAL, PARAMETER :: Rconstant = 287.16   ! The gas constant
      REAL, PARAMETER :: EpsWaterVapor = 0.61 ! The mass ratio of 
                                              ! water vapor to dry air
      REAL, PARAMETER :: AFMRV = 8.1301E-4
      REAL, PARAMETER :: BFMRV = -2.7830
      REAL, PARAMETER :: EPSUS = 1.0E-6
      REAL, PARAMETER :: EPSU10 = 1.0E-3

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
!     *YINVEPS*   REAL      1./YEPS.
!     *EGRCRV*    REAL      PARAMETER OF THE NON DIMENSIONAL ENERGY 
!                           GROWTH CURVE.
!                           ESTAR=EGRCRV*(TSTAR/(AGRCRV+TSTAR))**BGRCRV
!     *AGRCRV*    REAL      PARAMETER OF THE NON DIMENSIONAL ENERGY 
!                           GROWTH CURVE.
!     *BGRCRV*    REAL      PARAMETER OF THE NON DIMENSIONAL ENERGY 
!                           GROWTH CURVE.
!     *AFCRV*     REAL      PARAMETER OF THE NON DIMENSIONAL FREQUENCY 
!                           GROWTH CURVE. Based on f * F(f)
!                           FSTAR=(EGRCRV/AFCRV)**(1/BFCRV)
!     *BFCRV*     REAL      PARAMETER OF THE NON DIMENSIONAL FREQUENCY
!                           GROWTH CURVES.
!     *AFMRV*     REAL      PARAMETER OF THE NON DIMENSIONAL FREQUENCY 
!                           GROWTH CURVE. Based on 1/f F(f)
!                           FSTAR=(EGRCRV/AFCRV)**(1/BFCRV)
!     *BFMRV*     REAL      PARAMETER OF THE NON DIMENSIONAL FREQUENCY
!                           GROWTH CURVE.
!     *EPSUS*     REAL      SMALL NUMBER ADDED TO THE SQUARE OF USTAR
!                           WHEN COMPUTING CD OR THE CHARNOCK COEF.
!     *EPSU10*    REAL      SMALL NUMBER ADDED TO THE SQUARE OF U10
!                           WHEN COMPUTING CD.
!----------------------------------------------------------------------

      END MODULE YOWPCONS

