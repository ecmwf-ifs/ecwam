      MODULE YOWPCONS

      IMPLICIT NONE

!*    *PARAMETER* OF GLOBAL CONSTANTS.
!     THE NON PARAMETRIC CONSTANTS WILL BE RESET BY CALLING INIWCST !!!!

      REAL :: G = 9.806
      REAL, PARAMETER :: OLDPI = 3.1415927 
      REAL :: PI = OLDPI 
      REAL, PARAMETER :: CIRC = 40000000.
      REAL :: ZPI = 6.2831854
      REAL :: G2ZPI4 = 0.987152585
      REAL :: ZPISQRT = 1.7724539
      REAL :: ZCONST = 0.0281349 
      REAL :: RAD = 0.017453293 
      REAL :: DEG = 57.295778667 
      REAL :: R = 6366198. 
      REAL, PARAMETER :: EPSMIN = 0.1E-32
      REAL, PARAMETER :: DKMAX = 40.
      REAL, PARAMETER :: BF2MIN = -10.
      REAL, PARAMETER :: BF2MAX = 10.
      REAL, PARAMETER :: C4MIN = -0.33
      REAL, PARAMETER :: C4MAX = 1.
      REAL, PARAMETER :: PHIEPSMIN = -3276.80
      REAL, PARAMETER :: PHIEPSMAX = -0.05
      REAL :: ZMISS=-999.0
      REAL, PARAMETER :: ROAIR = 1.225
      REAL, PARAMETER :: ROWATER = 1000.
      REAL, PARAMETER :: YEPS = ROAIR/ROWATER
      REAL, PARAMETER :: YINVEPS = 1./YEPS
      REAL, PARAMETER :: WSTAR0 = 0.
      REAL, PARAMETER :: EGRCRV = 1877
      REAL, PARAMETER :: AGRCRV = 0.544E+6
      REAL, PARAMETER :: BGRCRV = 1.90
      REAL, PARAMETER :: AFCRV = 5.054E-4
      REAL, PARAMETER :: BFCRV = -2.959
      REAL, PARAMETER :: ESH = 1711.
      REAL, PARAMETER :: ASH = 8.0E-4 
      REAL, PARAMETER :: BSH = 0.96 
      REAL, PARAMETER :: ASWKM=0.0981
      REAL, PARAMETER :: BSWKM=0.425
      REAL, PARAMETER :: Rconstant = 287.16   ! The gas constant
      REAL, PARAMETER :: EpsWaterVapor = 0.61 ! The mass ratio of 
                                              ! water vapor to dry air
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
!     *G2ZPI4*    REAL      G**2/ZPI**4 
!     *ZPISQRT*   REAL      SQRT(PI)
!     *ZCONST*    REAL      1./(8.*PI*SQRT(2.))
!     *R*         REAL      EARTH RADIUS        (METRES).
!     *EPSMIN*    REAL      SMALL NUMBER
!     *DKMAX*     REAL      MAXIMUM VALUE OF DEPTH*WAVENUMBER BEFORE
!                           SIMPLIFYING TO DEEP WATER EXPRESIONS.
!     *BF2MAX*    REAL      MAXIMUM VALUE ALLOWED FOR BFI SQUARED
!     *BF2MIN*    REAL      MINIMUM VALUE ALLOWED FOR BFI SQUARED
!     *C4MAX*     REAL      MAXIMUM VALUE ALLOWED FOR KURTOSIS
!     *C4MIN*     REAL      MINIMUM VALUE ALLOWED FOR KUTOSIS 
!     *PHIEPSMIN* REAL      MINIMUM VALUE ALLOWED FOR NORMALISED ENERGY FLUX INTO OCEAN
!     *PHIEPSMAX* REAL      MAXIMUM VALUE ALLOWED FOR NORMALISED ENERGY FLUX INTO OCEAN
!     *ZMISS*     REAL      MISSING DATA INDICATOR 
!                          (SET IN CHIEF OR VIA THE IFS).
!     *ROAIR*     REAL      AIR DENSITY.
!     *ROWATER*   REAL      WATER DENSITY.
!     *YEPS*      REAL      ROAIR/ROWATER.
!     *YINVEPS*   REAL      1./YEPS.
!     *WSTAR0*    REAL      DEFAULT VALUE FOR w*.
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
!     *ESH*       REAL      PARAMETER OF THE NON DIMENSIONAL ENERGY LIMIT
!                           AS A FUNCTION OF NON-DIMENSIONAL DEPTH DSTAR
!                           ESTAR_LIMIT=ESH*TANH(ASH*DSTAR**BSH), where
!                           DSTAR=DEPTH*G/USTAR**2
!     *ASH*       REAL      PARAMETER OF THE NON DIMENSIONAL ENERGY LIMIT
!                           SEE ESH.
!     *BSH*       REAL      PARAMETER OF THE NON DIMENSIONAL ENERGY LIMIT
!                           SEE ESH.
!     *ASWKM*     REAL      PARAMETER TO EXPRESS MINIMUM KMEAN AS A FUNCTION
!                           OF DEPTH: KMEAN >= ASWKM/DEPTH**BSWKM 
!     *BSWKM*     REAL      PARAMETER TO EXPRESS MINIMUM KMEAN AS A FUNCTION
!                           OF DEPTH: KMEAN >= ASWKM/DEPTH**BSWKM 
!     *EPSUS*     REAL      SMALL NUMBER ADDED TO THE SQUARE OF USTAR
!                           WHEN COMPUTING CD OR THE CHARNOCK COEF.
!     *EPSU10*    REAL      SMALL NUMBER ADDED TO THE SQUARE OF U10
!                           WHEN COMPUTING CD.
!----------------------------------------------------------------------

      END MODULE YOWPCONS

