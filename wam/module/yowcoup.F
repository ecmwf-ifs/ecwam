      MODULE YOWCOUP

      IMPLICIT NONE

!*    ** *COUPL* - PARAMETERS FOR COUPLING.

      INTEGER, PARAMETER :: JPLEVC=1 
      LOGICAL LWCOU 
      LOGICAL LWCOU2W 
      REAL BETAMAX 
      REAL ZALP 
      REAL ALPHA 
      REAL XKAPPA 
      REAL XNLEV(JPLEVC) 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *LWCOU*     LOGICAL   CONTROLS COUPLING WITH ATMOSPHERIC MODEL
!      *LWCOU2W*   LOGICAL   CONTROLS 1-WAY OR 2-WAY COUPLING
!      *BETAMAX*   REAL      PARAMETER FOR WIND INPUT.
!      *ZALP*      REAL      SHIFTS GROWTH CURVE.
!      *ALPHA*     REAL      CHARNOCK CONSTANT.
!      *XKAPPA*    REAL      VON KARMAN CONSTANT.
!      *XNLEV*     REAL      WIND SPEED REFERENCE HEIGHT

! ----------------------------------------------------------------------
      END MODULE YOWCOUP
