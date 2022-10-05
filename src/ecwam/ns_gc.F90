INTEGER FUNCTION NS_GC (USTAR)

! ----------------------------------------------------------------------

!**** *NS_GC* - FUNCTION TO DETERMINE THE CUT-OFF ANGULAR FREQUENCY INDEX 
!               FOR THE GRAVITY-CAPILLARY MODEL
!               !!!! rounded to the closest index of XK_GC  !!!!!

!**   INTERFACE.
!     ----------

!       *FUNCTION* *NS_GC (USTAR)*

!       *USTAR*  - FRICTION VELOCITY. 

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWFRED  , ONLY : XLOGKRATIOM1_GC, NWAV_GC,  XKM_GC
USE YOWPCONS , ONLY : SQRTGOSURFT

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN) :: USTAR

      REAL(KIND=JWRB) :: Y, XKS

! ----------------------------------------------------------------------


!!!Y = 1.0_JWRB/(1.48_JWRB+2.05_JWRB*UST)
!!!Y = (1.0_JWRB + UST**2)/(1.0_JWRB+10.0_JWRB*UST**2)

XKS = SQRTGOSURFT/(1.48_JWRB+2.05_JWRB*USTAR)

NS_GC = MIN( INT(LOG(XKS*XKM_GC(1))*XLOGKRATIOM1_GC) + 1, NWAV_GC-1)


END FUNCTION NS_GC
