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

USE YOWFRED  , ONLY : XLOGKRATIO_GC, NWAV_GC,  XK_GC
USE YOWPCONS , ONLY : SQRTGOSURFT

USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN) :: USTAR

      REAL(KIND=JWRB) :: Y, XKS
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('NS_GC',0,ZHOOK_HANDLE)

!!!Y = 1.0_JWRB/(1.48_JWRB+2.05_JWRB*UST)
!!!Y = (1.0_JWRB + UST**2)/(1.0_JWRB+10.0_JWRB*UST**2)

XKS = SQRTGOSURFT/(1.48_JWRB+2.05_JWRB*USTAR)

NS_GC = MIN( INT(LOG(XKS/XK_GC(1))/XLOGKRATIO_GC) + 1, NWAV_GC-1)

IF (LHOOK) CALL DR_HOOK('NS_GC',1,ZHOOK_HANDLE)

END FUNCTION NS_GC
