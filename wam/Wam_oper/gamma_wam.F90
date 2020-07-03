FUNCTION GAMMA_WAM(XK, C, OM3GMKM, USTAR, Z0)
 
!---------------------------------------------------------------------
 
!**** *GAMMA_WAM* - COMPUTATION OF GROWTHRATE  !!! but without ratio of air density to water density
 
!     P.A.E.M. JANSSEN
 
!     PURPOSE.
!     ---------
 
!     COMPUTES GROWTHRATE BY WIND WITHOUT the multiplicative representing the ratio of air density to water density (eps)
 
!**   INTERFACE.
!     ----------
 
!          *FUNCTION CALL*
 
!     REFERENCES.
!     -----------
 
!          FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
 
!-------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWPCONS , ONLY : G
USE YOWPHYS  , ONLY : XKAPPA, ZALP,   BETAMAXOXKAPPA2
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JWRB) :: GAMMA_WAM 

REAL(KIND=JWRB), INTENT(IN) :: XK     ! wave number
REAL(KIND=JWRB), INTENT(IN) :: C      ! phase speed
REAL(KIND=JWRB), INTENT(IN) :: OM3GMKM   ! OMEGA * OMEGA**2/(g*xk) 
REAL(KIND=JWRB), INTENT(IN) :: USTAR  ! friction velocity
REAL(KIND=JWRB), INTENT(IN) :: Z0     ! roughness length

REAL(KIND=JWRB) :: X, XLOG, ZLOG, ZLOG2X, ZBETA
REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GAMMA_WAM',0,ZHOOK_HANDLE)

!* 1. DETERMINE GROWTH ACCORDING TO JANSSEN-MILES
!     -------------------------------------------
 
X       = USTAR/C
XLOG    = LOG(XK*Z0) + XKAPPA/(X + ZALP) 
ZLOG    = MIN(XLOG,0.0_JWRB)
ZLOG2X  = ZLOG*ZLOG*X
ZBETA   = BETAMAXOXKAPPA2*EXP(ZLOG)*ZLOG2X**2

GAMMA_WAM = ZBETA*OM3GMKM
 
! -----------------------------------------------------------------
 
IF (LHOOK) CALL DR_HOOK('GAMMA_WAM',1,ZHOOK_HANDLE)

END FUNCTION GAMMA_WAM
