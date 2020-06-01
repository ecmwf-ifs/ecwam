FUNCTION GAMMA_WAM(OMEGA, XK, USTAR, Z0, EPS)
 
!---------------------------------------------------------------------
 
!**** *GAMMA_WAM* - COMPUTATION OF GROWTHRATE
 
!     P.A.E.M. JANSSEN
 
!     PURPOSE.
!     ---------
 
!     COMPUTES GROWTHRATE BY WIND
 
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

REAL(KIND=JWRB), INTENT(IN) :: OMEGA  ! angular frequency
REAL(KIND=JWRB), INTENT(IN) :: XK     ! wave number
REAL(KIND=JWRB), INTENT(IN) :: USTAR  ! friction velocity
REAL(KIND=JWRB), INTENT(IN) :: Z0     ! roughness length
REAL(KIND=JWRB), INTENT(IN) :: EPS    ! ratio of air density to water density

REAL(KIND=JWRB) :: CM, ZFAK, X, XLOG, ZLOG, ZLOG2X, ZBETA
REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GAMMA_WAM',0,ZHOOK_HANDLE)

!* 1. DETERMINE GROWTH ACCORDING TO JANSSEN-MILES
!     -------------------------------------------
 
CM   = XK/OMEGA
ZFAK = OMEGA**2/(G*XK)

X       = USTAR*CM 
XLOG    = LOG(XK*Z0) + XKAPPA/(X + ZALP) 
ZLOG    = MIN(XLOG,0.0_JWRB)
ZLOG2X  = ZLOG*ZLOG*X
ZBETA   = EPS*BETAMAXOXKAPPA2*EXP(ZLOG)*ZLOG2X**2

GAMMA_WAM = ZBETA*ZFAK*OMEGA
 
! -----------------------------------------------------------------
 
IF (LHOOK) CALL DR_HOOK('GAMMA_WAM',1,ZHOOK_HANDLE)

END FUNCTION GAMMA_WAM 
