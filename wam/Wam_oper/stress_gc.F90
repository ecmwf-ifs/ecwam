      SUBROUTINE STRESS_GC(UST, EPS, Z0, ALPHAP, OMS, XMSS, TAUW)

!***  DETERMINE MSS AND WAVE INDUCED STRESS FOR GRAV-CAP WAVES

!     AUTHOR: PETER JANSSEN
!     ------

!     REFERENCES:
!     ----------

!     VIERS PAPER EQ.(29)

!----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LLGCBZ0
      USE YOWFRED  , ONLY : NWAV_GC, KRATIO_GC, OMEGA_GC, XK_GC, DELK_GC
      USE YOWPCONS , ONLY : G, ZPI, SURFT

      USE YOMHOOK  ,ONLY : LHOOK,   DR_HOOK

!----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN) :: UST, EPS, Z0, ALPHAP
      REAL(KIND=JWRB), INTENT(OUT) :: OMS, XMSS, TAUW

      INTEGER(KIND=JWIM) :: I, NS 

      REAL(KIND=JWRB), PARAMETER :: ANG = 0.5_JWRB

      REAL(KIND=JWRB) :: GAMMA_WAM

      REAL(KIND=JWRB) :: X, X3, Y, BB, XKS, XKP, EPS0
      REAL(KIND=JWRB) :: SUMT, SUMS, GAM_W
      REAL(KIND=JWRB) :: BS, ALPHA3, ALPHA3M, XK0, OM, F1
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
   
!     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
#include "gc_dispersion.h"

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('STRESS_GC',0,ZHOOK_HANDLE)


!*    1.0  DETERMINE GRAV_CAP SPECTRUM, MSS AND TAUWHF.
!          -------------------------------------------

      ALPHA3 = 3._JWRB*ZPI
      ALPHA3M = 1.0_JWRB / ALPHA3
      XK0 = SQRT(G/SURFT)    

      Y = 1.0_JWRB/(1.48_JWRB+2.05_JWRB*UST)
      XKS = XK0*Y
      NS = LOG(XKS/XK_GC(1))/LOG(KRATIO_GC)+1
      OMS = OMEG_GC(XKS)

      BS  = 0.5_JWRB*ALPHAP
      EPS0 = ALPHA3*C(XKS)**4/VG_GC(XKS)*BS**2   
      SUMS = 0.0_JWRB
      SUMT = 0.0_JWRB

      DO I = NS, NWAV_GC
        X     = XK_GC(I)
        X3    = X**3
!       ANALYTICAL FORM INERTIAL SUB RANGE
        BB    = SQRT(EPS0*VG_GC(X)*ALPHA3M)/C_GC(X)**2
        F1 = BB/(X*X3)
        SUMS   = SUMS + X3 * F1 * DELK_GC(I)
        GAM_W = ANG * GAMMA_WAM(OMEGA_GC(I), X, UST, Z0, EPS)
        SUMT  = SUMT + OMEG_GC(X) * GAM_W * BB * DELK_GC(I) /X3 /EPS
      ENDDO
 
      XMSS = SUMS
      TAUW = SUMT

      IF (LHOOK) CALL DR_HOOK('STRESS_GC',1,ZHOOK_HANDLE)
 
      END SUBROUTINE STRESS_GC
