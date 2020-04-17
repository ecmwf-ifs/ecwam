      SUBROUTINE STRESS_GC(UST, EPS, Z0, ALPHAP, XMSSCG, TAUWCG)

!***  DETERMINE MSS AND WAVE INDUCED STRESS FOR GRAV-CAP WAVES

!     AUTHOR: PETER JANSSEN
!     ------

!     REFERENCES:
!     ----------

!     VIERS PAPER EQ.(29)

!----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LLGCBZ0
      USE YOWFRED  , ONLY : NWAV_GC, KRATIO_GC, OMEGA_GC, XK_GC, VG_GC, C_GC
      USE YOWPCONS , ONLY : G, ZPI

      USE YOMHOOK  ,ONLY : LHOOK,   DR_HOOK

!----------------------------------------------------------------------

      IMPLICIT NONE

#include "omegagc.intfb.h"

      REAL(KIND=JWRB), INTENT(IN) :: UST, EPS, Z0, ALPHAP
      REAL(KIND=JWRB), INTENT(OUT) :: XMSSCG, TAUWCG

      INTEGER(KIND=JWIM) :: NS
      INTEGER(KIND=JWIM) :: I

      REAL(KIND=JWRB), PARAMETER :: ANG = 0.5_JWRB

      REAL(KIND=JWRB) :: GAMMA_WAM

      REAL(KIND=JWRB) :: XKS, OMS
      REAL(KIND=JWRB) :: X, X3, BB, EPS0ALPHA3M
      REAL(KIND=JWRB) :: SUMT, SUMS, GAM_W, BS, ALPHA3, OM, F1
      REAL(KIND=JWRB) :: DELK_GC(NWAV_GC)
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
   
!     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
#include "gc_dispersion.h"

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('STRESS_GC',0,ZHOOK_HANDLE)


!*    1.0  DETERMINE GRAV_CAP SPECTRUM, MSS AND TAUWHF.
!          -------------------------------------------

      ALPHA3 = 3._JWRB*ZPI

      CALL OMEGAGC(UST, NS, XKS, OMS)

      BS  = 0.5_JWRB*ALPHAP
      EPS0ALPHA3M = (ALPHA3*FC_GC(XKS)**4/FVG_GC(XKS)*BS**2)/ALPHA3   
      SUMS = 0.0_JWRB
      SUMT = 0.0_JWRB

      DELK_GC(NS) = 0.5*(XK_GC(NS+1)-XK_GC(NS))
      DO I=NS+1,NWAV_GC-1
        DELK_GC(N) = 0.5_JWRB*(XK_GC(I+1)-XK_GC(I-1))
      ENDDO
      DELK_GC(NWAV_GC) = 0.5_JWRB*(XK_GC(NWAV_GC)-XK_GC(NWAV_GC-1))

      DO I = NS, NWAV_GC
        X     = XK_GC(I)
        X3    = X**3
!       ANALYTICAL FORM INERTIAL SUB RANGE
        BB    = SQRT(EPS0ALPHA3M*VG_GC(I))/C_GC(I)**2
        F1    = BB/(X*X3)
        SUMS  = SUMS + X3 * F1 * DELK_GC(I)
        GAM_W = ANG * GAMMA_WAM(OMEGA_GC(I), X, UST, Z0, EPS)
        SUMT  = SUMT + ( OMEGA_GC(I) * GAM_W * BB * DELK_GC(I) ) / (X3*EPS)
      ENDDO
 
      XMSSCG = SUMS
      TAUWCG = SUMT

      IF (LHOOK) CALL DR_HOOK('STRESS_GC',1,ZHOOK_HANDLE)
 
      END SUBROUTINE STRESS_GC
