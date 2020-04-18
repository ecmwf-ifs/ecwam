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
      USE YOWPCONS , ONLY : G, THREEZPI, SURFT

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
      REAL(KIND=JWRB) :: XM, BBDELK, EPS0THREEZPIM
      REAL(KIND=JWRB) :: SUMT, SUMS, GAM_W, BS, OM, EPSM1
      REAL(KIND=JWRB) :: DELK_GC(NWAV_GC)
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
   
!     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
#include "gc_dispersion.h"

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('STRESS_GC',0,ZHOOK_HANDLE)


!*    1.0  DETERMINE GRAV_CAP SPECTRUM, MSS AND TAUWHF.
!          -------------------------------------------

      CALL OMEGAGC(UST, NS, XKS, OMS)

      EPSM1 = 1.0_JWRB/EPS
      BS  = 0.5_JWRB*ALPHAP
      EPS0THREEZPIM = (THREEZPI*FC_GC(XKS)**4/FVG_GC(XKS)*BS**2)/THREEZPI   
      SUMS = 0.0_JWRB
      SUMT = 0.0_JWRB

      DELK_GC(NS) = 0.5*(XK_GC(NS+1)-XK_GC(NS))
      DO I=NS+1,NWAV_GC-1
        DELK_GC(I) = 0.5_JWRB*(XK_GC(I+1)-XK_GC(I-1))
      ENDDO
      DELK_GC(NWAV_GC) = 0.5_JWRB*(XK_GC(NWAV_GC)-XK_GC(NWAV_GC-1))

      DO I = NS, NWAV_GC
        XM    = 1.0_JWRB/XK_GC(I)
!       ANALYTICAL FORM INERTIAL SUB RANGE F1(N) = X**(-4)*BB
        BBDELK    = DELK_GC(I)*SQRT(EPS0THREEZPIM*VG_GC(I))/C_GC(I)**2
        SUMS  = SUMS + BBDELK * XM
        GAM_W = ANG * GAMMA_WAM(OMEGA_GC(I), XK_GC(I), UST, Z0, EPS)
        SUMT  = SUMT + EPSM1 * OMEGA_GC(I) * GAM_W * BBDELK * XM**3
      ENDDO
 
      XMSSCG = SUMS
      TAUWCG = SUMT

      IF (LHOOK) CALL DR_HOOK('STRESS_GC',1,ZHOOK_HANDLE)
 
      END SUBROUTINE STRESS_GC
