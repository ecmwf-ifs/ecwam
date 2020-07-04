      SUBROUTINE STRESS_GC(ANG_GC, USTAR, Z0, HALPHAP, TAUWCG)

!***  DETERMINE WAVE INDUCED STRESS FOR GRAV-CAP WAVES

!     AUTHOR: PETER JANSSEN
!     ------

!     REFERENCES:
!     ----------

!     VIERS PAPER EQ.(29)
!     FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

!----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : NWAV_GC, OMEGA_GC, XK_GC, &
     &                      OMXKM3_GC, CM_GC, C2OSQRTVG_GC, OM3GMKM_GC, &
     &                      DELKCC_GC
      USE YOWPCONS , ONLY : G,      SURFT
      USE YOWPHYS  , ONLY : XKAPPA, ZALP,   BETAMAXOXKAPPA2

      USE YOMHOOK  ,ONLY : LHOOK,   DR_HOOK

      USE YOWTEST  , ONLY : IU06

!----------------------------------------------------------------------

      IMPLICIT NONE

#include "omegagc.intfb.h"

      REAL(KIND=JWRB), INTENT(IN) :: ANG_GC  ! factor to account for angular spreading of the input.
      REAL(KIND=JWRB), INTENT(IN) :: USTAR ! friction velocity
      REAL(KIND=JWRB), INTENT(IN) :: Z0 !  surface roughness
      REAL(KIND=JWRB), INTENT(IN) :: HALPHAP  ! 1/2 Phillips parameter
      REAL(KIND=JWRB), INTENT(OUT) :: TAUWCG ! wave induced kinematic stress for gravity-capillary waves

      INTEGER(KIND=JWIM) :: NS
      INTEGER(KIND=JWIM) :: I

      REAL, PARAMETER :: TAUWCG_MIN = 0.00001_JWRB
      REAL(KIND=JWRB) :: XKS, OMS
      REAL(KIND=JWRB) :: X, XLOG, ZLOG, ZLOG2X
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NWAV_GC) :: GAM_W
   
!     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
#include "gc_dispersion.h"

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('STRESS_GC',0,ZHOOK_HANDLE)

!*    1.0  DETERMINE GRAV_CAP SPECTRUM, TAUWHF.
!          ------------------------------------

!     FIND NS:
      CALL OMEGAGC(USTAR, NS, XKS, OMS)

      DO I = NS, NWAV_GC
!       GROWTHRATE BY WIND WITHOUT the multiplicative representing the ratio of air density to water density (eps)
!       and BETAMAXOXKAPPA2
        X       = USTAR*CM_GC(I)
        XLOG    = LOG(XK_GC(I)*Z0) + XKAPPA/(X + ZALP) 
        ZLOG    = MIN(XLOG,0.0_JWRB)
        ZLOG2X  = ZLOG*ZLOG*X
        GAM_W(I)= ZLOG2X*EXP(ZLOG)*ZLOG2X*OM3GMKM_GC(I)
      ENDDO

      TAUWCG = 0.0_JWRB
      DO I = NS, NWAV_GC
!       ANALYTICAL FORM INERTIAL SUB RANGE F(k) = k**(-4)*BB
!       BB = HALPHAP * C2OSQRTVG_GC(NS)*SQRT(VG_GC(I))/C_GC(I)**2
!       Tauwcg : (rhow * g /rhoa) * integral of (1/c) * gammma * F(k)  k dk 
!       with omega=g*k and omega=k*c,  then
!       Tauwcg : (rhow /rhoa) * integral of omega * gammma * F(k)  k dk
!       but gamma is computed wihtout the rhoa/rhow factor so
!       Tauwcg : integral of omega * gammma_wam * F(k)  k dk
!       It should be done in vector form with actual directional spreading information
!       It simplified here by using the ANG_GC factor.
        TAUWCG = TAUWCG + GAM_W(I) * DELKCC_GC(I) * OMXKM3_GC(I) 
      ENDDO
      TAUWCG = MAX(ANG_GC * BETAMAXOXKAPPA2 * HALPHAP * C2OSQRTVG_GC(NS) * TAUWCG, TAUWCG_MIN)

      IF (LHOOK) CALL DR_HOOK('STRESS_GC',1,ZHOOK_HANDLE)
 
      END SUBROUTINE STRESS_GC
