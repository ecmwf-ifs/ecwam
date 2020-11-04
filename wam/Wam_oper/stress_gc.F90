      SUBROUTINE STRESS_GC(ANG_GC, USTAR, Z0, HALP, ZBREDUC, TAUWCG)

!***  DETERMINE WAVE INDUCED STRESS FOR GRAV-CAP WAVES

!     AUTHOR: PETER JANSSEN
!     ------

!     REFERENCES:
!     ----------

!     VIERS PAPER EQ.(29)
!     FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

!----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LLNORMAGAM
      USE YOWFRED  , ONLY : NWAV_GC, OMEGA_GC, XK_GC, &
     &                      OMXKM3_GC, CM_GC, C2OSQRTVG_GC, XKMSQRTVGOC2_GC, &
     &                      OM3GMKM_GC, DELKCC_GC, DELKCC_GC_NS
      USE YOWPCONS , ONLY : G, EPSUS, SURFT
      USE YOWPHYS  , ONLY : XKAPPA, ZALP,   BETAMAXOXKAPPA2, DELTA_THETA_RN

      USE YOMHOOK  ,ONLY : LHOOK,   DR_HOOK

!----------------------------------------------------------------------

      IMPLICIT NONE

#include "omegagc.intfb.h"

      REAL(KIND=JWRB), INTENT(IN) :: ANG_GC  ! factor to account for angular spreading of the input.
      REAL(KIND=JWRB), INTENT(IN) :: USTAR ! friction velocity
      REAL(KIND=JWRB), INTENT(IN) :: Z0 !  surface roughness
      REAL(KIND=JWRB), INTENT(IN) :: HALP  ! 1/2 Phillips parameter
      REAL(KIND=JWRB), INTENT(IN) :: ZBREDUC  ! reduction factor for the gravity-capillary stress
      REAL(KIND=JWRB), INTENT(OUT) :: TAUWCG ! wave induced kinematic stress for gravity-capillary waves

      INTEGER(KIND=JWIM) :: NS
      INTEGER(KIND=JWIM) :: I

      REAL, PARAMETER :: TAUWCG_MIN = 0.00001_JWRB
      REAL(KIND=JWRB) :: XKS, OMS, ZABHRC
      REAL(KIND=JWRB) :: X, XLOG, ZLOG, ZLOG2X
      REAL(KIND=JWRB) :: CONST, USTM 
      REAL(KIND=JWRB) :: GAMNORMA ! RENORMALISATION FACTOR OF THE GROWTH RATE
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

      ZABHRC = ANG_GC * BETAMAXOXKAPPA2 * HALP * ZBREDUC * C2OSQRTVG_GC(NS)
      IF(LLNORMAGAM) THEN
!!        CONST = DELTA_THETA_RN * BETAMAXOXKAPPA2 * HALP * ZBREDUC * C2OSQRTVG_GC(NS) /XKAPPA
        CONST = ZABHRC /XKAPPA
      ELSE
        CONST = 0.0_JWRB
      ENDIF

      USTM = 1.0_JWRB / MAX(USTAR,EPSUS)

      DO I = NS, NWAV_GC
!       GROWTHRATE BY WIND WITHOUT the multiplicative representing the ratio of air density to water density (eps)
!       and BETAMAXOXKAPPA2
        X       = USTAR*CM_GC(I)
        XLOG    = LOG(XK_GC(I)*Z0) + XKAPPA/(X + ZALP) 
        ZLOG    = MIN(XLOG,0.0_JWRB)
        ZLOG2X  = ZLOG*ZLOG*X
        GAM_W(I)= ZLOG2X*EXP(ZLOG)*ZLOG2X*OM3GMKM_GC(I)
      ENDDO

      GAMNORMA  = 1.0_JWRB + CONST*XKMSQRTVGOC2_GC(NS)*GAM_W(NS)*USTM
      TAUWCG = GAM_W(NS) * DELKCC_GC_NS(NS) * OMXKM3_GC(NS) / GAMNORMA
      DO I = NS+1, NWAV_GC
!       ANALYTICAL FORM INERTIAL SUB RANGE F(k) = k**(-4)*BB
!       BB = HALP * C2OSQRTVG_GC(NS)*SQRT(VG_GC(I))/C_GC(I)**2
!       Tauwcg : (rhow * g /rhoa) * integral of (1/c) * gammma * F(k)  k dk 
!       with omega=g*k and omega=k*c,  then
!       Tauwcg : (rhow /rhoa) * integral of omega * gammma * F(k)  k dk
!       but gamma is computed wihtout the rhoa/rhow factor so
!       Tauwcg : integral of omega * gammma_wam * F(k)  k dk
!       It should be done in vector form with actual directional spreading information
!       It simplified here by using the ANG_GC factor.
        GAMNORMA  = 1.0_JWRB + CONST*XKMSQRTVGOC2_GC(I)*GAM_W(I)*USTM
        TAUWCG = TAUWCG + GAM_W(I) * DELKCC_GC(I) * OMXKM3_GC(I) / GAMNORMA
      ENDDO
      TAUWCG = MAX(ZABHRC * TAUWCG, TAUWCG_MIN)

      IF (LHOOK) CALL DR_HOOK('STRESS_GC',1,ZHOOK_HANDLE)
 
      END SUBROUTINE STRESS_GC
