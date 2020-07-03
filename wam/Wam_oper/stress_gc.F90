      SUBROUTINE STRESS_GC(ANG_GC, USTAR, Z0, HALPHAP, TAUWCG)

!***  DETERMINE WAVE INDUCED STRESS FOR GRAV-CAP WAVES

!     AUTHOR: PETER JANSSEN
!     ------

!     REFERENCES:
!     ----------

!     VIERS PAPER EQ.(29)

!----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : NWAV_GC, KRATIO_GC, OMEGA_GC, XK_GC, XKM_GC, &
     &                      OMXKM3_GC, VG_GC, C_GC, C2OSQRTVG_GC, OM3GMKM_GC, &
     &                      DELKCC_GC
      USE YOWPCONS , ONLY : G,      SURFT

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

      REAL(KIND=JWRB) :: GAMMA_WAM

      REAL(KIND=JWRB) :: XKS, OMS
      REAL(KIND=JWRB) :: COEF
      REAL(KIND=JWRB) :: BS, OM
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NWAV_GC) :: GAM_W, BBDELK
   
!     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
#include "gc_dispersion.h"

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('STRESS_GC',0,ZHOOK_HANDLE)

!*    1.0  DETERMINE GRAV_CAP SPECTRUM, TAUWHF.
!          ------------------------------------

      CALL OMEGAGC(USTAR, NS, XKS, OMS)

      DO I = NS, NWAV_GC
        GAM_W(I) = ANG_GC * GAMMA_WAM(XK_GC(I), C_GC(I), OM3GMKM_GC(I), USTAR, Z0)
      ENDDO

      COEF = C2OSQRTVG_GC(NS)*HALPHAP
      TAUWCG = 0.0_JWRB
      DO I = NS, NWAV_GC
!       ANALYTICAL FORM INERTIAL SUB RANGE F(k) = k**(-4)*BB
!       BB = COEF*SQRT(VG_GC(I))/C_GC(I)**2
        BBDELK(I) = COEF*DELKCC_GC(I)
!       Tauwcg : (rhow * g /rhoa) * integral of (1/c) * gammma * F(k)  k dk 
!       with omega=g*k and omega=k*c,  then
!       Tauwcg : (rhow /rhoa) * integral of omega * gammma * F(k)  k dk
!       but gamma is computed wihtout the rhoa/rhow factor so
!       Tauwcg : integral of omega * gammma_wam * F(k)  k dk
!       It should be done in vector form with actual directional spreading information
!       It simplfied here by using the ANG_GC factor.
        TAUWCG = TAUWCG + GAM_W(I) * BBDELK(I) * OMXKM3_GC(I) 
      ENDDO

      IF (LHOOK) CALL DR_HOOK('STRESS_GC',1,ZHOOK_HANDLE)
 
      END SUBROUTINE STRESS_GC
