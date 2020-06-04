      SUBROUTINE STRESS_GC(USTAR, Z0, ALPHAP, XMSSCG, TAUWCG)

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
      USE YOWPCONS , ONLY : G,      SURFT
      USE YOWPHYS  , ONLY : TAUWSHELTER,  ANG_GC

      USE YOMHOOK  ,ONLY : LHOOK,   DR_HOOK

      USE YOWTEST  , ONLY : IU06

!----------------------------------------------------------------------

      IMPLICIT NONE

#include "omegagc.intfb.h"

      REAL(KIND=JWRB), INTENT(IN) :: USTAR ! friction velocity
      REAL(KIND=JWRB), INTENT(IN) :: Z0 !  surface roughness
      REAL(KIND=JWRB), INTENT(IN) :: ALPHAP  ! Phillips parameter
      REAL(KIND=JWRB), INTENT(OUT) :: XMSSCG  ! mean squre slope for gravity-capillary waves
      REAL(KIND=JWRB), INTENT(OUT) :: TAUWCG ! wave induced kinematic stress for gravity-capillary waves

      INTEGER(KIND=JWIM) :: NS
      INTEGER(KIND=JWIM) :: I

      REAL(KIND=JWRB) :: GAMMA_WAM

      REAL(KIND=JWRB) :: TAU, TAUCT, UST, XKS, OMS
      REAL(KIND=JWRB) :: XM, BBDELK, COEF
      REAL(KIND=JWRB) :: GAM_W, BS, OM
      REAL(KIND=JWRB) :: DELK_GC(NWAV_GC)
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
   
!     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
#include "gc_dispersion.h"

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('STRESS_GC',0,ZHOOK_HANDLE)


!*    1.0  DETERMINE GRAV_CAP SPECTRUM, MSS AND TAUWHF.
!          -------------------------------------------

      CALL OMEGAGC(USTAR, NS, XKS, OMS)

      BS  = 0.5_JWRB*ALPHAP
      COEF = FC_GC(XKS)**4/FVG_GC(XKS)*BS**2
      XMSSCG = 0.0_JWRB
      TAUWCG = 0.0_JWRB

      DELK_GC(NS) = 0.5*(XK_GC(NS+1)-XK_GC(NS))
      DO I=NS+1,NWAV_GC-1
        DELK_GC(I) = 0.5_JWRB*(XK_GC(I+1)-XK_GC(I-1))
      ENDDO
      DELK_GC(NWAV_GC) = 0.5_JWRB*(XK_GC(NWAV_GC)-XK_GC(NWAV_GC-1))

      UST = USTAR
      TAU = USTAR**2
      DO I = NS, NWAV_GC
        XM    = 1.0_JWRB/XK_GC(I)
!       ANALYTICAL FORM INERTIAL SUB RANGE F(k) = k**(-4)*BB
!        BB = SQRT(COEF*VG_GC(I))/C_GC(I)**2
        BBDELK    = DELK_GC(I)*SQRT(COEF*VG_GC(I))/C_GC(I)**2
!       mss :  integral of k**2 F(k)  k dk
        XMSSCG  = XMSSCG + BBDELK * XM
!       Tauwcg : (rhow * g /rhoa) * integral of (1/c) * gammma * F(k)  k dk 
!       with omega=g*k and omega=k*c,  then
!       Tauwcg : (rhow /rhoa) * integral of omega * gammma * F(k)  k dk
!       but gamma is computed wihtout the rhoa/rhow factor so
!       Tauwcg : integral of omega * gammma_wam * F(k)  k dk
!       It should be done in vector form with actual directional spreading information
!       It simplfied here by using the ANG_GC factor.
        GAM_W = ANG_GC * GAMMA_WAM(OMEGA_GC(I), XK_GC(I), UST, Z0)
        TAUCT = OMEGA_GC(I) * GAM_W * BBDELK * XM**3
        TAU = MAX(TAU - TAUWSHELTER*TAUCT, 0.0_JWRB)
        UST = SQRT(TAU)
        TAUWCG  = TAUWCG + TAUCT 
      ENDDO
 
      IF (LHOOK) CALL DR_HOOK('STRESS_GC',1,ZHOOK_HANDLE)
 
      END SUBROUTINE STRESS_GC
