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

      USE YOWTEST  , ONLY : IU06

!----------------------------------------------------------------------

      IMPLICIT NONE

#include "omegagc.intfb.h"

      REAL(KIND=JWRB), INTENT(IN) :: UST ! friction velocity
      REAL(KIND=JWRB), INTENT(IN) :: EPS ! ration of air density to water density
      REAL(KIND=JWRB), INTENT(IN) :: Z0 !  surface roughness
      REAL(KIND=JWRB), INTENT(IN) :: ALPHAP  ! Phillips parameter
      REAL(KIND=JWRB), INTENT(OUT) :: XMSSCG  ! mean squre slope for gravity-capillary waves
      REAL(KIND=JWRB), INTENT(OUT) :: TAUWCG ! wave induced kinematic stress for gravity-capillary waves

      INTEGER(KIND=JWIM) :: NS
      INTEGER(KIND=JWIM) :: I

      REAL(KIND=JWRB), PARAMETER :: ANG = 0.5_JWRB   ! factor to account for angular spreading.

      REAL(KIND=JWRB) :: GAMMA_WAM

      REAL(KIND=JWRB) :: XKS, OMS
      REAL(KIND=JWRB) :: XM, BBDELK, EPS0THREEZPIM
      REAL(KIND=JWRB) :: SUMT, SUMS, GAM_W, BS, OM
      REAL(KIND=JWRB) :: DELK_GC(NWAV_GC)
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
   
!     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
#include "gc_dispersion.h"

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('STRESS_GC',0,ZHOOK_HANDLE)


!*    1.0  DETERMINE GRAV_CAP SPECTRUM, MSS AND TAUWHF.
!          -------------------------------------------

      CALL OMEGAGC(UST, NS, XKS, OMS)

!!!! debile 
        write(IU06,*) 'OMS in stress_gc ',OMS



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
!       ANALYTICAL FORM INERTIAL SUB RANGE F(k) = k**(-4)*BB
!        BB = SQRT(EPS0THREEZPIM*VG_GC(I))/C_GC(I)**2
        BBDELK    = DELK_GC(I)*SQRT(EPS0THREEZPIM*VG_GC(I))/C_GC(I)**2
!       mss :  integral of k**2 F(k)  k dk
        SUMS  = SUMS + BBDELK * XM
!       Tauwcg : (rhow * g /rhoa) * integral of (1/c) * gammma * F(k)  k dk 
!       with omega=g*k and omega=k*c,  then
!       Tauwcg : (rhow /rhoa) * integral of omega * gammma * F(k)  k dk
!       It should be done in vector form with actual directional spreading information
!       It simplfied here by using the ANG factor.
        GAM_W = ANG * GAMMA_WAM(OMEGA_GC(I), XK_GC(I), UST, Z0, EPS)
!       the factor 1/eps is applied after the loop.
        SUMT  = SUMT + OMEGA_GC(I) * GAM_W * BBDELK * XM**3
      ENDDO
 
      XMSSCG = SUMS
      TAUWCG = SUMT/EPS

      IF (LHOOK) CALL DR_HOOK('STRESS_GC',1,ZHOOK_HANDLE)
 
      END SUBROUTINE STRESS_GC
