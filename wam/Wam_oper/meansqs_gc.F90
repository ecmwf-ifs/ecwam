      SUBROUTINE MEANSQS_GC(ANG_GC, USTAR, Z0, HALPHAP, XMSSCG)

!***  DETERMINE MSS FOR GRAV-CAP WAVES

!     AUTHOR: PETER JANSSEN
!     ------

!     REFERENCES:
!     ----------

!     VIERS PAPER EQ.(29)

!----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : NWAV_GC, KRATIO_GC, OMEGA_GC, XK_GC, XKM_GC, &
     &                      VG_GC, C2OSQRTVG_GC, DELKCC_GC
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
      REAL(KIND=JWRB), INTENT(OUT) :: XMSSCG  ! mean square slope for gravity-capillary waves

      INTEGER(KIND=JWIM) :: NS
      INTEGER(KIND=JWIM) :: I


      REAL(KIND=JWRB) :: XKS, OMS
      REAL(KIND=JWRB) :: COEF
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NWAV_GC) :: BBDELK
   
!     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
#include "gc_dispersion.h"

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('MEANSQS_GC',0,ZHOOK_HANDLE)

!*    1.0  DETERMINE GRAV_CAP SPECTRUM, MSS AND TAUWHF.
!          -------------------------------------------

      CALL OMEGAGC(USTAR, NS, XKS, OMS)

      COEF = C2OSQRTVG_GC(NS)*HALPHAP
      XMSSCG = 0.0_JWRB
      DO I = NS, NWAV_GC
!       ANALYTICAL FORM INERTIAL SUB RANGE F(k) = k**(-4)*BB
!       BB = COEF*SQRT(VG_GC(I))/C_GC(I)**2
        BBDELK(I) = COEF*DELKCC_GC(I)
!       mss :  integral of k**2 F(k)  k dk
        XMSSCG = XMSSCG + BBDELK(I) * XKM_GC(I) 
      ENDDO

      IF (LHOOK) CALL DR_HOOK('MEANSQS_GC',1,ZHOOK_HANDLE)
 
      END SUBROUTINE MEANSQS_GC
