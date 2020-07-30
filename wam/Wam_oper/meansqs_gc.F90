      SUBROUTINE MEANSQS_GC(IJS, IJL, HALPHAP, USTAR, XMSSCG, FRGC)

!***  DETERMINE MSS FOR GRAV-CAP WAVES

!     AUTHOR: PETER JANSSEN
!     ------

!     REFERENCES:
!     ----------

!     VIERS PAPER EQ.(29)

!----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : NWAV_GC, KRATIO_GC, OMEGA_GC, XK_GC, XKM_GC, &
     &                      VG_GC, C2OSQRTVG_GC, DELKCC_GC, DELKCC_GC_NS
      USE YOWPCONS , ONLY : G, ZPI,  SURFT

      USE YOMHOOK  ,ONLY : LHOOK,   DR_HOOK

!----------------------------------------------------------------------

      IMPLICIT NONE

#include "omegagc.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL

      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: HALPHAP  ! 1/2 Phillips parameter
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: USTAR ! friction velocity
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: XMSSCG  ! mean square slope for gravity-capillary waves
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: FRGC  ! Frequency from which the gravity-capillary spectrum is approximated

      INTEGER(KIND=JWIM) :: IJ, I, NS
      REAL(KIND=JWRB) :: XKS, OMS, COEF
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
   
!     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
#include "gc_dispersion.h"

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('MEANSQS_GC',0,ZHOOK_HANDLE)

      DO IJ = IJS, IJL
        CALL OMEGAGC(USTAR(IJ), NS, XKS, OMS)
        FRGC(IJ) = OMS/ZPI
        COEF = C2OSQRTVG_GC(NS)*HALPHAP(IJ)
        XMSSCG(IJ) = DELKCC_GC_NS(NS) * XKM_GC(NS) 
        DO I = NS+1, NWAV_GC
!         ANALYTICAL FORM INERTIAL SUB RANGE F(k) = k**(-4)*BB
!         BB = COEF*SQRT(VG_GC(I))/C_GC(I)**2
!         mss :  integral of k**2 F(k)  k dk
          XMSSCG(IJ) = XMSSCG(IJ) + DELKCC_GC(I) * XKM_GC(I) 
        ENDDO
        XMSSCG(IJ) = XMSSCG(IJ)*COEF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('MEANSQS_GC',1,ZHOOK_HANDLE)
 
      END SUBROUTINE MEANSQS_GC
