      SUBROUTINE MEANSQS_GC(XKMSS, IJS, IJL, HALPHAP, U10, USTAR, XMSSCG, FRGC)

!***  DETERMINE MSS FOR GRAV-CAP WAVES UP TO WAVE NUMBER XKMSS

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

      REAL(KIND=JWRB), INTENT(IN) :: XKMSS ! WAVE NUMBER CUT-OFF
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: HALPHAP  ! 1/2 Phillips parameter
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: U10 ! 10m wind speed
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: USTAR ! friction velocity
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: XMSSCG  ! mean square slope for gravity-capillary waves
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: FRGC  ! Frequency from which the gravity-capillary spectrum is approximated

      INTEGER(KIND=JWIM) :: IJ, I, NE
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL) :: NS
      REAL(KIND=JWRB) :: DIRSPRD_GC
      REAL(KIND=JWRB) :: XKS, OMS, COEF
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
   
!     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
#include "gc_dispersion.h"

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('MEANSQS_GC',0,ZHOOK_HANDLE)

      NE = MIN(MAX(NINT(LOG(XKMSS/XK_GC(1))/(LOG(KRATIO_GC))),1),NWAV_GC)

      DO IJ = IJS, IJL
        CALL OMEGAGC(USTAR(IJ), NS(IJ), XKS, OMS)
        FRGC(IJ) = OMS/ZPI
        IF(XKS > XKMSS) THEN
          NS(IJ) = NE
          XMSSCG(IJ) = 0.0_JWRB
        ELSE
          XMSSCG(IJ) = DELKCC_GC_NS(NS(IJ)) * XKM_GC(NS(IJ)) 
        ENDIF
      ENDDO

      DO IJ = IJS, IJL
        DO I = NS(IJ)+1, NE 
!         ANALYTICAL FORM INERTIAL SUB RANGE F(k) = k**(-4)*BB
!         BB = COEF*SQRT(VG_GC(I))/C_GC(I)**2
!         mss :  integral of k**2 F(k)  k dk
          XMSSCG(IJ) = XMSSCG(IJ) + DELKCC_GC(I) * XKM_GC(I) 
        ENDDO
        COEF = C2OSQRTVG_GC(NS(IJ))*HALPHAP(IJ)*DIRSPRD_GC(U10(IJ))
        XMSSCG(IJ) = XMSSCG(IJ)*COEF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('MEANSQS_GC',1,ZHOOK_HANDLE)
 
      END SUBROUTINE MEANSQS_GC
