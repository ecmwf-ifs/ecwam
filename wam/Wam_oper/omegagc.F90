      SUBROUTINE OMEGAGC(UST, NS, XKS, OMS)

!***  DETERMINE THE CUT-OFF ANGULAR FREQUENCY FOR THE GRAV-CAPILLARY WAVES
!     !!!! rounded to the closest index of XK_GC  !!!!!

!     AUTHOR: PETER JANSSEN
!     ------

!     REFERENCES:
!     ----------

!     VIERS PAPER EQ.(29)

!----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : KRATIO_GC, NWAV_GC, OMEGA_GC, XK_GC
      USE YOWPCONS , ONLY : G, SURFT

      USE YOMHOOK  ,ONLY : LHOOK,   DR_HOOK

!----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN) :: UST
      INTEGER(KIND=JWIM), INTENT(OUT) :: NS ! index in array XK_GC corresponding to XKS and OMS
      REAL(KIND=JWRB), INTENT(OUT) :: XKS   ! cut-off wave number
      REAL(KIND=JWRB), INTENT(OUT) :: OMS   ! cut-off angular frequency

      REAL(KIND=JWRB) :: Y
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
   
!     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
#include "gc_dispersion.h"

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('OMEGAGC',0,ZHOOK_HANDLE)

!!!      Y = 1.0_JWRB/(1.48_JWRB+2.05_JWRB*UST)
      Y = (1.0_JWRB + UST**2)/(1.0_JWRB+10.0_JWRB*UST**2)
      XKS = Y*SQRT(G/SURFT)
      NS = MIN(INT(LOG(XKS/XK_GC(1))/LOG(KRATIO_GC))+1,NWAV_GC-1)
      XKS = XK_GC(NS)
      OMS = OMEGA_GC(NS)

      IF (LHOOK) CALL DR_HOOK('OMEGAGC',1,ZHOOK_HANDLE)
 
      END SUBROUTINE OMEGAGC
