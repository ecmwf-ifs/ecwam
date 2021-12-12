SUBROUTINE OMEGAGC(KIJS, KIJL, UST, NS, XKS, OMS)

!***  DETERMINE THE CUT-OFF ANGULAR FREQUENCY FOR THE GRAV-CAPILLARY WAVES
!     !!!! rounded to the closest index of XK_GC  !!!!!

!     AUTHOR: PETER JANSSEN
!     ------

!     REFERENCES:
!     ----------

!     VIERS PAPER EQ.(29)

!----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWFRED  , ONLY : OMEGA_GC, XK_GC

USE YOMHOOK  ,ONLY : LHOOK,   DR_HOOK

!----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: UST
INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(OUT) :: NS ! index in array XK_GC corresponding to XKS and OMS
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: XKS   ! cut-off wave number
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: OMS   ! cut-off angular frequency


INTEGER(KIND=JWIM) :: IJ
REAL(KIND=JWRB) :: ZHOOK_HANDLE
   
#include "ns_gc.intfb.h"

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('OMEGAGC',0,ZHOOK_HANDLE)

DO IJ = KIJS, KIJL
  NS(IJ) = NS_GC(UST(IJ))
  XKS(IJ) = XK_GC(NS(IJ))
  OMS(IJ) = OMEGA_GC(NS(IJ))
ENDDO

IF (LHOOK) CALL DR_HOOK('OMEGAGC',1,ZHOOK_HANDLE)
 
END SUBROUTINE OMEGAGC
