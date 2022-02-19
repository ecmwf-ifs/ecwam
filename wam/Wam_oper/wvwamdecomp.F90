      SUBROUTINE WVWAMDECOMP

! ----------------------------------------------------------------------

!**** *WVWAMDECOMP* - WAVE MODEL MODEL DECOMPOSITION

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : LOUTMDLDCP
      USE YOWGRID  , ONLY : IJS      ,IJL
      USE YOWMPP   , ONLY : NPROC    ,MPMAXLENGTH
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "mpdecomp.intfb.h"
#include "outmdlcdp.intfb.h"

      INTEGER(KIND=JWIM) :: NPR, MAXLEN

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      LOGICAL :: LLIRANK

! ----------------------------------------------------------------------
 
      IF (LHOOK) CALL DR_HOOK('WVWAMDECOMP',0,ZHOOK_HANDLE)

      NPR=NPROC
      LLIRANK=.FALSE.
      CALL MPDECOMP(NPR,MAXLEN,LLIRANK)
      MPMAXLENGTH=MAXLEN

      ! Write model decomposition to a file
      IF ( LOUTMDLDCP ) CALL OUTMDLDCP(IJS(1), IJL(1))

      IF (LHOOK) CALL DR_HOOK('WVWAMDECOMP',1,ZHOOK_HANDLE)
 
      END SUBROUTINE WVWAMDECOMP
