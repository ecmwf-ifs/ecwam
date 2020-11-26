      SUBROUTINE WVWAMDECOMP

! ----------------------------------------------------------------------

!**** *WVWAMDECOMP* - WAVE MODEL MODEL DECOMPOSITION

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMPP   , ONLY : NPROC    ,MPMAXLENGTH
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "mpdecomp.intfb.h"

      INTEGER(KIND=JWIM) :: NPR, MAXLEN

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      LOGICAL :: LLIRANK

! ----------------------------------------------------------------------
 
      IF (LHOOK) CALL DR_HOOK('WVWAMDECOMP',0,ZHOOK_HANDLE)

      NPR=NPROC
      LLIRANK=.FALSE.
      CALL MPDECOMP(NPR,MAXLEN,LLIRANK)
      MPMAXLENGTH=MAXLEN

      IF (LHOOK) CALL DR_HOOK('WVWAMDECOMP',1,ZHOOK_HANDLE)
 
      END SUBROUTINE WVWAMDECOMP
