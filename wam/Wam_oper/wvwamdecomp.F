      SUBROUTINE WVWAMDECOMP

! ----------------------------------------------------------------------

!**** *WVWAMDECOMP* - WAVE MODEL MODEL DECOMPOSITION

! ----------------------------------------------------------------------

      USE YOWMPP   , ONLY : NPROC    ,MPMAXLENGTH
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "mpdecomp.intfb.h"

      INTEGER :: NPR, MAXLEN

      REAL :: ZHOOK_HANDLE

      LOGICAL :: LLIRANK

! ----------------------------------------------------------------------
#ifdef ECMWF 
      IF (LHOOK) CALL DR_HOOK('WVWAMDECOMP',0,ZHOOK_HANDLE)
#endif

      NPR=NPROC
      LLIRANK=.FALSE.
      CALL MPDECOMP(NPR,MAXLEN,LLIRANK)
      MPMAXLENGTH=MAXLEN

#ifdef ECMWF 
      IF (LHOOK) CALL DR_HOOK('WVWAMDECOMP',1,ZHOOK_HANDLE)
#endif
 
      END SUBROUTINE WVWAMDECOMP
