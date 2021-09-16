      SUBROUTINE WVWAMINIT1(LDWCOUIFS, LDWCOU2W, LDWCOURNW, LDWCOUHMF, LDWFLUX, LFDBOPIFS)

! ----------------------------------------------------------------------

!**** *WVWAMINIT1* - WAVE MODEL CONTINUED INITIALISATION

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU    ,LWCOU2W  ,LWCOURNW, LWCOUHMF, LWFLUX
      USE YOWCOUT  , ONLY : LFDB
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPREVIOUS,    NNEXT

      USE MPL_MODULE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      LOGICAL, INTENT(IN) :: LDWCOUIFS, LDWCOU2W, LDWCOURNW, LDWCOUHMF, LDWFLUX, LFDBOPIFS

! ----------------------------------------------------------------------
 
      IF (LHOOK) CALL DR_HOOK('WVWAMINIT1',0,ZHOOK_HANDLE)

! RE-INITIALIZE LOGICALS IF COUPLED TO IFS

      IF (LDWCOUIFS) THEN
        LWCOU=LDWCOUIFS
        LWCOU2W=LDWCOU2W
        LWCOURNW=LDWCOURNW
        LWCOUHMF=LDWCOUHMF
        LWFLUX=LDWFLUX
        LFDB=LFDBOPIFS
      ENDIF

      NPREVIOUS=IRANK-1
      IF (IRANK == NPROC) THEN
        NNEXT=0
      ELSE
        NNEXT=IRANK+1
      ENDIF

      IF (LHOOK) CALL DR_HOOK('WVWAMINIT1',1,ZHOOK_HANDLE)
 
      END SUBROUTINE WVWAMINIT1
