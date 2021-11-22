      SUBROUTINE WVWAMINIT1(LDWCOUIFS,LDWCOU2W,LDWFLUX,LWCUR,LFDBOPIFS)

! ----------------------------------------------------------------------

!**** *WVWAMINIT1* - WAVE MODEL CONTINUED INITIALISATION

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU    ,LWCOU2W  ,LWFLUX
      USE YOWCOUT  , ONLY : LFDB
      USE YOWMESPAS, ONLY : LFDBIOOUT
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPREVIOUS,NNEXT    ,    &
     &         MPMAXLENGTH
      USE MPL_MODULE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal, MyRankLocal

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      LOGICAL :: LFDBOPIFS, LDWCOUIFS, LDWCOU2W, LDWFLUX
      LOGICAL :: LWCUR

! ----------------------------------------------------------------------
 
      IF (LHOOK) CALL DR_HOOK('WVWAMINIT1',0,ZHOOK_HANDLE)

! RE-INITIALIZE LOGICALS IF COUPLED TO IFS

      IF(LDWCOUIFS) THEN
        LWCOU=LDWCOUIFS
        LWCOU2W=LDWCOU2W
        LWFLUX=LDWFLUX
        LFDB=LFDBOPIFS
        LFDBIOOUT=LFDBOPIFS
      ENDIF

#if !defined MODEL_COUPLING_ATM_WAV && !defined MODEL_COUPLING_OCN_WAV
        MyRankGlobal=IRANK-1
        MyRankLocal =IRANK-1
#endif

      NPREVIOUS=IRANK-1
      IF(IRANK.EQ.NPROC) THEN
        NNEXT=0
      ELSE
        NNEXT=IRANK+1
      ENDIF

      IF (LHOOK) CALL DR_HOOK('WVWAMINIT1',1,ZHOOK_HANDLE)
 
      END SUBROUTINE WVWAMINIT1
