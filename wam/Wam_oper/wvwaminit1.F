      SUBROUTINE WVWAMINIT1(LDWCOUIFS,LDWCOU2W,LDWFLUX,LWCUR,LFDBOPIFS)

! ----------------------------------------------------------------------

!**** *WVWAMINIT1* - WAVE MODEL CONTINUED INITIALISATION

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU    ,LWCOU2W  ,LWFLUX
      USE YOWCOUT  , ONLY : LFDB
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPREVIOUS,NNEXT    ,
     &         MPMAXLENGTH
      USE MPL_MODULE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL :: ZHOOK_HANDLE

      LOGICAL :: LFDBOPIFS,LDWCOUIFS, LDWCOU2W, LDWFLUX
      LOGICAL :: LWCUR

! ----------------------------------------------------------------------
#ifdef ECMWF 
      IF (LHOOK) CALL DR_HOOK('WVWAMINIT1',0,ZHOOK_HANDLE)
#endif

! RE-INITIALIZE LOGICALS IF COUPLED TO IFS

      IF(LDWCOUIFS) THEN
        LFDB=LFDBOPIFS
        LWCOU=LDWCOUIFS
        LWCOU2W=LDWCOU2W
        LWFLUX=LDWFLUX
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

#ifdef ECMWF 
      IF (LHOOK) CALL DR_HOOK('WVWAMINIT1',1,ZHOOK_HANDLE)
#endif
 
      END SUBROUTINE WVWAMINIT1
