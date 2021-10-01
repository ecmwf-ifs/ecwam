      SUBROUTINE WVALLOC

! ----------------------------------------------------------------------

!**** *WVALLOC* - WAVE MODEL ARRAY ALLOCATION

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWMEAN  , ONLY : INTFLDS
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSHAL  , ONLY : WVPRPT
      USE YOWSPEC  , ONLY : NSTART   ,NEND     ,FF_NOW   ,FL1
      USE YOWSTAT  , ONLY : IPROPAGS ,LSUBGRID ,IREFRA   ,IDELPRO
      USE YOWTEST  , ONLY : IU06
      USE YOWWIND  , ONLY : FF_NEXT

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: NPR, MAXLEN

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
 
      IF (LHOOK) CALL DR_HOOK('WVALLOC',0,ZHOOK_HANDLE)

!     1.  ALLOCATE NECESSARY ARRAYS
!         -------------------------

      IF (.NOT.ALLOCATED(WVPRPT)) ALLOCATE(WVPRPT(NSTART(IRANK):NEND(IRANK),NFRE))

      IF (.NOT.ALLOCATED(FF_NOW)) ALLOCATE(FF_NOW(NSTART(IRANK):NEND(IRANK)))


      IF (.NOT.ALLOCATED(FL1)) THEN
        ALLOCATE(FL1(NSTART(IRANK):NEND(IRANK),NANG,NFRE))
        FL1(:,:,:) = 0.0_JWRB
      ENDIF

      IF (.NOT.ALLOCATED(INTFLDS)) THEN 
        ALLOCATE(INTFLDS(NSTART(IRANK):NEND(IRANK)))

        INTFLDS(:)%PHIEPS  = 0.0_JWRB
        INTFLDS(:)%PHIAW   = 0.0_JWRB
        INTFLDS(:)%TAUOC   = 0.0_JWRB
        INTFLDS(:)%STRNMS  = 0.0_JWRB
        INTFLDS(:)%ALTWH   = ZMISS 
        INTFLDS(:)%CALTWH  = ZMISS 
        INTFLDS(:)%RALTCOR = ZMISS 
      ENDIF


      IF (.NOT.ALLOCATED(FF_NEXT)) ALLOCATE(FF_NEXT(NSTART(IRANK):NEND(IRANK)))


      IF (LHOOK) CALL DR_HOOK('WVALLOC',1,ZHOOK_HANDLE)
 
      END SUBROUTINE WVALLOC
