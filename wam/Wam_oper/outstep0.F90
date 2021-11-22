SUBROUTINE OUTSTEP0 (WVENVI, WVPRPT, FF_NOW, INTFLDS,  &
 &                   WAM2NEMO, NEMO2WAM, FL1)

! ----------------------------------------------------------------------

!**** *OUTSTEP0* - OUTPUTS INITIAL CONDITION AND/OR FORECAST STEP 0

!**   INTERFACE.
!     ----------

!     *CALL* *OUTSTEP0 (WVENVI, WVPRPT, FF_NOW, INTFLDS,
!    &                  WAM2NEMO, NEMO2WAM, FL1)
!        *WVENVI*    WAVE ENVIRONMENT FIELDS
!        *WVPRPT*    WAVE PROPERTIES FIELDS
!        *FF_NOW*    FORCING FIELDS AT CURRENT TIME.
!        *INTFLDS*   INTEGRATED/DERIVED PARAMETERS
!        *WAM2NEMO*  WAVE FIELDS PASSED TO NEMO
!        *NEMO2WAM*  FIELDS FRON OCEAN MODEL to WAM

! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY, FORCING_FIELDS,  &
     &                         INTGT_PARAM_FIELDS, WAVE2OCEAN, OCEAN2WAVE

      USE YOWCOUT  , ONLY : JPPFLAG  ,FFLAG    ,GFLAG    ,              & 
     &                      IRCD     ,IRU10    , IRALTHS  ,IRALTHSC  ,IRALTRC ,   &
     &                      NGOUT    ,LLOUTERS ,                                  &
     &                      NIPRMOUT ,                                            &
     &                      LFDB     ,                                            &
     &                      LRSTST0  ,LWAMANOUT
      USE YOWGRIBHD, ONLY : LGRHDIFS 
      USE YOWGRID  , ONLY : IJS      ,IJL       ,IJSLOC   ,IJLLOC  ,NPROMA_WAM, NCHNK
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE
      USE YOWMESPAS, ONLY : LGRIBOUT ,LNOCDIN  ,LWAVEWIND 
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : CDTPRO   ,CDTINTT  ,IREST   , MARSTYPE ,    &
     &                      LLSOURCE , LANAONLY ,LFRSTFLD
      USE YOWTEXT  , ONLY : LRESTARTED
      USE UNWAM, ONLY : EXCHANGE_FOR_FL1
      USE YOWUNPOOL ,ONLY : LLUNSTR

      USE FDBSUBS_MOD, ONLY : IFLUSHFDBSUBS
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
      
! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "mpcrtbl.intfb.h"
#include "outwint.intfb.h"
#include "outwpsp.intfb.h"
#include "outbs.intfb.h"
#include "outint.intfb.h"
#include "outspec.intfb.h"
#include "setice.intfb.h"
#include "wdfluxes.intfb.h"

      TYPE(ENVIRONMENT), DIMENSION(IJS:IJL), INTENT(INOUT) :: WVENVI
      TYPE(FREQUENCY), DIMENSION(IJS:IJL,NFRE), INTENT(INOUT) :: WVPRPT
      TYPE(FORCING_FIELDS), DIMENSION(IJS:IJL), INTENT(INOUT) :: FF_NOW
      TYPE(INTGT_PARAM_FIELDS), DIMENSION(IJS:IJL), INTENT(INOUT) :: INTFLDS
      TYPE(WAVE2OCEAN), DIMENSION(IJS:IJL), INTENT(INOUT) :: WAM2NEMO
      TYPE(OCEAN2WAVE), DIMENSION(IJS:IJL), INTENT(IN) :: NEMO2WAM
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL1


      INTEGER(KIND=JWIM) :: IJ
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL) :: MIJ

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:) :: BOUTST0
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: XLLWS

      CHARACTER(LEN= 2) :: MARSTYPEBAK
      CHARACTER(LEN=14) :: CDTINTTBAK

      LOGICAL :: LLFLUSH
      LOGICAL :: FFLAGBAK(JPPFLAG), GFLAGBAK(JPPFLAG)

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('OUTSTEP0',0,ZHOOK_HANDLE)

ASSOCIATE(WSWAVE => FF_NOW%WSWAVE, &
 &        WDWAVE => FF_NOW%WDWAVE, &
 &        CICOVER => FF_NOW%CICOVER)


      LLFLUSH = .FALSE.
      LRSTST0=.FALSE.

      NPROMA=NPROMA_WAM


      CDTINTTBAK=CDTINTT
      CDTINTT=CDTPRO


!     1.0  WAVES FLUXES AT INITIAL TIME
!          ----------------------------
      IF (LLSOURCE) THEN

!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          CALL WDFLUXES (KIJS, KIJL,                                 &
     &                   MIJ(KIJS),                                  &
     &                   FL1(KIJS:KIJL,:,:), XLLWS(KIJS:KIJL,:,:),   &
     &                   WVPRPT(KIJS:KIJL,:),                        &
     &                   WVENVI(KIJS), FF_NOW(KIJS),                 &
     &                   INTFLDS(KIJS), WAM2NEMO(KIJS) )
        ENDDO
!$OMP   END PARALLEL DO

        IF (LLUNSTR) CALL EXCHANGE_FOR_FL1(FL1)

      ELSE
        MIJ(:) = NFRE
        XLLWS(:,:,:) = 0.0_JWRB
      ENDIF


!     1.2 SET SPECTRA ON ICE POINTS TO ZERO
!         ---------------------------------
      IF (LICERUN .AND. LMASKICE .AND. LLSOURCE) THEN
        CALL GSTATS(1439,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          CALL SETICE(KIJS, KIJL, FL1(KIJS:KIJL,:,:) ,            &
     &                CICOVER(KIJS:KIJL), WSWAVE(KIJS:KIJL), WDWAVE(KIJS:KIJL))
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1439,1)
      ENDIF

!     2.0 OUTPUT POINT SPECTRA (not usually used at ECMWF)
!         -----------------------------------------------
      IF (NGOUT > 0 .OR. LLOUTERS) CALL OUTWPSP (IJS, IJL, FL1, FF_NOW)


!*    3.0 SAVE INITIAL INTEGRATED FIELDS (if needed)
!         ------------------------------
      IF ((MARSTYPE == 'cf' .OR. MARSTYPE == 'pf' .OR.                 &
     &     MARSTYPE == 'fc' .OR. MARSTYPE == '4v' .OR.                 &
     &     LANAONLY         .OR. LFRSTFLD             )                &
     &   .AND. LWAMANOUT) THEN

        IF (LFRSTFLD) THEN
          FFLAGBAK=FFLAG
          GFLAGBAK=GFLAG
          MARSTYPEBAK=MARSTYPE
          IF (MARSTYPE == 'fg'.OR.MARSTYPE == 'an') THEN
            IF (.NOT.LNOCDIN) THEN
              FFLAG(IRCD)=.FALSE.
              GFLAG(IRCD)=.FALSE.
            ENDIF
            IF (LWAVEWIND) THEN
              FFLAG(IRU10)=.FALSE.
              GFLAG(IRU10)=.FALSE.
            ENDIF
            GFLAG(IRALTHS)=.FALSE.
            GFLAG(IRALTHSC)=.FALSE.
            GFLAG(IRALTRC)=.FALSE.
            MARSTYPE='an'
          ENDIF

!         change the output flags need to reset mapping
            CALL MPCRTBL

        ENDIF

!       NEED TO TEMPORARLY RESET THE IFS FORECAST STEP
!       IF THE IFS GRIB HEADER IS USED SUCH THAT IT POINTS TO
!       THE START OF THE RUN.
        IF (LGRHDIFS) LRSTST0=.TRUE.

!       COMPUTE OUTPUT PARAMETERS AND PRINT OUT NORMS
        IF (NIPRMOUT > 0) THEN
          ALLOCATE(BOUTST0(IJS:IJL,NIPRMOUT))
          CALL OUTBS (IJS, IJL, MIJ, FL1, XLLWS,                   &
     &                WVPRPT, WVENVI, FF_NOW, INTFLDS, NEMO2WAM,   &
     &                BOUTST0)
        ENDIF

        IF ( .NOT. LRESTARTED ) THEN
          IF (IREST == 1 .AND. MARSTYPE /= 'an' .AND. LGRIBOUT) THEN
            CALL OUTSPEC(IJS, IJL, FL1, CICOVER)
            LLFLUSH = .TRUE.
          ENDIF

          IF (NIPRMOUT > 0 ) THEN
            CALL OUTWINT(IJS, IJL, BOUTST0)
            LLFLUSH = .TRUE.
          ENDIF

        ENDIF

        IF (ALLOCATED(BOUTST0)) DEALLOCATE(BOUTST0)

        IF (LFDB .AND. LLFLUSH) THEN
           CALL GSTATS(1976,0)
           CALL IFLUSHFDBSUBS()
           CALL GSTATS(1976,1)
           LLFLUSH=.FALSE.
        ENDIF

        IF (LFRSTFLD) THEN
          MARSTYPE=MARSTYPEBAK
          FFLAG=FFLAGBAK
          GFLAG=GFLAGBAK
          LFRSTFLD=.FALSE.
!         change the output flags need to reset mapping
          CALL MPCRTBL
        ENDIF

        CDTINTT=CDTINTTBAK

        IF (LGRHDIFS) LRSTST0=.FALSE.

      ENDIF


END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('OUTSTEP0',1,ZHOOK_HANDLE)

END SUBROUTINE OUTSTEP0
