      SUBROUTINE WDFLUXES (KIJS, KIJL,                        &
     &                     MIJ,                               &
     &                     FL1, XLLWS,                        &
     &                     WVPRPT,                            &
     &                     DEPTH,                             &
     &                     FF_NOW,                            &
     &                     INTFLDS)

! ----------------------------------------------------------------------

!**** *WDFLUXES* - IF NEEDED, IT EVALUATES SINPUT AND SDISSIP
!****              FOR THE CALCULATION OF OUTPUT WAVE DEPENDANT FLUXES.
!****              THIS IS NORMALLY DONE IN *IMPLSCH* BUT AT INITILISATION,
!****              THE SOURCES TERMS ARE NOT YET EVALUATED.
!****              ALSO DETERMINE THE SURFACE STOKES DRIFT ANd STRIN IN THE ICE.

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *WDFLUXES (KIJS, KIJL,
!    &                    MIJ,
!    &                    FL1, XLLWS,
!    &                    WVPRPT,
!    &                    FF_NOW,
!    &                    INTFLDS)

!          *KIJS*   - INDEX OF FIRST GRIDPOINT.
!          *KIJL*   - INDEX OF LAST GRIDPOINT.
!          *FL1*    - SPECTRUM(INPUT).
!          *XLLWS*  - TOTAL WINDSEA MASK FROM INPUT SOURCE TERM
!          *WVPRPT* - WAVE PROPERTIES FIELDS
!          *DEPTH*  - WATER DEPTH.
!          *FF_NOW* - FORCING FIELDS AT CURRENT TIME.
!          *INTFLDS*-  INTEGRATED/DERIVED PARAMETERS

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : FREQUENCY, FORCING_FIELDS, INTGT_PARAM_FIELDS

      USE YOWCOUP  , ONLY : LWFLUX   ,LWVFLX_SNL, LWNEMOCOUSTRN
      USE YOWCOUT  , ONLY : LWFLUXOUT 
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : WSEMEAN_MIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "femeanws.intfb.h"
#include "fkmean.intfb.h"
#include "sdissip.intfb.h"
#include "sinflx.intfb.h"
#include "snonlin.intfb.h"
#include "stokestrn.intfb.h"
#include "wnfluxes.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM),  DIMENSION(KIJS:KIJL), INTENT(OUT) :: MIJ

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: FL1
      TYPE(FREQUENCY), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WVPRPT
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH
      TYPE(FORCING_FIELDS), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: FF_NOW
      TYPE(INTGT_PARAM_FIELDS), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: INTFLDS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      INTEGER(KIND=JWIM) :: IJ, K, M
      INTEGER(KIND=JWIM) :: ICALL, NCALL

      REAL(KIND=JWRB) :: TAU, XN, PHIDIAG, TAUO
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TAUW_LOC  ! TAUW should not be updated do use a local array
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TAUWDIR_LOC  ! TAUW should not be updated do use a local array
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EMEAN, FMEAN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EMEANWS, FMEANWS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: F1MEAN, AKMEAN, XKMEAN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: PHIWA
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: FLM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: RHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: FLD, SL, SPOS

      LOGICAL :: LCFLX
      LOGICAL :: LUPDTUS

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WDFLUXES',0,ZHOOK_HANDLE)

ASSOCIATE(WAVNUM => WVPRPT%WAVNUM, &
 &        CINV => WVPRPT%CINV, &
 &        CGROUP => WVPRPT%CGROUP, &
 &        STOKFAC => WVPRPT%STOKFAC, &
 &        WSWAVE => FF_NOW%WSWAVE, &
 &        WDWAVE => FF_NOW%WDWAVE, &
 &        UFRIC => FF_NOW%UFRIC, &
 &        Z0M => FF_NOW%Z0M, &
 &        Z0B => FF_NOW%Z0B, &
 &        AIRD => FF_NOW%AIRD, &
 &        WSTAR => FF_NOW%WSTAR, &
 &        CICOVER => FF_NOW%CICOVER, &
 &        WSEMEAN => INTFLDS%WSEMEAN, &
 &        WSFMEAN => INTFLDS%WSFMEAN)

!*    1. INITIALISATION.
!        ---------------

      LCFLX=LWFLUX.OR.LWFLUXOUT
! ----------------------------------------------------------------------

!*    1.2 COMPUTATION OF RELEVANT SOURCE FUNCTIONS.
!         -----------------------------------------

      CALL FKMEAN(KIJS, KIJL, FL1, WAVNUM,                         &
     &            EMEAN, FMEAN, F1MEAN, AKMEAN, XKMEAN)

      TAUW_LOC(:) = 0.0_JWRB
      TAUWDIR_LOC(:) = WDWAVE(:)

      LUPDTUS = .FALSE.
      FMEANWS(:) = FMEAN(:)
      FLM(:,:) = 0.0_JWRB
      NCALL = 1
      ICALL = 1
      CALL SINFLX (ICALL, NCALL, KIJS, KIJL,                        &
     &             LUPDTUS,                                         &
     &             FL1,                                             &
     &             WAVNUM, CINV, CGROUP,                            &
     &             WSWAVE, WDWAVE, AIRD, WSTAR, CICOVER,            &
     &             FMEAN, FMEANWS,                                  &
     &             FLM,                                             &
     &             UFRIC, TAUW_LOC, TAUWDIR_LOC, Z0M, Z0B, PHIWA, &
     &             FLD, SL, SPOS,                                   &
     &             MIJ, RHOWGDFTH, XLLWS)

      IF (LCFLX) THEN

        CALL SDISSIP (KIJS, KIJL, FL1 ,FLD, SL,    &
     &                WAVNUM, CGROUP,              &
     &                EMEAN, F1MEAN, XKMEAN,       &
     &                UFRIC, WDWAVE, AIRD)

        IF (.NOT. LWVFLX_SNL) THEN
          CALL WNFLUXES (KIJS, KIJL,                        &
     &                   MIJ, RHOWGDFTH,                    &
     &                   CINV,                              &
     &                   SL, CICOVER,                       &
     &                   PHIWA,                             &
     &                   EMEAN, F1MEAN, WSWAVE, WDWAVE,     &
     &                   UFRIC, AIRD, INTFLDS, .FALSE.)
        ENDIF

        CALL SNONLIN (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN)

        IF (LWVFLX_SNL) THEN
          CALL WNFLUXES (KIJS, KIJL,                        &
     &                   MIJ, RHOWGDFTH,                    &
     &                   CINV,                              &
     &                   SL, CICOVER,                       &
     &                   PHIWA,                             &
     &                   EMEAN, F1MEAN, WSWAVE, WDWAVE,     &
     &                   UFRIC, AIRD, INTFLDS, .FALSE.)
        ENDIF

        IF (LWFLUX) THEN
         CALL FEMEANWS(KIJS, KIJL, FL1, XLLWS, EMEANWS, FMEANWS)

          DO IJ=KIJS,KIJL
            IF (EMEANWS(IJ) < WSEMEAN_MIN) THEN
              WSEMEAN(IJ) = WSEMEAN_MIN 
              WSFMEAN(IJ) = 2._JWRB*FR(NFRE)
            ELSE
              WSEMEAN(IJ) = EMEANWS(IJ)
              WSFMEAN(IJ) = FMEANWS(IJ) 
            ENDIF
          ENDDO
        ENDIF

        CALL STOKESTRN(KIJS, KIJL, FL1, WAVNUM, STOKFAC, DEPTH, FF_NOW, INTFLDS)

      ENDIF
! ----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('WDFLUXES',1,ZHOOK_HANDLE)

END SUBROUTINE WDFLUXES

