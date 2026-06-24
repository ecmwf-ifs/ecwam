SUBROUTINE UPDATEWD (KIJS, KIJL, FL1,         &
 &                   WAVNUM, CINV, XK2CG,     &
 &                   USA, WSWAVE, WDWAVE, UFRIC, &
 &                   Z0M, Z0B, CHRNCK, TAUW, TAUWDIR, &
 &                   AIRD, WSTAR, CICOVER, &
 &                   MIJ, XLLWS)         

! ----------------------------------------------------------------------

!**** *UPDATEWD* - COMPUTATION OF WINDSPEED UPDATE AFTER ASSIMILATION

!     J. BIDLOT             ECMWF  APRIL 2001

!*    PURPOSE.
!     --------

!       UPDATE THE WIND SPEED AND WAVE STRESS USING NEW UFRIC AND
!       SPRECTRUM.

!**   INTERFACE.
!     ----------

!       *CALL* *UPDATEWD (KIJS, KIJL, FL1,
!                         WVPRPT,
!                         USA, FF_NOW,
!                         MIJ, XLLWS)         

!          *KIJS*     FIRST INDEX IN BLOCK.
!          *KIJL*     LAST  INDEX IN BLOCK.
!          *FL1*      SPECTRUM.
!          *WVPRPT* - WAVE PROPERTIES
!          *FF_NOW  - FORCING FIELDS.

!     METHOD.
!     -------

!     FIND THE NEW CHARNOCK PARAMETER BASED ON THE NEW USTAR AND TAUW
!     AND THUS UPDATE THE ROUGHNESS LENGTH. WITH THIS NEW ROUGHNESS
!     LENGTH AND THE NEW USTAR FIND THE NEW U10 USING THE LOGARITHMIC
!     WIND PROFILE.

!     EXTERNALS.
!     -----------

!       NONE.

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ROWATERM1
      USE YOWFRED  , ONLY : TH
      USE YOWTABL  , ONLY : UMAX
      USE YOWWIND  , ONLY : WSPMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! -----------------------------------------------------------------------

      IMPLICIT NONE
#include "fkmean.intfb.h"
#include "sinflx.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM, CINV, XK2CG
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: USA  ! SPEUDO ANALYSED FRICTION VELOCITY
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: WSWAVE, WDWAVE, UFRIC, CICOVER, WSTAR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: Z0M, Z0B, CHRNCK, TAUW, TAUWDIR, AIRD
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(OUT) :: MIJ
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      INTEGER(KIND=JWIM) :: IJ, K
      INTEGER(KIND=JWIM) :: ICALL, NCALL

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: RAORW 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: F1MEAN, AKMEAN, XKMEAN 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EMEAN, FMEAN, FMEANWS, HALP
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: PHIWA
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: COSWDIF, SINWDIF2 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: FLM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: RHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: FLD, SL, SPOS

      LOGICAL :: LUPDTUS

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('UPDATEWD',0,ZHOOK_HANDLE)

!     UPDATE WIND SPEED:
      DO IJ=KIJS,KIJL
        WSWAVE(IJ) = (USA(IJ)/UFRIC(IJ))*WSWAVE(IJ)
        WSWAVE(IJ) = MIN(WSWAVE(IJ),UMAX)
        WSWAVE(IJ) = MAX(WSWAVE(IJ),WSPMIN)
      ENDDO

!     ADJUST UFRIC ZONEW AND TAUW

      CALL FKMEAN(KIJS, KIJL, FL1, WAVNUM,               &
     &            EMEAN, FMEAN, F1MEAN, AKMEAN, XKMEAN)

      LUPDTUS = .TRUE.
      NCALL = 2
      FMEANWS(KIJS:KIJL) = 0.0_JWRB    ! setting to zero so that it is not used in the first call to sinflx
      FLM(KIJS:KIJL,:) = 0.0_JWRB

      DO IJ=KIJS,KIJL
        RAORW(IJ) = MAX(AIRD(IJ), 1.0_JWRB) * ROWATERM1
      ENDDO

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          COSWDIF(IJ,K) = COS(TH(K)-WDWAVE(IJ))
          SINWDIF2(IJ,K) = SIN(TH(K)-WDWAVE(IJ))**2
        ENDDO
      ENDDO

      DO ICALL = 1, NCALL 
        CALL SINFLX (ICALL, NCALL, KIJS, KIJL,                      &
     &               LUPDTUS,                                       &
     &               FL1,                                           &
     &               WAVNUM, CINV, XK2CG,                           &
     &               WSWAVE, WDWAVE, AIRD, RAORW, WSTAR, CICOVER,   &
     &               COSWDIF, SINWDIF2,                             &
     &               FMEAN, HALP, FMEANWS,                          &
     &               FLM,                                           &
     &               UFRIC, TAUW, TAUWDIR, Z0M, Z0B, CHRNCK, PHIWA, &
     &               FLD, SL, SPOS,                                 &
     &               MIJ, RHOWGDFTH, XLLWS)
      ENDDO

IF (LHOOK) CALL DR_HOOK('UPDATEWD',1,ZHOOK_HANDLE)

END SUBROUTINE UPDATEWD
