      SUBROUTINE IMPLSCH (KIJS, KIJL, FL1,                         &
     &                    WAVNUM, CINV,                            &
     &                    CGROUP, STOKFAC,                         &
     &                    DEPTH, EMAXDPT,                          &
     &                    TAUW, TAUWDIR,                           &
     &                    CICOVER, CIWA,                           &
     &                    WSWAVE, WDWAVE, UFRIC,                   &
     &                    Z0M, Z0B,                                &
     &                    AIRD, WSTAR,                             &
     &                    WSEMEAN, WSFMEAN,                        &
     &                    USTOKES, VSTOKES, STRNMS,                &
     &                    MIJ, XLLWS)

! ----------------------------------------------------------------------

!**** *IMPLSCH* - IMPLICIT SCHEME FOR TIME INTEGRATION OF SOURCE
!****             FUNCTIONS.


!*    PURPOSE.
!     --------

!       THE IMPLICIT SCHEME ENABLES THE USE OF A TIMESTEP WHICH IS
!       LARGE COMPARED WITH THE CHARACTERISTIC DYNAMIC TIME SCALE.
!       THE SCHEME IS REQUIRED FOR THE HIGH FREQUENCIES WHICH
!       RAPIDLY ADJUST TO A QUASI-EQUILIBRIUM.

!**   INTERFACE.
!     ----------

!       *CALL* *IMPLSCH (KIJS, KIJL, FL1,
!    &                   WAVNUM, CINV, CGROUP, STOKFAC,
!    &                   DEPTH, EMAXDPT,
!    &                   TAUW,TAUWDIR,
!    &                   CICOVER, CIWA,
!    &                   WSWAVE,WDWAVE,UFRIC,Z0M,AIRD,WSTAR,
!    &                   USTOKES, VSTOKES, STRNMS,
!    &                   MIJ,  XLLWS)
!      *KIJS*    - LOCAL INDEX OF FIRST GRIDPOINT
!      *KIJL*    - LOCAL INDEX OF LAST GRIDPOINT
!      *FL1*     - FREQUENCY SPECTRUM(INPUT AND OUTPUT).
!      *WAVNUM*  - WAVE NUMBER.
!      *CINV*    - INVERSE PHASE VELOCITY.
!      *CGROUP*  - GROUP SPPED.
!      *STOKFAC* - FACTOR TO COMPUTE THE STOKES SURFACE DRIFT
!      *DEPTH*     WATER DEPTH
!      *EMAXDPT*   MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
!      *WSWAVE*    WIND SPEED IN M/S.
!      *WDWAVE*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!      *UFRIC*     FRICTION VELOCITY IN M/S.
!      *Z0M*       ROUGHNESS LENGTH IN M.
!      *Z0B*       BACKGROUND ROUGHNESS LENGTH.
!      *TAUW*      WAVE STRESS IN (M/S)**2
!      *TAUWDIR*   WAVE STRESS DIRECTION. 
!      *AIRD*      AIR DENSITY IN KG/M3.
!      *WSTAR*  FREE CONVECTION VELOCITY SCALE (M/S)
!      *CICOVER*   SEA ICE COVER.
!      *CIWA*      SEA ICE WAVE ATTENUATION.
!      *WSEMEAN*   WINDSEA VARIANCE.
!      *WSFMEAN*   WINDSEA MEAN FREQUENCY.
!      *USTOKES*   U-COMP SURFACE STOKES DRIFT.
!      *VSTOKES*   V-COMP SURFACE STOKES DRIFT.
!      *STRNMS*    MEAN SQUARE STRAIN INTO THE SEA ICE (only if LWNEMOCOUSTRN).
!      *MIJ*       LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!      *XLLWS*     TOTAL WINDSEA MASK FROM INPUT SOURCE TERM


!     METHOD.
!     -------

!       THE SPECTRUM AT TIME (TN+1) IS COMPUTED AS
!       FN+1=FN+DELT*(SN+SN+1)/2., WHERE SN IS THE TOTAL SOURCE
!       FUNCTION AT TIME TN, SN+1=SN+(DS/DF)*DF - ONLY THE DIAGONAL
!       TERMS OF THE FUNCTIONAL MATRIX DS/DF ARE YOWPUTED, THE
!       NONDIAGONAL TERMS ARE NEGLIGIBLE.
!       THE ROUTINE IS CALLED AFTER PROPAGATION FOR TIME PERIOD
!       BETWEEN TWO PROPAGATION CALLS - ARRAY FL1 CONTAINS THE
!       SPECTRUM AND FL IS USED AS AN INTERMEDIATE STORAGE FOR THE
!       DIAGONAL TERM OF THE FUNCTIONAL MATRIX.


!     REFERENCE.
!     ----------

!       S. HASSELMANN AND K. HASSELMANN, "A GLOBAL WAVE MODEL",
!       30/6/85 (UNPUBLISHED NOTE)

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWFLUX   , LWVFLX_SNL , LWNEMOCOU, LWNEMOCOUSTRN 
      USE YOWCOUT  , ONLY : LWFLUXOUT 
      USE YOWFRED  , ONLY : FR       ,TH       ,COFRM4    ,FLMAX
      USE YOWICE   , ONLY : FLMIN    ,LCIWABR  ,LICERUN   ,LMASKICE
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : WSEMEAN_MIN 
      USE YOWSHAL  , ONLY : IODP     ,IOBND
      USE YOWSTAT  , ONLY : IDELT    ,CDTPRO   ,LBIWBK
      USE YOWUNPOOL, ONLY : LLUNSTR
      USE YOWWNDG  , ONLY : ICODE    ,ICODE_CPL

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "imphftail.intfb.h"
#include "sdepthlim.intfb.h"
#include "cimsstrn.intfb.h"
#include "ciwabr.intfb.h"
#include "femeanws.intfb.h"
#include "fkmean.intfb.h"
#include "sbottom.intfb.h"
#include "sdissip.intfb.h"
#include "sdiwbk.intfb.h"
#include "setice.intfb.h"
#include "sinflx.intfb.h"
#include "snonlin.intfb.h"
#include "stokesdrift.intfb.h"
#include "wnfluxes.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: FL1
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(OUT) :: MIJ
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM, CINV, CGROUP, STOKFAC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH, EMAXDPT
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: TAUW, TAUWDIR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: CICOVER
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: WSWAVE, UFRIC 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: WDWAVE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: Z0M
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: Z0B
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: AIRD, WSTAR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: WSEMEAN, WSFMEAN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: USTOKES, VSTOKES, STRNMS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: CIWA
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      INTEGER(KIND=JWIM) :: IJ, K, M
      INTEGER(KIND=JWIM) :: ICALL, NCALL

      REAL(KIND=JWRB) :: DELT, DELTM, XIMP, DELT5
      REAL(KIND=JWRB) :: GTEMP1, GTEMP2, FLHAB
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: DELFL(NFRE)
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EMEAN, FMEAN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EMEANWS, FMEANWS, USFM, GADIAG 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: F1MEAN, AKMEAN, XKMEAN 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: PHIWA
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: DPTHREDUC

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: FLM 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: TEMP
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: RHOWGDFTH
!     *FLD* DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!     *SL*  TOTAL SOURCE FUNCTION ARRAY.
!     *SPOS* : POSITIVE SINPUT ONLY
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: FLD, SL, SPOS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: CIREDUC 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: SSOURCE 

      LOGICAL :: LCFLX
      LOGICAL :: LUPDTUS

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('IMPLSCH',0,ZHOOK_HANDLE)


!*    1. INITIALISATION.
!        ---------------

      DELT = IDELT
      DELTM = 1.0_JWRB/DELT
      XIMP = 1.0_JWRB
      DELT5 = XIMP*DELT

      LCFLX=LWFLUX.OR.LWFLUXOUT.OR.LWNEMOCOU

! ----------------------------------------------------------------------

!*    2. COMPUTATION OF IMPLICIT INTEGRATION.
!        ------------------------------------

!         INTEGRATION IS DONE FROM CDATE UNTIL CDTPRO FOR A BLOCK
!         OF LATITUDES BETWEEN PROPAGATION CALLS.


!     REDUCE WAVE ENERGY IF LARGER THAN DEPTH LIMITED WAVE HEIGHT
      IF (LBIWBK) THEN
         CALL SDEPTHLIM(KIJS, KIJL, EMAXDPT, FL1)
      ENDIF

!*    2.2 COMPUTE MEAN PARAMETERS.
!        ------------------------

      CALL FKMEAN(KIJS, KIJL, FL1, WAVNUM,                    &
     &            EMEAN, FMEAN, F1MEAN, AKMEAN, XKMEAN)

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          FLM(IJ,K)=FLMIN*MAX(0.0_JWRB,COS(TH(K)-WDWAVE(IJ)))**2
        ENDDO
      ENDDO

!     COMPUTE DAMPING COEFFICIENT DUE TO FRICTION ON BOTTOM OF THE SEA ICE.
!!! testing sea ice attenuation (might need to restrict usage when needed)
      IF (LCIWABR) THEN
        CALL CIWABR(KIJS, KIJL, CICOVER, FL1, WAVNUM, CGROUP, CIREDUC)
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              CIREDUC(IJ,K,M)=CIWA(IJ,M)*CIREDUC(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              CIREDUC(IJ,K,M)=CIWA(IJ,M)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    2.3 COMPUTATION OF SOURCE FUNCTIONS.
!         --------------------------------

!*    2.3.1 ITERATIVELY UPDATE STRESS AND COMPUTE WIND INPUT TERMS. 
!           -------------------------------------------------------

      LUPDTUS = .TRUE.
      NCALL = 2
      DO ICALL = 1, NCALL 
        CALL SINFLX (ICALL, NCALL, KIJS, KIJL,                 &
     &               LUPDTUS,                                  &
     &               FL1,                                      &
     &               WAVNUM, CINV, CGROUP,                     &
     &               WSWAVE, WDWAVE, AIRD, WSTAR, CICOVER,     &
     &               FMEAN, FMEANWS,                        &
     &               FLM,                                      &
     &               UFRIC, TAUW, TAUWDIR, Z0M, Z0B, PHIWA,    &
     &               FLD, SL, SPOS,                            &
     &               MIJ, RHOWGDFTH, XLLWS)

      ENDDO

!     2.3.3 ADD THE OTHER SOURCE TERMS.
!           ---------------------------

      CALL SDISSIP (KIJS, KIJL, FL1 ,FLD, SL,     &
     &              WAVNUM, CGROUP,               &
     &              EMEAN, F1MEAN, XKMEAN,     &
     &              UFRIC, WDWAVE, AIRD)

!     Save source term contributions relevant for the calculation of ocean fluxes
      IF (LCFLX .AND. .NOT.LWVFLX_SNL) THEN
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              SSOURCE(IJ,K,M) = SL(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      CALL SNONLIN (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN)

      IF (LCFLX .AND. LWVFLX_SNL) THEN
!     Save source term contributions relevant for the calculation of ocean fluxes
!!!!!!  SL must only contain contributions contributed to fluxes into the oceans
!       MODULATE SL BY IMPLICIT FACTOR
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              GTEMP1 = MAX((1.0_JWRB-DELT5*FLD(IJ,K,M)),1.0_JWRB)
              SSOURCE(IJ,K,M) = SL(IJ,K,M)/GTEMP1
            ENDDO
          ENDDO
        ENDDO
      ENDIF


      CALL SDIWBK(KIJS, KIJL, FL1 ,FLD, SL, DEPTH, EMAXDPT, EMEAN, F1MEAN)

      CALL SBOTTOM (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH)

! ----------------------------------------------------------------------

!*    2.4 COMPUTATION OF NEW SPECTRA.
!         ---------------------------

!     INCREASE OF SPECTRUM IN A TIME STEP IS LIMITED TO A FINITE
!     FRACTION OF A TYPICAL F**(-4) EQUILIBRIUM SPECTRUM.

      DO M=1,NFRE
        DELFL(M) = COFRM4(M)*DELT
      ENDDO
      DO IJ=KIJS,KIJL
        USFM(IJ) = UFRIC(IJ)*MAX(FMEANWS(IJ),FMEAN(IJ))
      ENDDO

      DO M=1,NFRE
        DO IJ=KIJS,KIJL
          TEMP(IJ,M) = USFM(IJ)*DELFL(M)
        ENDDO
      ENDDO

      IF (LLUNSTR) THEN
      DO K=1,NANG
        DO M=1,NFRE
          DO IJ=KIJS,KIJL
            GTEMP1 = MAX((1.0_JWRB-DELT5*FLD(IJ,K,M)),1.0_JWRB)
            GTEMP2 = DELT*SL(IJ,K,M)/GTEMP1
            FLHAB = ABS(GTEMP2)
            FLHAB = MIN(FLHAB,TEMP(IJ,M))
            FL1(IJ,K,M) = FL1(IJ,K,M) + IOBND(IJ)*SIGN(FLHAB,GTEMP2)
            FL1(IJ,K,M) = MAX(IODP(IJ)*CIREDUC(IJ,K,M)*FL1(IJ,K,M),FLM(IJ,K))
            SSOURCE(IJ,K,M) = SSOURCE(IJ,K,M) + DELTM * MIN(FLMAX(M)-FL1(IJ,K,M),0.0_JWRB)
            FL1(IJ,K,M) = MIN(FL1(IJ,K,M),FLMAX(M))
          ENDDO
        ENDDO
      ENDDO
      ELSE
      DO K=1,NANG
        DO M=1,NFRE
          DO IJ=KIJS,KIJL
            GTEMP1 = MAX((1.0_JWRB-DELT5*FLD(IJ,K,M)),1.0_JWRB)
            GTEMP2 = DELT*SL(IJ,K,M)/GTEMP1
            FLHAB = ABS(GTEMP2)
            FLHAB = MIN(FLHAB,TEMP(IJ,M))
            FL1(IJ,K,M) = FL1(IJ,K,M) + SIGN(FLHAB,GTEMP2)
            FL1(IJ,K,M) = MAX(CIREDUC(IJ,K,M)*FL1(IJ,K,M),FLM(IJ,K))
            SSOURCE(IJ,K,M) = SSOURCE(IJ,K,M) + DELTM * MIN(FLMAX(M)-FL1(IJ,K,M),0.0_JWRB)
            FL1(IJ,K,M) = MIN(FL1(IJ,K,M),FLMAX(M))
          ENDDO
        ENDDO
      ENDDO
      ENDIF

      IF (LCFLX) THEN
        CALL WNFLUXES (KIJS, KIJL,                              &
     &                 MIJ, RHOWGDFTH,                          &
     &                 CINV,                                    &
     &                 SSOURCE, CICOVER,                        &
     &                 PHIWA,                                   &
     &                 EMEAN, F1MEAN, WSWAVE, WDWAVE,        &
     &                 UFRIC, AIRD, .TRUE.)
      ENDIF
! ----------------------------------------------------------------------

!*    2.5 REPLACE DIAGNOSTIC PART OF SPECTRA BY A F**(-5) TAIL.
!         -----------------------------------------------------

      CALL FKMEAN(KIJS, KIJL, FL1, WAVNUM,                      &
     &            EMEAN, FMEAN, F1MEAN, AKMEAN, XKMEAN)

!     MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
      CALL FEMEANWS(KIJS, KIJL, FL1, XLLWS, EMEANWS, FMEANWS)

      CALL IMPHFTAIL(KIJS, KIJL, MIJ, FLM, WAVNUM, CGROUP, FL1)


!     UPDATE WINDSEA VARIANCE AND MEAN FREQUENCY IF PASSED TO ATMOSPHERE
!     ------------------------------------------------------------------
      IF (LWFLUX) THEN
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


!*    2.6 SET FL1 ON ICE POINTS TO ZERO
!         -----------------------------

      IF (LICERUN .AND. LMASKICE) THEN
        CALL SETICE(KIJS, KIJL, FL1, CICOVER, WSWAVE, WDWAVE)
      ENDIF


!*    2.7 SURFACE STOKES DRIFT AND STRAIN IN SEA ICE
!         ------------------------------------------

      CALL STOKESDRIFT(KIJS, KIJL, FL1, STOKFAC, WSWAVE, WDWAVE, CICOVER, USTOKES, VSTOKES)

      IF (LWNEMOCOUSTRN) CALL CIMSSTRN(KIJS, KIJL, FL1, WAVNUM, DEPTH, STRNMS)


! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('IMPLSCH',1,ZHOOK_HANDLE)

      END SUBROUTINE IMPLSCH
