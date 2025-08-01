! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE IMPLSCH (KIJS, KIJL, FL1,                         &
 &                  WAVNUM, CGROUP, CIWA, CINV, XK2CG, STOKFAC, &
 &                  EMAXDPT, DEPTH, IOBND, IODP,IBRMEM,      &
 &                  AIRD, WDWAVE, CICOVER, WSWAVE, WSTAR, USTRA, VSTRA, &
 &                  UFRIC, TAUW, TAUWDIR, Z0M, Z0B, CHRNCK, CITHICK, &
 &                  NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN, &
 &                  NPHIEPS, NTAUOC, NSWH, NMWP, NEMOTAUX, &
 &                  NEMOTAUY, NEMOTAUICX, NEMOTAUICY, &
 &                  NEMOWSWAVE, NEMOPHIF, &
 &                  WSEMEAN, WSFMEAN, USTOKES, VSTOKES, STRNMS, &
 &                  TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUOC, &
 &                  TAUICX, TAUICY, &
 &                  PHIOCD, PHIEPS, PHIAW, &
 &                  MIJ, XLLWS)

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
!    &                   WVPRPT,
!    &                   WVENVI, FF_NOW,
!    &                   INTFLDS, WAM2NEMO,
!    &                   MIJ,  XLLWS)
!      *KIJS*    - LOCAL INDEX OF FIRST GRIDPOINT
!      *KIJL*    - LOCAL INDEX OF LAST GRIDPOINT
!      *FL1*     - FREQUENCY SPECTRUM(INPUT AND OUTPUT).
!      *WVPRPT*  - WAVE PROPERTIES FIELDS
!      *WVENVI*  - WAVE ENVIRONMENT  
!      *FF_NOW*    FORCING FIELDS
!      *INTFLDS*   INTEGRATED/DERIVED PARAMETERS
!      *WAM2NEMO*  WAVE FIELDS PASSED TO NEMO
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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY, FORCING_FIELDS,   &
 &                             INTGT_PARAM_FIELDS, WAVE2OCEAN

      USE YOWCOUP  , ONLY : LWFLUX   , LWVFLX_SNL , LWNEMOCOU,           &
                            LWNEMOCOUSTRN, LWNEMOCOUIBR, LWNEMOCOUWRS
      USE YOWCOUT  , ONLY : LWFLUXOUT 
      USE YOWFRED  , ONLY : FR       ,TH       ,COFRM4    ,FLMAX
      USE YOWICE   , ONLY : FLMIN    ,LICERUN   ,LMASKICE ,              &
                            LCIWA1   ,LCIWA2    ,LCIWA3   ,LCISCAL   ,   &
 &                          ZALPFACX
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,LLUNSTR
      USE YOWPCONS , ONLY : WSEMEAN_MIN, ROWATERM1
      USE YOWSTAT  , ONLY : IDELT    ,LBIWBK   ,XIMP
      USE YOWWNDG  , ONLY : ICODE    ,ICODE_CPL
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "femeanws.intfb.h"
#include "fkmean.intfb.h"
#include "sbottom.intfb.h"
#include "imphftail.intfb.h"
#include "sdepthlim.intfb.h"
#include "sdissip.intfb.h"
#include "sdiwbk.intfb.h"
#include "sdice.intfb.h"
#include "icebreak_modify_attenuation.intfb.h"
#include "setice.intfb.h"
#include "sinflx.intfb.h"
#include "snonlin.intfb.h"
#include "stokestrn.intfb.h"
#include "wnfluxes.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: WAVNUM
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: CGROUP
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: CIWA
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: CINV
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: XK2CG
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: STOKFAC

      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: EMAXDPT
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: DEPTH
      INTEGER(KIND=JWIM), DIMENSION(KIJL), INTENT(IN) :: IODP
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: IBRMEM      
      INTEGER(KIND=JWIM), DIMENSION(KIJL), INTENT(IN) :: IOBND

      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: AIRD, WDWAVE, CICOVER, WSWAVE, WSTAR, USTRA, VSTRA
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: UFRIC, TAUW, TAUWDIR, Z0M, Z0B, CHRNCK, CITHICK
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: WSEMEAN, WSFMEAN, USTOKES, VSTOKES, STRNMS
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUOC, PHIOCD
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: TAUICX, TAUICY
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: PHIEPS, PHIAW
      REAL(KIND=JWRO), DIMENSION(KIJL), INTENT(INOUT) :: NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN
      REAL(KIND=JWRO), DIMENSION(KIJL), INTENT(INOUT) :: NPHIEPS, NTAUOC, NSWH, NMWP, NEMOTAUX
      REAL(KIND=JWRO), DIMENSION(KIJL), INTENT(INOUT) :: NEMOTAUY, NEMOWSWAVE, NEMOPHIF
      REAL(KIND=JWRO), DIMENSION(KIJL), INTENT(INOUT) :: NEMOTAUICX, NEMOTAUICY
      INTEGER(KIND=JWIM), DIMENSION(KIJL), INTENT(OUT) :: MIJ
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      INTEGER(KIND=JWIM) :: IJ, K, M
      INTEGER(KIND=JWIM) :: ICALL, NCALL
      
      REAL(KIND=JWRB) :: DELTM, DELT5, DELT
      REAL(KIND=JWRB) :: GTEMP1, GTEMP2, FLHAB, BETA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: DELFL(NFRE)
      REAL(KIND=JWRB), DIMENSION(KIJL) :: RAORW
      REAL(KIND=JWRB), DIMENSION(KIJL) :: EMEAN, FMEAN, HALP
      REAL(KIND=JWRB), DIMENSION(KIJL) :: EMEANWS, FMEANWS, USFM
      REAL(KIND=JWRB), DIMENSION(KIJL) :: F1MEAN, AKMEAN, XKMEAN 
      REAL(KIND=JWRB), DIMENSION(KIJL) :: PHIWA

      REAL(KIND=JWRB), DIMENSION(KIJL,NANG) :: FLM
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG) :: COSWDIF, SINWDIF2
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: TEMP
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: RHOWGDFTH
!     *FLD* DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!     *SL*  TOTAL SOURCE FUNCTION ARRAY.
!     *SPOS* : POSITIVE SINPUT ONLY
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE) :: FLD, SL, SPOS
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE) :: SSOURCE 
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE) :: SLICE

      REAL(KIND=JWRB),    DIMENSION(KIJL) :: ALPFAC

      LOGICAL :: LCFLX
      LOGICAL :: LUPDTUS

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('IMPLSCH',0,ZHOOK_HANDLE)


!*    1. INITIALISATION.
!        ---------------

      DELT  = IDELT
      DELTM = 1.0_JWRB/DELT
      DELT5 = XIMP*DELT

      LCFLX=LWFLUX.OR.LWFLUXOUT.OR.LWNEMOCOU


      DO IJ=KIJS,KIJL
        RAORW(IJ) = MAX(AIRD(IJ), 1.0_JWRB) * ROWATERM1
      ENDDO

      DO IJ=KIJS,KIJL 
        ALPFAC(IJ)  = ZALPFACX ! <1=some reduction, 1=no reduction to attenuation
      ENDDO

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          COSWDIF(IJ,K) = COS(TH(K)-WDWAVE(IJ))
          SINWDIF2(IJ,K) = SIN(TH(K)-WDWAVE(IJ))**2
        ENDDO
      ENDDO

      IF (LWNEMOCOUWRS .AND. (.NOT. (LCIWA1 .OR. LCIWA2 .OR. LCIWA3)))  THEN
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              SLICE(IJ,K,M) = 0.0_JWRB
            ENDDO
          ENDDO
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    2. COMPUTATION OF IMPLICIT INTEGRATION.
!        ------------------------------------

!         INTEGRATION IS DONE FROM CDATE UNTIL CDTPRO FOR A BLOCK
!         OF LATITUDES BETWEEN PROPAGATION CALLS.


!     REDUCE WAVE ENERGY IF LARGER THAN DEPTH LIMITED WAVE HEIGHT
      IF (LBIWBK) THEN
         !$loki inline
         CALL SDEPTHLIM(KIJS, KIJL, EMAXDPT, FL1)
      ENDIF

!*    2.2 COMPUTE MEAN PARAMETERS.
!        ------------------------

      !$loki inline
      CALL FKMEAN(KIJS, KIJL, FL1, WAVNUM,                    &
     &            EMEAN, FMEAN, F1MEAN, AKMEAN, XKMEAN)

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          ! still allow noise in full sea ice cover, but only ten percent
          FLM(IJ,K) = (1._JWRB - 0.9_JWRB*MIN(CICOVER(IJ),0.99_JWRB))*FLMIN*MAX(0.0_JWRB, COSWDIF(IJ,K))**2 
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    2.3 COMPUTATION OF SOURCE FUNCTIONS.
!         --------------------------------

!*    2.3.1 ITERATIVELY UPDATE STRESS AND COMPUTE WIND INPUT TERMS. 
!           -------------------------------------------------------

      LUPDTUS = .TRUE.
      NCALL = 2
      DO ICALL = 1, NCALL 
        !$loki inline
        CALL SINFLX (ICALL, NCALL, KIJS, KIJL,  &
     &               LUPDTUS,                   &
     &               FL1,                       &
     &               WAVNUM, CINV, XK2CG,       &
     &               WSWAVE, WDWAVE, AIRD,      &
     &               RAORW, WSTAR, CICOVER,     &
     &               COSWDIF, SINWDIF2,         &
     &               FMEAN, HALP, FMEANWS,      &
     &               FLM,                       &
     &               UFRIC, TAUW, TAUWDIR,      &
     &               Z0M, Z0B, CHRNCK, PHIWA,   &
     &               FLD, SL, SPOS,             &
     &               MIJ, RHOWGDFTH, XLLWS)

      ENDDO

!     2.3.3 ADD THE OTHER SOURCE TERMS.
!           ---------------------------

      !$loki inline
      CALL SDISSIP (KIJS, KIJL, FL1 ,FLD, SL,   &
     &              WAVNUM, CGROUP, XK2CG,      &
     &              EMEAN, F1MEAN, XKMEAN,      &
     &              UFRIC, COSWDIF, RAORW)

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


      !$loki inline
      CALL SDIWBK(KIJS, KIJL, FL1 ,FLD, SL, DEPTH, EMAXDPT, EMEAN, F1MEAN)

      IF ( LICERUN ) THEN

!        Use linear scaling of ALL proceeding source terms under sea ice (this is a complete unknown)
         IF (LCISCAL) THEN
           DO M = 1,NFRE
             DO K = 1,NANG
               DO IJ = KIJS,KIJL
                  BETA=1._JWRB-CICOVER(IJ)
                  SL(IJ,K,M)  = BETA*SL(IJ,K,M)  
                  FLD(IJ,K,M) = BETA*FLD(IJ,K,M) 
               END DO
             END DO
           END DO
         ENDIF

!        Coupling of waves and sea ice (type 1): wave-induced sea ice break up + reduced attenuation
         IF(LWNEMOCOUIBR) THEN 
          !$loki inline
          CALL ICEBREAK_MODIFY_ATTENUATION (KIJS,KIJL,IBRMEM,ALPFAC)           
         ENDIF

!        Attenuation of waves in ice
         IF(LCIWA1 .OR. LCIWA2 .OR. LCIWA3) THEN
            !$loki inline
            CALL SDICE (KIJS, KIJL, FL1, FLD, SL, SLICE, WAVNUM, CGROUP, CICOVER, CITHICK, ALPFAC)
         ENDIF

      ENDIF

      !$loki inline
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
        USFM(IJ) = UFRIC(IJ)*MAX(FMEANWS(IJ), FMEAN(IJ))
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
              FL1(IJ,K,M) =  MAX(IODP(IJ)*FL1(IJ,K,M),FLM(IJ,K)) 
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
              FL1(IJ,K,M)  = MAX(FL1(IJ,K,M),FLM(IJ,K))
              SSOURCE(IJ,K,M) = SSOURCE(IJ,K,M) + DELTM * MIN(FLMAX(M)-FL1(IJ,K,M),0.0_JWRB)
              FL1(IJ,K,M) = MIN(FL1(IJ,K,M),FLMAX(M))
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF (LCFLX) THEN
        !$loki inline
        CALL WNFLUXES (KIJS, KIJL,                       &
     &                 MIJ, RHOWGDFTH,                   &
     &                 CINV,                             &
     &                 SSOURCE, SLICE, CICOVER,          &
     &                 PHIWA,                            &
     &                 EMEAN, F1MEAN, WSWAVE, WDWAVE,    &
     &                 USTRA, VSTRA,                     &
     &                 UFRIC, AIRD,                      &
     &                 NPHIEPS, NTAUOC, NSWH, NMWP,      &
     &                 NEMOTAUX, NEMOTAUY,               &
     &                 NEMOTAUICX, NEMOTAUICY,           &
     &                 NEMOWSWAVE, NEMOPHIF,             &
     &                 TAUXD, TAUYD, TAUOCXD, TAUOCYD,   &
     &                 TAUOC, TAUICX, TAUICY,            &
     &                 PHIOCD, PHIEPS, PHIAW,            &
     &                 .TRUE.)
      ENDIF
! ----------------------------------------------------------------------

!*    2.5 REPLACE DIAGNOSTIC PART OF SPECTRA BY A F**(-5) TAIL.
!         -----------------------------------------------------

      !$loki inline
      CALL FKMEAN(KIJS, KIJL, FL1, WAVNUM,                      &
     &            EMEAN, FMEAN, F1MEAN, AKMEAN, XKMEAN)

!     MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
     !$loki inline
      CALL FEMEANWS(KIJS, KIJL, FL1, XLLWS, FMEANWS, EMEANWS)

      !$loki inline
      CALL IMPHFTAIL(KIJS, KIJL, MIJ, FLM, WAVNUM, XK2CG, FL1)


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
        !$loki inline
        CALL SETICE(KIJS, KIJL, FL1, CICOVER, COSWDIF)
      ENDIF


!*    2.7 SURFACE STOKES DRIFT AND STRAIN IN SEA ICE
!         ------------------------------------------

      !$loki inline
      CALL STOKESTRN(KIJS, KIJL, FL1, WAVNUM, STOKFAC, DEPTH, WSWAVE, WDWAVE, CICOVER, CITHICK, &
 &                   USTOKES, VSTOKES, STRNMS, NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN)

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('IMPLSCH',1,ZHOOK_HANDLE)


END SUBROUTINE IMPLSCH
