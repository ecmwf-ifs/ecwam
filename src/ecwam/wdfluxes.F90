! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WDFLUXES (KIJS, KIJL,                         &
     &                     MIJ,                                &
     &                     FL1, XLLWS,                         &
     &                     WAVNUM, CINV, CGROUP,               &
     &                     XK2CG, STOKFAC,                     &
     &                     DEPTH, IBRMEM,                      &
     &                     WSWAVE, WDWAVE,                     &
     &                     AIRD, WSTAR,                        &
     &                     USTRA, VSTRA,                       &
     &                     CICOVER,                            &
     &                     UFRIC, Z0M,                         &
     &                     Z0B, CHRNCK, CITHICK,               &
     &                     WSEMEAN, WSFMEAN,                   &
     &                     USTOKES, VSTOKES, STRNMS,           &
     &                     TAUXD, TAUYD, TAUOCXD,              &
     &                     TAUOCYD, TAUOC, TAUICX, TAUICY,     & 
     &                     PHIOCD, PHIEPS, PHIAW,              &
     &                     NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN, &
     &                     NPHIEPS, NTAUOC, NSWH,              &
     &                     NMWP,NEMOTAUX, NEMOTAUY,            &
     &                     NEMOTAUICX, NEMOTAUICY,             &
     &                     NEMOWSWAVE, NEMOPHIF)

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
!    &                    WVENVI, FF_NOW, INTFLDS, WAM2NEMO)

!          *KIJS*   - INDEX OF FIRST GRIDPOINT.
!          *KIJL*   - INDEX OF LAST GRIDPOINT.
!          *FL1*    - SPECTRUM(INPUT).
!          *XLLWS*  - TOTAL WINDSEA MASK FROM INPUT SOURCE TERM
!          *WVPRPT* - WAVE PROPERTIES FIELDS
!          *WVENVI* - WAVE ENVIRONMENT
!          *FF_NOW* - FORCING FIELDS AT CURRENT TIME.
!          *INTFLDS*-  INTEGRATED/DERIVED PARAMETERS
!          *WAM2NEMO* FIELDS PASSED FROM WAM TO NEMO

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY, FORCING_FIELDS,  &
     &                         INTGT_PARAM_FIELDS, WAVE2OCEAN

      USE YOWCOUP  , ONLY : LWFLUX   ,LWVFLX_SNL, LWNEMOCOUSTRN,        &
                            LWNEMOCOUWRS, LWNEMOCOUIBR
      USE YOWCOUT  , ONLY : LWFLUXOUT 
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWICE   , ONLY : LICERUN  ,              &
                            LCIWA1   ,LCIWA2    ,LCIWA3   ,LCISCAL   ,   &
 &                          ZALPFACX
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : WSEMEAN_MIN, ROWATERM1
      USE YOWSTAT  , ONLY : IDELT    

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "femeanws.intfb.h"
#include "fkmean.intfb.h"
#include "sdissip.intfb.h"
#include "sdice.intfb.h"
#include "icebreak_modify_attenuation.intfb.h"
#include "sinflx.intfb.h"
#include "snonlin.intfb.h"
#include "stokestrn.intfb.h"
#include "wnfluxes.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM),  DIMENSION(KIJS:KIJL), INTENT(OUT) :: MIJ

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NFRE), INTENT(IN) :: WAVNUM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NFRE), INTENT(IN) :: CINV
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NFRE), INTENT(IN) :: CGROUP
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NFRE), INTENT(IN) :: XK2CG
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NFRE), INTENT(IN) :: STOKFAC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: IBRMEM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: WSWAVE, WDWAVE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: AIRD, WSTAR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: USTRA, VSTRA
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: CICOVER
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: UFRIC, Z0M, Z0B, CHRNCK, CITHICK
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUOC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: TAUICX, TAUICY
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: PHIOCD, PHIEPS, PHIAW, USTOKES, VSTOKES
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: WSEMEAN, WSFMEAN, STRNMS
      REAL(KIND=JWRO), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN
      REAL(KIND=JWRO), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: NPHIEPS, NTAUOC, NSWH, NMWP, NEMOTAUX
      REAL(KIND=JWRO), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: NEMOTAUY, NEMOWSWAVE, NEMOPHIF
      REAL(KIND=JWRO), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: NEMOTAUICX, NEMOTAUICY
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      INTEGER(KIND=JWIM) :: IJ, K, M
      INTEGER(KIND=JWIM) :: ICALL, NCALL

      REAL(KIND=JWRB) :: DELT, GTEMP1, XIMP, DELT5
      REAL(KIND=JWRB) :: TAU, XN, PHIDIAG, TAUO, BETA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TAUW_LOC  ! TAUW should not be updated do use a local array
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TAUWDIR_LOC  ! TAUW should not be updated do use a local array
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: RAORW
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EMEAN, FMEAN, HALP
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EMEANWS, FMEANWS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: F1MEAN, AKMEAN, XKMEAN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: PHIWA
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: COSWDIF, SINWDIF2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: FLM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: RHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: FLD, SL, SPOS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: SLICE, SLTEMP

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ALPFAC

      LOGICAL :: LCFLX
      LOGICAL :: LUPDTUS

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WDFLUXES',0,ZHOOK_HANDLE)

!*    1. INITIALISATION.
!        ---------------

      DELT = IDELT
      XIMP = 1.0_JWRB
      DELT5 = XIMP*DELT

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

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            SLICE(IJ,K,M)  = 0.0_JWRB
            SLTEMP(IJ,K,M) = 0.0_JWRB
          ENDDO
        ENDDO
      ENDDO

      NCALL = 1
      ICALL = 1
      CALL SINFLX (ICALL, NCALL, KIJS, KIJL,        &
     &             LUPDTUS,                         &
     &             FL1,                             &
     &             WAVNUM, CINV, XK2CG,             &
     &             WSWAVE, WDWAVE, AIRD,            &
     &             RAORW, WSTAR, CICOVER,           &
     &             COSWDIF, SINWDIF2,               &
     &             FMEAN, HALP, FMEANWS,            &
     &             FLM,                             &
     &             UFRIC, TAUW_LOC, TAUWDIR_LOC,    &
     &             Z0M, Z0B, CHRNCK, PHIWA,         &
     &             FLD, SL, SPOS,                   &
     &             MIJ, RHOWGDFTH, XLLWS)

      IF (LCFLX) THEN

        CALL SDISSIP (KIJS, KIJL, FL1 ,FLD, SL,     &
     &                WAVNUM, CGROUP, XK2CG,        &
     &                EMEAN, F1MEAN, XKMEAN,        &
     &                UFRIC, COSWDIF, RAORW) 

        IF ( LICERUN ) THEN      

!         Use linear scaling of ALL proceeding source terms under sea ice (this is a complete unknown)
          IF (LCISCAL) THEN
            DO M = 1,NFRE
              DO K = 1,NANG
                DO IJ = KIJS,KIJL
                  BETA        = 1._JWRB - CICOVER(IJ)
                  SL(IJ,K,M)  = BETA*SL(IJ,K,M)
                  FLD(IJ,K,M) = BETA*FLD(IJ,K,M)
                END DO
              END DO
            END DO
          ENDIF

!        Coupling of waves and sea ice (type 1): wave-induced sea ice break up + reduced attenuation
          IF(LWNEMOCOUIBR) THEN 
            CALL ICEBREAK_MODIFY_ATTENUATION (KIJS,KIJL,IBRMEM,ALPFAC)           
          ENDIF


!        Save source term contributions relevant for the calculation of ice fluxes
          IF(LWNEMOCOUWRS) THEN 
            DO M=1,NFRE
              DO K=1,NANG
                DO IJ=KIJS,KIJL
                  SLTEMP(IJ,K,M) = SL(IJ,K,M)
                ENDDO
              ENDDO
            ENDDO
         ENDIF

!        Attenuation of waves in ice
          IF(LCIWA1 .OR. LCIWA2 .OR. LCIWA3) THEN
            CALL SDICE (KIJS, KIJL, FL1, FLD, SL, WAVNUM, CGROUP, CICOVER, CITHICK, ALPFAC)
         ENDIF
         

!        Save source term contributions relevant for the calculation of ice fluxes
         IF (LWNEMOCOUWRS) THEN
           DO M=1,NFRE
             DO K=1,NANG
               DO IJ=KIJS,KIJL
                 SLICE(IJ,K,M) = SL(IJ,K,M) - SLTEMP(IJ,K,M)
               ENDDO
             ENDDO
           ENDDO
         ENDIF
        
        ENDIF
            
        IF (.NOT. LWVFLX_SNL) THEN
          CALL WNFLUXES (KIJS, KIJL,                        &
     &                   MIJ, RHOWGDFTH,                    &
     &                   CINV,                              &
     &                   SL, SLICE, CICOVER,                &
     &                   PHIWA,                             &
     &                   EMEAN, F1MEAN, WSWAVE, WDWAVE,     &
     &                   USTRA, VSTRA,                      &
     &                   UFRIC, AIRD,                       &
     &                   NPHIEPS, NTAUOC, NSWH, NMWP,       &
     &                   NEMOTAUX, NEMOTAUY,                &
     &                   NEMOTAUICX, NEMOTAUICY,            &
     &                   NEMOWSWAVE, NEMOPHIF,              &
     &                   TAUXD, TAUYD, TAUOCXD, TAUOCYD,    &
     &                   TAUOC, TAUICX, TAUICY,             &
     &                   PHIOCD, PHIEPS, PHIAW,             &
     &                  .FALSE.)
        ENDIF

        CALL SNONLIN (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN)

        IF (LWVFLX_SNL) THEN
          CALL WNFLUXES (KIJS, KIJL,                        &
     &                   MIJ, RHOWGDFTH,                    &
     &                   CINV,                              &
     &                   SL, SLICE, CICOVER,                &
     &                   PHIWA,                             &
     &                   EMEAN, F1MEAN, WSWAVE, WDWAVE,     &
     &                   USTRA, VSTRA,                      &
     &                   UFRIC, AIRD,                       &
     &                   NPHIEPS, NTAUOC, NSWH, NMWP,       &
     &                   NEMOTAUX, NEMOTAUY,                &
     &                   NEMOTAUICX, NEMOTAUICY,            &
     &                   NEMOWSWAVE, NEMOPHIF,              &
     &                   TAUXD, TAUYD, TAUOCXD, TAUOCYD,    &
     &                   TAUOC, TAUICX, TAUICY,             &
     &                   PHIOCD, PHIEPS, PHIAW,             &
     &                  .FALSE.)
        ENDIF

        IF (LWFLUX) THEN
         CALL FEMEANWS(KIJS, KIJL, FL1, XLLWS, FMEANWS, EMEANWS)

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

        CALL STOKESTRN(KIJS, KIJL, FL1, WAVNUM, STOKFAC, DEPTH, WSWAVE, WDWAVE, CICOVER, CITHICK, &
&                      USTOKES, VSTOKES, STRNMS, NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN)

      ENDIF
! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WDFLUXES',1,ZHOOK_HANDLE)

END SUBROUTINE WDFLUXES
