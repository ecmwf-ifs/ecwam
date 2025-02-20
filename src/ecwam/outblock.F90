! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE OUTBLOCK (KIJS, KIJL, MIJ,                 &
 &                   FL1, XLLWS,                      & 
 &                   WAVNUM, CINV, CGROUP,            &
 &                   DEPTH, UCUR, VCUR, IODP,IBRMEM,  &
 &                   ALTWH, CALTWH, RALTCOR,          &
 &                   USTOKES, VSTOKES, STRNMS,        &
 &                   TAUXD, TAUYD, TAUOCXD,           &
 &                   TAUOCYD, TAUOC,                  &
 &                   TAUICX, TAUICY, PHIOCD,          &
 &                   PHIEPS, PHIAW,                   &
 &                   AIRD, WDWAVE, CICOVER,           &
 &                   WSWAVE, WSTAR,                   &
 &                   UFRIC, TAUW,                     &
 &                   Z0M, Z0B, CHRNCK,                &
 &                   CITHICK,                         &
 &                   NEMOCICOVER,                     &
 &                   NEMOCITHICK, NEMOUCUR, NEMOVCUR, &
 &                   BOUT)

! ----------------------------------------------------------------------

!**** *OUTBLOCK* - MODEL OUTPUT

!*    PURPOSE.
!     --------

!       PREPARE OUTPUT INTEGRATED PARAMETERS.


!**   INTERFACE.
!     ----------
!      *CALL*OUTBLOCK (KIJS, KIJL, MIJ,
!     &                FL1, XLLWS,
!     &                WVPRPT,
!     &                WVENVI, FF_NOW, INTFLDS, NEMO2WAM,
!     &                BOUT)
!      *KIJS*    - INDEX OF FIRST LOCAL GRIDPOINT
!      *KIJL*    - INDEX OF LAST LOCAL GRIDPOINT
!      *MIJ*     - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!      *FL1*     - INPUT SPECTRUM.
!      *XLLWS*   - WINDSEA MASK FROM INPUT SOURCE TERM
!      *WVPRPT*  - WAVE PROPERTIES
!      *WVENVI*  - WAVE ENVORONEMENT
!      *FF_NOW*  - FORCING FIELDS
!      *INTFLDS* - INTEGRATED/DERIVED PARAMETERS
!      *NEMO2WAM*- FIELDS FRON OCEAN MODEL to WAM
!      *BOUT*    - OUTPUT PARAMETER BUFFER

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY, FORCING_FIELDS,  &
     &                         INTGT_PARAM_FIELDS, OCEAN2WAVE

      USE YOWCOUT  , ONLY : NTRAIN   ,IPFGTBL  ,LSECONDORDER,           &
     &            NIPRMOUT, ITOBOUT  ,NTEWH    ,IPRMINFO
      USE YOWCOUP  , ONLY : LWNEMOCOUSTRN
      USE YOWFRED  , ONLY : FR, TH , DFIM, DELTH, COSTH, SINTH, XKMSS_CUTOFF
      USE YOWICE   , ONLY : FLMIN    ,LICERUN  ,LMASKICE
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZMISS    ,DEG      ,EPSUS    ,EPSU10, G, ZPI
      USE YOWSTAT  , ONLY : IREFRA

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "cal_second_order_spec.intfb.h"
#include "cimsstrn.intfb.h"
#include "ctcor.intfb.h"
#include "femean.intfb.h"
#include "intpol.intfb.h"
#include "kurtosis.intfb.h"
#include "meansqs.intfb.h"
#include "mwp1.intfb.h"
#include "mwp2.intfb.h"
#include "outbeta.intfb.h"
#include "outsetwmask.intfb.h"
#include "dominant_period.intfb.h"
#include "se10mean.intfb.h"
#include "sebtmean.intfb.h"
#include "sepwisw.intfb.h"
#include "sthq.intfb.h"
#include "wdirspread.intfb.h"
#include "weflux.intfb.h"
#include "w_maxh.intfb.h"
#include "ibrmemout.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), DIMENSION(KIJL), INTENT(IN) :: MIJ
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FL1, XLLWS
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: WAVNUM
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: CINV
      REAL(KIND=JWRB), DIMENSION(KIJL, NFRE), INTENT(IN) :: CGROUP

      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: DEPTH
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: UCUR
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: VCUR
      INTEGER(KIND=JWIM), DIMENSION(KIJL), INTENT(IN) :: IODP
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: IBRMEM

      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: AIRD, WDWAVE, CICOVER, WSWAVE, WSTAR
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: UFRIC, TAUW, Z0M, Z0B, CHRNCK, CITHICK
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: ALTWH, CALTWH, RALTCOR, USTOKES, VSTOKES, STRNMS
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUOC, PHIOCD
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: TAUICX, TAUICY
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: PHIEPS, PHIAW
      REAL(KIND=JWRO), DIMENSION(KIJL), INTENT(IN) :: NEMOCICOVER, NEMOCITHICK, NEMOUCUR, NEMOVCUR
      REAL(KIND=JWRB), DIMENSION(KIJL,NIPRMOUT), INTENT(OUT) :: BOUT


      INTEGER(KIND=JWIM) :: IJ, K, M, ITG, ITR, IH
      INTEGER(KIND=JWIM) :: IRA
      
      REAL(KIND=JWRB) :: SIG
      REAL(KIND=JWRB) :: GOZPI 
      REAL(KIND=JWRB) :: XMODEL_CUTOFF
      REAL(KIND=JWRB) :: TEWHMIN, TEWHMAX
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJL) :: EM, FM, DP
      REAL(KIND=JWRB), DIMENSION(KIJL) :: C3, C4, BF, QP, HMAX, TMAX
      REAL(KIND=JWRB), DIMENSION(KIJL) :: CMAX_F, HMAX_N, CMAX_ST, HMAX_ST, PHIST
      REAL(KIND=JWRB), DIMENSION(KIJL) :: ETA_M, R, XNSLC, SIG_TH, EPS, XNU
      REAL(KIND=JWRB), DIMENSION(KIJL) :: FLD1, FLD2
      REAL(KIND=JWRB), DIMENSION(KIJL) :: ESWELL ,FSWELL ,THSWELL, P1SWELL, P2SWELL, SPRDSWELL
      REAL(KIND=JWRB), DIMENSION(KIJL) :: ESEA   ,FSEA   ,THWISEA, P1SEA  , P2SEA  , SPRDSEA
      REAL(KIND=JWRB), DIMENSION(KIJL) :: CHARNOCK, BETAHQ, CDATM
      REAL(KIND=JWRB), DIMENSION(KIJL) :: HALP
      REAL(KIND=JWRB), DIMENSION(KIJL) :: ZTHRS, ZRDUC 
      REAL(KIND=JWRB), DIMENSION(KIJL,NTRAIN) :: EMTRAIN
      REAL(KIND=JWRB), DIMENSION(KIJL,NTRAIN) :: THTRAIN, PMTRAIN
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG) :: COSWDIF

!     *FL2ND*  SPECTRUM with second order effect added if LSECONDORDER is true .
!            and in the absolute frame of reference if currents are used
!            and will have low frequency noise added if waves in sea-ice
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE) :: FL2ND

      LOGICAL :: LLPEAKF

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('OUTBLOCK',0,ZHOOK_HANDLE)

!
!*    1. COMPUTE MEAN PARAMETERS.
!        ------------------------

!     PREPARE THE WAVE SPECTRA THAT SHOULD BE USED FOR OUTPUT

      LLPEAKF = .FALSE.

      IRA=1
      SIG = 1._JWRB
      GOZPI=G/ZPI

      IF (IREFRA == 2 .OR. IREFRA == 3) THEN
        CALL INTPOL (KIJS, KIJL, FL1, FL2ND, WAVNUM, UCUR, VCUR, IRA)
      ELSE
        FL2ND(KIJS:KIJL,:,:) = FL1(KIJS:KIJL,:,:)
      ENDIF
      IF (LSECONDORDER) CALL CAL_SECOND_ORDER_SPEC(KIJS, KIJL, FL2ND, WAVNUM, DEPTH, SIG)

      ! Adapting the noise level structure to be more consistent in sea ice conditions
      IF (LICERUN .AND. .NOT. LMASKICE) THEN
        DO IJ=KIJS,KIJL
          ZTHRS(IJ) = (1._JWRB - 0.9_JWRB*MIN(CICOVER(IJ),0.99_JWRB))*FLMIN
        ENDDO

        DO M=1,NFRE
          DO IJ=KIJS,KIJL
            ZRDUC(IJ) = EXP(-10.0_JWRB*FR(M)**2/SQRT(MAX(WSWAVE(IJ),1.0_JWRB)))
          ENDDO

          DO K=1,NANG
            DO IJ=KIJS,KIJL
              IF (FL2ND(IJ,K,M) <= ZTHRS(IJ)) THEN
                FL2ND(IJ,K,M) = MAX(ZRDUC(IJ) * FL2ND(IJ,K,M), ZTHRS(IJ)*ZRDUC(IJ)**2)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!     COMPUTE MEAN PARAMETERS

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          COSWDIF(IJ,K) = COS(TH(K)-WDWAVE(IJ))
        ENDDO
      ENDDO

      CALL FEMEAN (KIJS, KIJL, FL2ND, EM, FM)

      CALL DOMINANT_PERIOD (KIJS, KIJL, FL2ND, DP)

      CALL KURTOSIS(KIJS, KIJL, FL1,                     &
     &              DEPTH,                               &
     &              C3, C4, BF, QP, HMAX, TMAX,          &
     &              ETA_M, R, XNSLC, SIG_TH, EPS, XNU)

!     WIND/SWELL PARAMETERS
      CALL SEPWISW (KIJS, KIJL, MIJ, FL1, XLLWS, CINV,             &
     &              UFRIC, WSWAVE, WDWAVE, COSWDIF,  &
     &              ESWELL, FSWELL, THSWELL, P1SWELL, P2SWELL, SPRDSWELL, &
     &              ESEA, FSEA, THWISEA, P1SEA, P2SEA, SPRDSEA,           &
     &              EMTRAIN, THTRAIN, PMTRAIN)


!     LOAD THE OUTPUT BUFFER:

      IF (IPFGTBL(1) /= 0) THEN
!       SIGNIFICANT WAVE HEIGHT CONVERSION
        BOUT(KIJS:KIJL,ITOBOUT(1))=4._JWRB*SQRT(MAX(EM(KIJS:KIJL),0._JWRB))
      ENDIF

      IF (IPFGTBL(2) /= 0) THEN
        ITG=ITOBOUT(2)
        CALL STHQ (KIJS, KIJL, FL2ND, BOUT(KIJS,ITG))
!       CONVERT DIRECTIONS TO DEGREES AND METEOROLOGICAL CONVENTION
        BOUT(KIJS:KIJL,ITG)=MOD(DEG*BOUT(KIJS:KIJL,ITG)+180._JWRB,360._JWRB)
      ENDIF

      IF (IPFGTBL(3) /= 0) THEN
!       CONVERSION TO PERIOD
        DO IJ=KIJS,KIJL
          IF (FM(IJ) > 0._JWRB) THEN
            BOUT(IJ,ITOBOUT(3))=1._JWRB/FM(IJ)
          ELSE
            BOUT(IJ,ITOBOUT(3))=ZMISS
          ENDIF
        ENDDO
      ENDIF

      IF (IPFGTBL(4) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(4))=UFRIC(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(5) /= 0) THEN
!       CONVERT DIRECTIONS TO DEGREES AND METEOROLOGICAL CONVENTION
        BOUT(KIJS:KIJL,ITOBOUT(5))=MOD(DEG*WDWAVE(KIJS:KIJL)+180._JWRB,360._JWRB)
      ENDIF

      IF (IPFGTBL(6) /= 0) THEN
!       CONVERSION TO PERIOD
        DO IJ=KIJS,KIJL
          IF (DP(IJ) > 0.0_JWRB) THEN
            BOUT(IJ,ITOBOUT(6))=DP(IJ)
          ELSE
            BOUT(IJ,ITOBOUT(6))=ZMISS
          ENDIF
        ENDDO
      ENDIF

      IF (IPFGTBL(7) /= 0) THEN
!!      if the numerical computation of TAU and CD changes, a similar
!!      modification has to be put in buildstress where the friction
!!      velocity is determined from U10 and CD.
!!      Because of the limited numerical resolution when encoding in grib, the maximum value for Cd is set to 0.01
!!!!!!!!!!        BOUT(KIJS:KIJL,ITOBOUT(IR))=MIN(MAX(UFRIC(KIJS:KIJL)**2,EPSUS)/MAX(WSWAVE(KIJS:KIJL)**2,EPSU10**2), 0.01_JWRB)
!!! output the drag coeffient that is consistent with the Charnock parameter that is returned to the atmosphere model:
        CALL OUTBETA (KIJS, KIJL,                      &
     &                WSWAVE, UFRIC, Z0M, Z0B, CHRNCK, &
     &                CHARNOCK, BETAHQ, CD=CDATM)

        BOUT(KIJS:KIJL,ITOBOUT(7))=MIN(CDATM(KIJS:KIJL), 0.01_JWRB)
      ENDIF

      IF (IPFGTBL(8) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(8))=TAUW(KIJS:KIJL)/MAX(UFRIC(KIJS:KIJL)**2,EPSUS)
      ENDIF

      IF (IPFGTBL(9) /= 0) THEN
        CALL MEANSQS (XKMSS_CUTOFF, KIJS, KIJL, FL1, WAVNUM, UFRIC, COSWDIF, BOUT(:,ITOBOUT(9)))
      ENDIF

      IF (IPFGTBL(10) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(10))=WSWAVE(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(11) /= 0) THEN
!       WINDSEA SIGNIFICANT WAVE HEIGHT CONVERSION
        BOUT(KIJS:KIJL,ITOBOUT(11))=4._JWRB*SQRT(MAX(ESEA(KIJS:KIJL),0._JWRB))
      ENDIF

      IF (IPFGTBL(12) /= 0) THEN
!       TOTAL SWELL SIGNIFICANT WAVE HEIGHT CONVERSION
        BOUT(KIJS:KIJL,ITOBOUT(12))=4._JWRB*SQRT(MAX(ESWELL(KIJS:KIJL),0._JWRB))
      ENDIF

      IF (IPFGTBL(13) /= 0) THEN
!       CONVERT DIRECTIONS TO DEGREES AND METEOROLOGICAL CONVENTION
        BOUT(KIJS:KIJL,ITOBOUT(13))=MOD(DEG*THWISEA(KIJS:KIJL)+180._JWRB,360._JWRB)
      ENDIF

      IF (IPFGTBL(14) /= 0) THEN
!       CONVERT DIRECTIONS TO DEGREES AND METEOROLOGICAL CONVENTION
        BOUT(KIJS:KIJL,ITOBOUT(14))=MOD(DEG*THSWELL(KIJS:KIJL)+180._JWRB,360._JWRB)
      ENDIF

      IF (IPFGTBL(15) /= 0) THEN
!       CONVERSION TO PERIOD
        DO IJ=KIJS,KIJL
          IF (FSEA(IJ) > 0._JWRB) THEN
            BOUT(IJ,ITOBOUT(15))=1._JWRB/FSEA(IJ)
          ELSE
            BOUT(IJ,ITOBOUT(15))=ZMISS
          ENDIF
        ENDDO
      ENDIF

      IF (IPFGTBL(16) /= 0) THEN
!       CONVERSION TO PERIOD
        DO IJ=KIJS,KIJL
          IF (FSWELL(IJ) > 0._JWRB) THEN
            BOUT(IJ,ITOBOUT(16))=1._JWRB/FSWELL(IJ)
          ELSE
            BOUT(IJ,ITOBOUT(16))=ZMISS
          ENDIF
        ENDDO
      ENDIF

      IF (IPFGTBL(17) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(17))=ALTWH(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(18) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(18))=CALTWH(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(19) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(19))=RALTCOR(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(20) /= 0) THEN
        CALL MWP1 (KIJS, KIJL, FL2ND, BOUT(:,ITOBOUT(20)))
      ENDIF

      IF (IPFGTBL(21) /= 0) THEN
        CALL MWP2 (KIJS, KIJL, FL2ND, BOUT(:,ITOBOUT(21)))
      ENDIF

      IF (IPFGTBL(22) /= 0) THEN
        CALL WDIRSPREAD (KIJS, KIJL, FL2ND, EM, LLPEAKF, BOUT(:,ITOBOUT(22)))
      ENDIF

      IF (IPFGTBL(23) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(23))=P1SEA(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(24) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(24))=P1SWELL(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(25) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(25))=P2SEA(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(26) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(26))=P2SWELL(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(27) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(27))=SPRDSEA(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(28) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(28))=SPRDSWELL(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(29) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(29))=C4(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(30) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(30))=BF(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(31) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(31))=QP(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(32) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(32))=DEPTH(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(33) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(33))=HMAX(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(34) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(34))=TMAX(KIJS:KIJL)
      ENDIF

!     SURFACE STOKES DRIFT U and V
      IF (IPFGTBL(35) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(35))=USTOKES(KIJS:KIJL)
      ENDIF
      IF (IPFGTBL(36) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(36))=VSTOKES(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(37) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(37))=UCUR(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(38) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(38))=VCUR(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(39) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(39))=PHIEPS(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(40) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(40))=PHIAW(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(41) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(41))=TAUOC(KIJS:KIJL)
      ENDIF

      DO ITR=1,NTRAIN
        IF (IPFGTBL(42 + (ITR-1)*3) /= 0) THEN
          BOUT(KIJS:KIJL,ITOBOUT(42 + (ITR-1)*3))=4._JWRB*SQRT(MAX(EMTRAIN(KIJS:KIJL,ITR),0._JWRB))
        ENDIF

        IF (IPFGTBL(43 + (ITR-1)*3) /= 0) THEN
          BOUT(KIJS:KIJL,ITOBOUT(43 + (ITR-1)*3))=MOD(DEG*THTRAIN(KIJS:KIJL,ITR)+180._JWRB,360._JWRB)
        ENDIF

        IF (IPFGTBL(44 + (ITR-1)*3) /= 0) THEN
          BOUT(KIJS:KIJL,ITOBOUT(44 + (ITR-1)*3))=PMTRAIN(KIJS:KIJL,ITR)
        ENDIF
      ENDDO

      IF (IPFGTBL(42 + 3*NTRAIN) /= 0) THEN
        IF (LWNEMOCOUSTRN) THEN
          BOUT(KIJS:KIJL,ITOBOUT(42 + 3*NTRAIN))=STRNMS(KIJS:KIJL)
        ELSE
          CALL CIMSSTRN (KIJS, KIJL, FL1, WAVNUM, DEPTH, CITHICK, BOUT(:,ITOBOUT(42 + 3*NTRAIN)))
        ENDIF
      ENDIF

      IF (IPFGTBL(43 + 3*NTRAIN) /= 0) THEN
        CALL SE10MEAN (KIJS, KIJL, FL2ND, FLD1)
        BOUT(KIJS:KIJL,ITOBOUT(43 + 3*NTRAIN))=4._JWRB*SQRT(MAX(FLD1(KIJS:KIJL),0._JWRB))
      ENDIF

      IF (IPFGTBL(44 + 3*NTRAIN) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(44 + 3*NTRAIN))=AIRD(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(45 + 3*NTRAIN) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(45 + 3*NTRAIN))=WSTAR(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(46 + 3*NTRAIN) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(46 + 3*NTRAIN))=CICOVER(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(47 + 3*NTRAIN) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(47 + 3*NTRAIN))=CITHICK(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(48 + 3*NTRAIN) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(48 + 3*NTRAIN))=C3(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(49 + 3*NTRAIN) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(49 + 3*NTRAIN))=NEMOCICOVER(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(50 + 3*NTRAIN) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(50 + 3*NTRAIN))=NEMOCITHICK(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(51 + 3*NTRAIN) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(51 + 3*NTRAIN))=NEMOUCUR(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(52 + 3*NTRAIN) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(52 + 3*NTRAIN))=NEMOVCUR(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(53 + 3*NTRAIN) /= 0 .OR. IPFGTBL(54 + 3*NTRAIN) /= 0) THEN
           CALL WEFLUX (KIJS, KIJL, FL1, CGROUP,    &
     &                  NFRE, NANG, DFIM, DELTH,    &
     &                  COSTH, SINTH,               &
     &                  FLD1, FLD2)
      ENDIF
      IF (IPFGTBL(53 + 3*NTRAIN) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(53 + 3*NTRAIN))=FLD1(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(54 + 3*NTRAIN) /= 0) THEN
!       CONVERT DIRECTIONS TO DEGREES AND METEOROLOGICAL CONVENTION
        BOUT(KIJS:KIJL,ITOBOUT(54 + 3*NTRAIN))=MOD(DEG*FLD2(KIJS:KIJL)+180._JWRB,360._JWRB)
      ENDIF

      DO IH=1,NTEWH
        IF (IPFGTBL(54 + 3*NTRAIN + IH) /= 0) THEN
          TEWHMIN = REAL(IPRMINFO(54 + 3*NTRAIN + IH,4),JWRB)
          TEWHMAX = REAL(IPRMINFO(54 + 3*NTRAIN + IH,5),JWRB) 
          CALL SEBTMEAN (KIJS, KIJL, FL2ND, TEWHMIN, TEWHMAX, BOUT(:,ITOBOUT(54 + 3*NTRAIN + IH)))
!         SIGNIFICANT WAVE HEIGHT CONVERSION
          BOUT(KIJS:KIJL,ITOBOUT(54 + 3*NTRAIN + IH))=4._JWRB*SQRT(MAX(BOUT(KIJS:KIJL,ITOBOUT(54 + 3*NTRAIN + IH)),0._JWRB))
        ENDIF
      ENDDO

      IF (IPFGTBL(55 + 3*NTRAIN + NTEWH) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(55 + 3*NTRAIN + NTEWH))=ETA_M(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(56 + 3*NTRAIN + NTEWH) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(56 + 3*NTRAIN + NTEWH))=R(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(57 + 3*NTRAIN + NTEWH) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(57 + 3*NTRAIN + NTEWH))=XNSLC(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(58 + 3*NTRAIN + NTEWH) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(58 + 3*NTRAIN + NTEWH))=TAUXD(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(59 + 3*NTRAIN + NTEWH) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(59 + 3*NTRAIN + NTEWH))=TAUYD(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(60 + 3*NTRAIN + NTEWH) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(60 + 3*NTRAIN + NTEWH))=TAUOCXD(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(61 + 3*NTRAIN + NTEWH) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(61 + 3*NTRAIN + NTEWH))=TAUOCYD(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(62 + 3*NTRAIN + NTEWH) /= 0) THEN
!      !!! make the energy flux positive
        BOUT(KIJS:KIJL,ITOBOUT(62 + 3*NTRAIN + NTEWH))=MAX(-PHIOCD(KIJS:KIJL),0.0_JWRB)
      ENDIF

!!    alternative ways to determine wave height extremes 
      IF (IPFGTBL(63 + 3*NTRAIN + NTEWH) /= 0 .OR. IPFGTBL(64 + 3*NTRAIN + NTEWH) /= 0  .OR. &
&         IPFGTBL(65 + 3*NTRAIN + NTEWH) /= 0 .OR. IPFGTBL(66 + 3*NTRAIN + NTEWH) /= 0 ) THEN
        CALL W_MAXH (KIJS, KIJL, FL1, DEPTH, WAVNUM,            &
     &               CMAX_F, HMAX_N, CMAX_ST, HMAX_ST, PHIST)
      ENDIF 

      IF (IPFGTBL(63 + 3*NTRAIN + NTEWH) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(63 + 3*NTRAIN + NTEWH))=CMAX_F(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(64 + 3*NTRAIN + NTEWH) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(64 + 3*NTRAIN + NTEWH))=HMAX_N(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(65 + 3*NTRAIN + NTEWH) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(65 + 3*NTRAIN + NTEWH))=CMAX_ST(KIJS:KIJL)
      ENDIF

      IF (IPFGTBL(66 + 3*NTRAIN + NTEWH) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(66 + 3*NTRAIN + NTEWH))=HMAX_ST(KIJS:KIJL)
      ENDIF
!!


!     COMPUTE OUTPUT EXTRA FIELDS
!     add necessary code to compute the extra output fields
!!!for testing
      IF (IPFGTBL(67 + 3*NTRAIN + NTEWH) /= 0) THEN
        CALL CTCOR (KIJS, KIJL, FL1, BOUT(:,ITOBOUT(67 + 3*NTRAIN + NTEWH)))
      ENDIF

      IF (IPFGTBL(68 + 3*NTRAIN + NTEWH) /= 0) THEN
        XMODEL_CUTOFF=(ZPI*FR(NFRE))**2/G
        CALL MEANSQS (XMODEL_CUTOFF, KIJS, KIJL, FL1, WAVNUM, UFRIC, COSWDIF, BOUT(:,ITOBOUT(68 + 3*NTRAIN + NTEWH)))
      ENDIF

      IF (IPFGTBL(69 + 3*NTRAIN + NTEWH) /= 0) THEN
        CALL IBRMEMOUT (KIJS, KIJL, IBRMEM, CICOVER, BOUT(:,ITOBOUT(69 + 3*NTRAIN + NTEWH)))
      ENDIF

      IF (IPFGTBL(70 + 3*NTRAIN + NTEWH) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(70 + 3*NTRAIN + NTEWH))=MAX(TAUICX(KIJS:KIJL),0.0_JWRB)
      ENDIF

      IF (IPFGTBL(71 + 3*NTRAIN + NTEWH) /= 0) THEN
        BOUT(KIJS:KIJL,ITOBOUT(71 + 3*NTRAIN + NTEWH))=MAX(TAUICY(KIJS:KIJL),0.0_JWRB)
      ENDIF



!     APPLY SEA ICE MASK AND SEA MASK IF NECESSARY
      CALL OUTSETWMASK (KIJS, KIJL, IODP, CICOVER, BOUT)

IF (LHOOK) CALL DR_HOOK('OUTBLOCK',1,ZHOOK_HANDLE)

END SUBROUTINE OUTBLOCK
