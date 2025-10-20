! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE SINFLX_ZBRY (ICALL, NCALL, NGST, KIJS, KIJL,  &
 &                 LUPDTUS,                   &
 &                 FL1,                       &
 &                 WAVNUM,CGROUP, CINV,       &
 &                 WSWAVE, WDWAVE, AIRD,      &
 &                 RAORW,  WSTAR, CICOVER,    &
 &                 COSWDIF, SINWDIF2,         &
 &                 FMEAN, HALP, FMEANWS,      &
 &                 FLM,                       &
 &                 UFRIC, UPROXY, TAUW, TAUWDIR,      & 
 &                 Z0M, Z0B, CHRNCK, PHIWA,   &
 &                 FLD, SL, SPOS,             &
 &                 MIJ, RHOWGDFTH, XLLWS)

! ----------------------------------------------------------------------

!**** *SINFLX_ZBRY* - COMPUTATION OF INPUT SOURCE FUNCTION AND STRESSES
!
!     JOSH KOUSAL & JEAN BIDLOT    ECMWF 2023
!
!*    PURPOSE.
!     ---------

!      Observation-based source term for wind input after Donelan, Babanin,
!      Young and Banner (Donelan et al ,2006) following the implementation
!      by Rogers et al. (2012).
!
!**   INTERFACE.
!     ----------

!     *CALL* *SINFLX_ZBRY (NGST, LLSNEG, KIJS, KIJL, FL1,
!    &                    WAVNUM, CGROUP, CINV,
!    &                    WSWAVE, WDWAVE, UFRIC, Z0M,
!    &                    COSWDIF, SINWDIF2,
!    &                    RAORW, WSTAR, RNFAC,
!    &                    FLD, SL, SPOS, XLLWS)
!         *NGST* - IF = 1 THEN NO GUSTINESS PARAMETERISATION
!                - IF = 2 THEN GUSTINESS PARAMETERISATION
!         *LLSNEG- IF TRUE THEN THE NEGATIVE SINPUT (SWELL DAMPING) WILL BE COMPUTED
!         *KIJS* - INDEX OF FIRST GRIDPOINT.
!         *KIJL* - INDEX OF LAST GRIDPOINT.
!          *FL1* - SPECTRUM.
!       *WAVNUM* - WAVE NUMBER.
!       *CGROUP* - GROUP SPEED
!         *CINV* - INVERSE PHASE VELOCITY.
!       *WDWAVE* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!        *UFRIC* - NEW FRICTION VELOCITY IN M/S.
!        *Z0M* - ROUGHNESS LENGTH IN M.
!      *COSWDIF* - COS(TH(K)-WDWAVE(IJ))
!     *SINWDIF2* - SIN(TH(K)-WDWAVE(IJ))**2
!        *RAORW* - RATIO AIR DENSITY TO WATER DENSITY.
!        *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
!        *RNFAC* - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
!        *CHRNCK*- CHARNOCK COEFFICIENT 
!          *FLD* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!           *SL* - TOTAL SOURCE FUNCTION ARRAY.
!         *SPOS* - POSITIVE SOURCE FUNCTION ARRAY.
!       *XLLWS*  - = 1 WHERE SINPUT IS POSITIVE

!     METHOD.
!     -------

!       SEE REFERENCE.

!     EXTERNALS.
!     ----------
!     TAU_WAVE_ATMOS
!     LFACTOR
!     IRANGE

!     ORIGIN.
!     ----------
!     Adapted from Babanin Young Donelan & Banner (ZBRY) physics 
!     as implemented as ST6 in WAVEWATCH-III
!     WW3 module:       W3SRC6MD
!     WW3 subroutine:   W3SIN6
!     Implementation into ECWAM DECEMBER 2021 by J. Kousal


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU    ,LLCAPCHNK

      USE YOWWNDG  , ONLY : ICODE    ,ICODE_CPL

      USE YOWFRED  , ONLY : FR       ,TH       ,DFIM     ,COSTH  ,SINTH, ZPIFR, DELTH, FRATIO, FRIC
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,GM1      ,EPSMIN, EPSUS, ZPI, ROWATER
      USE YOWPHYS  , ONLY : ZALP     ,TAUWSHELTER, XKAPPA, BETAMAXOXKAPPA2,    &
     &                      RNU      ,RNUM, &
     &                      SWELLF   ,SWELLF2  ,SWELLF3  ,SWELLF4  , SWELLF5, &
     &                      SWELLF6  ,SWELLF7  ,SWELLF7M1, Z0RAT   ,Z0TUBMAX , &
     &                      ABMIN  ,ABMAX, CDFAC, DTHRN_A  ,DTHRN_U, RNU_WATER, &
     &                      ZSIN6A0
      USE YOWTEST  , ONLY : IU06
      USE YOWTABL  , ONLY : IAB      ,SWELLFT
      USE YOWSTAT  , ONLY : IPHYS2_AIRSEA, LLLOWWINDS

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "airsea.intfb.h"
#include "femeanws.intfb.h"
#include "frcutindex.intfb.h"
#include "halphap.intfb.h"
#include "wsigstar.intfb.h"
#include "tau_wave_atmos.intfb.h"
#include "lfactor.intfb.h"
#include "irange.intfb.h"
#include "calcphiwa.intfb.h"


INTEGER(KIND=JWIM), INTENT(IN) :: ICALL  !! CALL NUMBER.
INTEGER(KIND=JWIM), INTENT(IN) :: NCALL  !! TOTAL NUMBER OF CALLS.
INTEGER(KIND=JWIM), INTENT(IN) :: NGST  !! GUSTINESS PARAMETERIZATION
INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL !! GRID POINT INDEXES.  

LOGICAL, INTENT(IN) :: LUPDTUS  !! IF TRUE UFRIC AND Z0M WILL BE UPDATED (CALLING AIRSEA).

REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FL1  !! WAVE SPECTRUM.
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM  !! WAVE NUMBER.
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: CGROUP  !! GROUP VELOCITY.
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: CINV    !! INVERSE PHASE VELOCITY.
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: WSWAVE !! WIND SPEED IN M/S.
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: WDWAVE !! WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC NOTATION.
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: AIRD  !! AIR DENSITY (KG/M**3).
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: RAORW  !! RATIO AIR DENSITY TO WATER DENSITY.
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: WSTAR !! FREE CONVECTION VELOCITY SCALE (M/S)
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: CICOVER  !! SEA ICE COVER.
REAL(KIND=JWRB), DIMENSION(KIJL, NANG), INTENT(IN) :: COSWDIF  !! COS(TH(K)-WDWAVE(IJ))
REAL(KIND=JWRB), DIMENSION(KIJL, NANG), INTENT(IN) :: SINWDIF2 !! SIN(TH(K)-WDWAVE(IJ))**2
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: FMEAN  !! MEAN FREQUENCY.
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: HALP   !! 1/2 PHILLIPS PARAMETER  
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(OUT) :: FMEANWS  !! MEAN FREQUENCY OF THE WINDSEA.
REAL(KIND=JWRB), DIMENSION(KIJL,NANG), INTENT(IN) :: FLM  !! SPECTAL DENSITY MINIMUM VALUE
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: UFRIC !! FRICTION VELOCITY IN M/S.
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: UPROXY !! FRICTION VELOCITY IN M/S.
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: TAUW  !! WAVE STRESS IN (M/S)**2
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: TAUWDIR  !! WAVE STRESS DIRECTION.
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: Z0M  !! ROUGHNESS LENGTH IN M.
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: Z0B  !! BACKGROUND ROUGHNESS LENGTH.
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: CHRNCK  !! CHARNOCK COEFFICIENT.

REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(OUT) :: PHIWA  !! ENERGY FLUX FROM WIND INTO WAVES.
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(OUT) :: FLD !! DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(OUT) :: SL !! TOTAL SOURCE FUNCTION ARRAY.
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(OUT) :: SPOS !!  POSITIVE SINPUT ONLY.

INTEGER(KIND=JWIM), INTENT(OUT) :: MIJ(KIJL)  !! LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.

REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(OUT) :: RHOWGDFTH  !! WATER DENSITY * G * DF * DTHETA

REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(OUT) :: XLLWS  !! TOTAL WINDSEA MASK FROM INPUT SOURCE TERM.

INTEGER(KIND=JWIM) :: IUSFG, ICODE_WND

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JWRB), DIMENSION(KIJL) :: RNFAC

LOGICAL            :: LLFACT
INTEGER(KIND=JWIM) :: IJ, K, M, IND, IGST
INTEGER(KIND=JWIM) :: NSPEC !num. of freqs, dirs, spec. bins
INTEGER(KIND=JWIM), DIMENSION(NANG) :: ITHN
INTEGER(KIND=JWIM), DIMENSION(NFRE) :: IKN

REAL(KIND=JWRB), DIMENSION(NANG*NFRE)      :: ECOS2, ESIN2, SIG2
REAL(KIND=JWRB), DIMENSION(NFRE)           :: DSII, SIG, DF
REAL(KIND=JWRB), DIMENSION(NANG,NFRE) :: SPOSDENSIG, SNEGDENSIG

REAL(KIND=JWRB), DIMENSION(KIJL,NANG*NFRE) :: CG2, WN2
REAL(KIND=JWRB), DIMENSION(KIJL,NANG*NFRE) :: SQRTBN2, CINV2, A
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE)      :: ADENSIG, KMAX, ANAR, SQRTBN, CINV1
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE) :: KK
REAL(KIND=JWRB), DIMENSION(KIJL,NANG*NFRE,NGST) :: W1, W2, S, D
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE,NGST)      :: LFACT
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE,NGST) :: SDENSIG, DINPOS, DINTOT


! REAL(KIND=JWRB), PARAMETER :: SIN6A0 = 9.0E-2_JWRB ! ST6 PARAM ! TODO, move to PHYS
REAL(KIND=JWRB), DIMENSION(KIJL,NGST) :: TAUWX, TAUWY ! Component of the wave-supported stress
REAL(KIND=JWRB), DIMENSION(KIJL,NGST) :: TAUNWX, TAUNWY ! Component of the neg. wave-supported stress
REAL(KIND=JWRB), DIMENSION(KIJL) :: COSU, SINU
REAL(KIND=JWRB), DIMENSION(KIJL,NGST) :: UPROXYGST

REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: XK, CGG_WAM, CM
REAL(KIND=JWRB), DIMENSION(NFRE)      :: SIGM1

! For USTAR, Z0, CHNK
REAL(KIND=JWRB), DIMENSION(KIJL,NGST) :: TAU
REAL(KIND=JWRB), PARAMETER :: ZRN=1.65E-6_JWRB  ! effective kinematic viscosity (0.11*1.5e-5)
REAL(KIND=JWRB), PARAMETER :: RKAP = 0.4_JWRB
REAL(KIND=JWRB)            :: ZNLEV, Z0, KUOUST, USTM1, USTM2
REAL(KIND=JWRB), PARAMETER :: XEPS=0.00001_JWRB
REAL(KIND=JWRB), PARAMETER :: USTMIN=0.000001_JWRB
REAL(KIND=JWRB), PARAMETER :: PCHARMAX=0.1_JWRB
REAL(KIND=JWRB), PARAMETER :: Z0FG=0.01_JWRB
INTEGER(KIND=JWIM) :: ITER
REAL(KIND=JWRB)    :: XZNLEV, PCHAROG, XKUTOP, XOLOGZ0
REAL(KIND=JWRB)    :: UST, USTOLD, Z0CH, Z0VIS, Z0TOT, FF, DELF
REAL(KIND=JWRB)    :: CHARNOCK_MIN,CHNKMIN     ! For Capping
INTEGER(KIND=JWIM), PARAMETER :: NITER=15
REAL(KIND=JWRB), PARAMETER :: ALPHAMAX=0.1_JWRB
REAL(KIND=JWRB), PARAMETER :: AMAX=0.02_JWRB
REAL(KIND=JWRB), PARAMETER :: BMAX=0.01_JWRB
REAL(KIND=JWRB) :: ALPHAOGMAXU10

REAL(KIND=JWRB), DIMENSION(KIJL)      :: ROAIRN, CHNKOG, TAUNW

! For GUSTINESS
REAL(KIND=JWRB)                          :: AVG_GST
REAL(KIND=JWRB), DIMENSION(KIJL)      :: SIG_N, SIG_U10, TAUWGST_AVG, TAUWDIRGST_AVG, TAUNWGST_AVG, USTARGST_AVG, UPROXYGST_AVG
REAL(KIND=JWRB), DIMENSION(KIJL,NGST) :: TAUWGST, TAUWDIRGST, TAUNWGST, UABSGST, USTARGST, Z0GST
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE) :: SLGST_AVG, SPOSGST_AVG, FLGST_AVG
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE,NGST) :: SLGST, SPOSGST, FLGST

! For PHIWA calculation
! REAL(KIND=JWRB),DIMENSION(KIJL,NFRE) :: RHOWGDFTH
! REAL(KIND=JWRB), DIMENSION(KIJL) :: SUMT


! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SINFLX',0,ZHOOK_HANDLE)

! UPDATE UFRIC AND Z0M
IF (ICALL == 1 ) THEN
  IUSFG = 0
  IF (LWCOU) THEN
    ICODE_WND = ICODE_CPL
  ELSE
    ICODE_WND = ICODE
  ENDIF
ELSE
  IUSFG = 1
  ICODE_WND = 3
ENDIF

IF(LLCAPCHNK) THEN
  RNFAC(KIJS:KIJL) = 1.0_JWRB+DTHRN_A*(1.0_JWRB+TANH(WSWAVE(KIJS:KIJL)-DTHRN_U))
ELSE
  RNFAC(KIJS:KIJL) = 1.0_JWRB
ENDIF


IF(LUPDTUS) THEN
  ! increase noise level in the tail
  IF (ICALL == 1 ) THEN
    DO K=1,NANG
      FL1(KIJS:KIJL,K,NFRE) = MAX(FL1(KIJS:KIJL,K,NFRE),FLM(KIJS:KIJL,K))
    ENDDO
    HALP(KIJS:KIJL) = 0.0_JWRB
  ENDIF

  !$loki inline
  CALL AIRSEA (KIJS, KIJL,                                  &
&              HALP, WSWAVE, WDWAVE, TAUW, TAUWDIR, RNFAC,  &
&              UFRIC, Z0M, Z0B, CHRNCK, ICODE_WND, IUSFG) 

ENDIF


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! input source term!!!! (start)

NSPEC = NANG * NFRE   ! NUMBER OF SPECTRAL BINS

!     Wind height
ZNLEV    = 10._JWRB

!     COMPUTE FREQUENCY INTERVALLS
DO M = 1,NFRE
  DF(M) = FR(M)*( FRATIO - 1.0_JWRB )
ENDDO

DO M = 1,NFRE
  SIG(M)   = ZPI*FR(M)
  DSII(M)  = ZPI*DF(M)
  SIGM1(M) = 1.0_JWRB/SIG(M)
END DO

DO M=1,NFRE
  DO IJ=KIJS,KIJL
    CM(IJ,M)     = WAVNUM(IJ,M)*SIGM1(M)
  ENDDO
ENDDO


ITHN   = IRANGE(1,NANG,1)    ! Index vector 1:NANG
DO M = 1, NFRE
  ECOS2 (ITHN+(M-1)*NANG) = COSTH
  ESIN2 (ITHN+(M-1)*NANG) = SINTH
END DO
!
IKN    = IRANGE(1,NSPEC,NANG)   ! Index vector for elements of 1 ... NFRE
!                               ! such that e.g. SIG(1:NFRE) = SIG2(IKN).

DO K = 1, NANG                    ! Apply to all directions 
  SIG2  (IKN+(K-1)) = SIG
END DO

!     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
CALL WSIGSTAR (KIJS, KIJL, WSWAVE, UFRIC, Z0M, WSTAR, SIG_N, SIG_U10)
AVG_GST = 1.0_JWRB/NGST

IF (NGST == 1) THEN
  DO IJ=KIJS,KIJL
    USTARGST(IJ,1) = UFRIC(IJ)
    UABSGST(IJ,1)  = WSWAVE(IJ)
  ENDDO
ELSE IF (NGST == 2) THEN
    DO IJ=KIJS,KIJL
      USTARGST(IJ,1)= UFRIC(IJ)*(1.0_JWRB+SIG_N(IJ))
      USTARGST(IJ,2)= UFRIC(IJ)*(1.0_JWRB-SIG_N(IJ))
      UABSGST(IJ,1)= WSWAVE(IJ)*(1.0_JWRB+SIG_U10(IJ))
      UABSGST(IJ,2)= WSWAVE(IJ)*(1.0_JWRB-SIG_U10(IJ))
    END DO
ELSE
      WRITE (IU06,*) '**************************************'
      WRITE (IU06,*) '*    FATAL ERROR                     *'
      WRITE (IU06,*) '*    ===========                     *'
      WRITE (IU06,*) '* IN SINFLX_ZBRY: NGST > 2            *'
      WRITE (IU06,*) '* NGST = ', NGST
      WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.  *'
      WRITE (IU06,*) '*                                    *'
      WRITE (IU06,*) '**************************************'
      CALL ABORT1
ENDIF

DO IJ=KIJS,KIJL
  CHNKOG(IJ)    = CHRNCK(IJ)*GM1
  ROAIRN(IJ)    = RAORW(IJ)*ROWATER        
END DO

! Define Z0GST associated with USTARGST (as in airsea_zbry)
DO IGST=1,NGST
  DO IJ=KIJS,KIJL
    UST            = USTARGST(IJ,IGST)
    PCHAROG        = MIN(CHNKOG(IJ),PCHARMAX/G)
    Z0CH           = PCHAROG*UST**2
    Z0VIS          = ZRN/UST
    Z0GST(IJ,IGST) = Z0CH+Z0VIS
  ENDDO
ENDDO

!/  --- Main loop over LOC ----------------------------------- /

DO K = 1, NANG
  DO IJ = KIJS,KIJL
      WN2   (IJ,IKN+(K-1)) = WAVNUM(IJ,:)  ! using WAM native WN,CG
      CG2   (IJ,IKN+(K-1)) = CGROUP(IJ,:)
  END DO
END DO

DO IJ = KIJS,KIJL
  CINV2(IJ,:)  = WN2(IJ,:) / SIG2            ! inverse phase speed
  CINV1(IJ,:)  = CINV2(IJ,IKN)
END DO

!/ 0) --- set up a basic variables ----------------------------------- /
DO IJ = KIJS,KIJL
  COSU(IJ)  = COS(WDWAVE(IJ))
  SINU(IJ)  = SIN(WDWAVE(IJ))
END DO  
!
DO IGST=1,NGST
  DO IJ = KIJS,KIJL
    TAUNWX(IJ,IGST) = 0.0_JWRB
    TAUNWY(IJ,IGST) = 0.0_JWRB
    TAUWX(IJ,IGST)  = 0.0_JWRB
    TAUWY(IJ,IGST)  = 0.0_JWRB
    TAU(IJ,IGST) = 0.0_JWRB
  ENDDO
END DO  

!
!/    --- scale  friction velocity to wind speed (10m) in
!/        the boundary layer ----------------------------------------- /
!/    Donelan et al. (2006) used U10 or U_{λ/2} in their S_{in}
!/    parameterization. To avoid some disadvantages of using U10 or
!/    U_{λ/2}, Rogers et al. (2012) used the following engineering
!/    conversion:
!/                    UPROXY = SIN6WS * UST
!/
!/    SIN6WS = FRIC = 28.0 following Komen et al. (1984) (developed seas)
!/    SIN6WS        = 32.0 suggested by E. Rogers (2014) (young seas)
!
DO IGST=1,NGST
  SELECT CASE (IPHYS2_AIRSEA)
  CASE(0)
    DO IJ = KIJS,KIJL
        UPROXYGST(IJ,IGST) = FRIC * CDFAC * USTARGST(IJ,IGST) ! original, suggested by E. Rogers (2014) (young seas)
    END DO
  CASE(1)
    DO IJ = KIJS,KIJL
        UPROXYGST(IJ,IGST) = FRIC * CDFAC * USTARGST(IJ,IGST) ! following Komen et al. (1984) (developed seas) (FRIC=28)
    END DO
  CASE(2)
    DO IJ = KIJS,KIJL  
        UPROXYGST(IJ,IGST) = UABSGST(IJ,IGST) * CDFAC ! because FRIC=1/sqrt(CD), then this turns to purely a wind dependence (USTARGST cancels out)
    ENDDO
  END SELECT
END DO    
!
  ! To reshape from 1D to 2D: 
  !    K = RESHAPE( A          , (/ NANG, NFRE /))
  ! To reshape from 2D to 1D:
  !    A = RESHAPE( F(IJ,:,:)  , (/NSPEC/)    )
DO IJ = KIJS,KIJL
  A(IJ,:) = RESHAPE( FL1(IJ,:,:) , (/NSPEC/)) * CG2(IJ,:) / ( ZPI * SIG2 )! ACTION DENSITY SPECTRUM
!
!/ 1) --- calculate 1d action density spectrum (A(sigma)) and
!/        zero-out values less than 1.0E-32 to avoid NaNs when
!/        computing directional narrowness in step 4). --------------- /
  KK(IJ,:,:)      = RESHAPE(A(IJ,:),(/ NANG, NFRE /))
  
  ADENSIG(IJ,:) = SUM(KK(IJ,:,:),1) * SIG * DELTH ! Integrate over directions.
  
  KMAX(IJ,:) = MAXVAL(KK(IJ,:,:),1)
END DO
!
!/ 2) --- calculate normalised directional spectrum K(theta,sigma) --- /
DO M = 1,NFRE
  DO IJ = KIJS,KIJL
      IF (KMAX(IJ,M).LT.1.0E-34_JWRB) THEN
        KK(IJ,1:NANG,M) = 1.0_JWRB
      ELSE
        KK(IJ,1:NANG,M) = KK(IJ,1:NANG,M)/KMAX(IJ,M)
      END IF
  END DO
END DO
!
!/ 3) --- calculate normalised spectral saturation BN(M) ------------ /
DO IJ = KIJS,KIJL
  ANAR(IJ,:)    = 1.0_JWRB/( SUM(KK(IJ,:,:),1) * DELTH )          ! directional narrowness
  !
  !        SQRTBN  = SQRT( ANAR * ADENSIG * WN(IJ,:)**3 )
  SQRTBN(IJ,:)  = SQRT( ANAR(IJ,:) * ADENSIG(IJ,:) * WAVNUM(IJ,:)**3 )
END DO

DO K  = 1, NANG
  DO IJ = KIJS,KIJL
    SQRTBN2(IJ,IKN+(K-1)) = SQRTBN(IJ,:)          ! Calculate SQRTBN for
  END DO                                    ! the entire spectrum.
END DO                                    ! the entire spectrum.
!
!/ 4) --- calculate growth rate GAMMA and S for all directions for
!/        following winds (U10/c - 1 is positive; W1) and in 7) for
!/        adverse winds (U10/c -1 is negative, W2). W1 and W2
!/        complement one another. ------------------------------------ /
DO IGST=1,NGST
  DO IJ = KIJS,KIJL
    W1(IJ,:,IGST)= MAX(0.0_JWRB,                   &
  &                 UPROXYGST(IJ,IGST)*CINV2(IJ,:)*(ECOS2*COSU(IJ) + ESIN2*SINU(IJ)) - 1.0_JWRB)**2
!
    D(IJ,:,IGST) = (RAORW(IJ) ) * SIG2 * &
                (2.8_JWRB-(1.0_JWRB+TANH(10.0_JWRB*SQRTBN2(IJ,:)*W1(IJ,:,IGST)-11.0_JWRB)))*&
  &             SQRTBN2(IJ,:)*W1(IJ,:,IGST)
!
    IF (LLLOWWINDS .AND. UABSGST(IJ,IGST)<=1.5_JWRB) THEN
      ! Reduce growth rates for low winds (following Muhammad Yasrab's work)
      D(IJ,:,IGST) = D(IJ,:,IGST) - (4._JWRB*(RNU_WATER)*(WAVNUM(IJ,:)**2))
      S(IJ,:,IGST) = D(IJ,:,IGST) * A(IJ,:)
    ELSE
      ! Update spectrum as per normal
      S(IJ,:,IGST) = D(IJ,:,IGST) * A(IJ,:)
    END IF

  ENDDO
ENDDO
!
!/ 5) --- calculate reduction factor LFACT using non-directional
!         spectral density of the wind input ------------------------- /

DO IGST=1,NGST
  DO IJ = KIJS,KIJL
    SDENSIG(IJ,:,:,IGST) = RESHAPE(S(IJ,:,IGST)*SIG2/CG2(IJ,:),(/ NANG, NFRE /))

    CALL LFACTOR(SDENSIG(IJ,:,:,IGST), CINV1(IJ,:), UABSGST(IJ,IGST), USTARGST(IJ,IGST), WDWAVE(IJ),    &
&                      ROAIRN(IJ), SIG, DSII, LFACT(IJ,:,IGST), TAUWX(IJ,IGST), TAUWY(IJ,IGST), TAU(IJ,IGST))
  ENDDO
ENDDO

!
!/ 6) --- apply reduction (LFACT) to the entire spectrum ------------- /

LLFACT = .TRUE.
IF (LLFACT) THEN
  DO IGST=1,NGST
    DO IJ = KIJS,KIJL ! TODO: how to make more efficient? is difficult...
      IF (SUM(LFACT(IJ,:,IGST)) .LT. NFRE) THEN
        DO K = 1, NANG
            D(IJ,IKN+K-1,IGST) = D(IJ,IKN+K-1,IGST) * LFACT(IJ,:,IGST)
        END DO
        S(IJ,:,IGST) = D(IJ,:,IGST) * A(IJ,:)
      END IF
      DINPOS(IJ,:,:,IGST)  = RESHAPE(D(IJ,:,IGST),(/ NANG, NFRE /))
    ENDDO
  ENDDO
END IF

!
!/ 7) --- compute negative wind input for adverse winds. negative
!/        growth is typically smaller by a factor of ~2.5 (=.28/.11)
!/        than those for the favourable winds [Donelan, 2006, Eq. (7)].
!/        the factor is adjustable with NAMELIST parameter in
!/        ww3_grid.inp: '&SIN6 SINA0 = 0.04 /' ----------------------- /
DO IGST=1,NGST
  IF (ZSIN6A0.GT.0.0_JWRB) THEN
    DO IJ = KIJS,KIJL
      W2(IJ,:,IGST)  = MIN( 0.0_JWRB,UPROXYGST(IJ,IGST) * CINV2(IJ,:) * &
  &                               (ECOS2*COSU(IJ) + ESIN2*SINU(IJ)) - 1.0_JWRB )**2
      D(IJ,:,IGST)   = D(IJ,:,IGST) - ( RAORW(IJ) * SIG2 * ZSIN6A0 * &
                    (2.8_JWRB-(1.0_JWRB+TANH(10.0_JWRB*SQRTBN2(IJ,:)*W2(IJ,:,IGST) - 11.0_JWRB)))&
  &                 *SQRTBN2(IJ,:)*W2(IJ,:,IGST) )

      DINTOT(IJ,:,:,IGST)= RESHAPE(D(IJ,:,IGST),(/NANG,NFRE/))
      S(IJ,:,IGST)       = D(IJ,:,IGST) * A(IJ,:)

! !     --- compute negative component of the wave supported stresses
! !         from negative part of the wind input  ---------------------- /
      SDENSIG(IJ,:,:,IGST) = RESHAPE(S(IJ,:,IGST)*SIG2/CG2(IJ,:),(/ NANG, NFRE /))
      CALL TAU_WAVE_ATMOS(SDENSIG(IJ,:,:,IGST), CINV1(IJ,:), SIG, DSII, TAUNWX(IJ,IGST), TAUNWY(IJ,IGST) )
    ENDDO
  ELSE
    DO IJ = KIJS,KIJL
      DINTOT(IJ,:,:,IGST)=DINPOS(IJ,:,:,IGST)
    ENDDO
  END IF
ENDDO
!
DO IGST=1,NGST
  DO IJ = KIJS,KIJL
    TAUWGST(IJ,IGST)     = SQRT(TAUWX(IJ,IGST)**2+TAUWY(IJ,IGST)**2) / ROAIRN(IJ) ! KINEMATIC TAUW
    TAUWDIRGST(IJ,IGST)  = ATAN2(TAUWX(IJ,IGST),TAUWY(IJ,IGST))
    TAUNWGST(IJ,IGST)    = SQRT(TAUNWX(IJ,IGST)**2+TAUNWY(IJ,IGST)**2) / ROAIRN(IJ) ! KINEMATIC TAUNW
    ENDDO
  ENDDO

SELECT CASE (IPHYS2_AIRSEA)
  ! CASE(0)
  !   USTARGST(IJ,IGST)    = USTARGST(IJ,IGST) ! i.e. do nothing here, don't update USTAR because it is not true to WW3_ST6
  CASE(1,2)
  DO IGST=1,NGST
    DO IJ = KIJS,KIJL
        USTARGST(IJ,IGST)    = SQRT(TAU(IJ,IGST) / ROAIRN(IJ) ) 
    ENDDO
  ENDDO
END SELECT

! 8) --- Calculate SL, FL and SPOS needed for ecWAM ------------- /

DO IGST=1,NGST
  DO M = 1,NFRE
    DO K = 1, NANG
      DO IJ = KIJS,KIJL
        SLGST(IJ,K,M,IGST)   = DINTOT(IJ,K,M,IGST)*FL1(IJ,K,M)
        SPOSGST(IJ,K,M,IGST) = DINPOS(IJ,K,M,IGST)*FL1(IJ,K,M)
      END DO
    END DO
  END DO
END DO

DO IGST=1,NGST
  DO IJ = KIJS,KIJL
    FLGST(IJ,:,:,IGST)   = DINTOT(IJ,:,:,IGST)
  ENDDO
END DO

! 9) --- Averaging over gust components ------------- /
IGST=1
  DO IJ = KIJS,KIJL
    TAUWGST_AVG(IJ)     = TAUWGST(IJ,IGST)
    TAUWDIRGST_AVG(IJ)  = TAUWDIRGST(IJ,IGST)
    TAUNWGST_AVG(IJ)    = TAUNWGST(IJ,IGST)
    USTARGST_AVG(IJ)    = USTARGST(IJ,IGST)
    UPROXYGST_AVG(IJ)   = UPROXYGST(IJ,IGST)
    SLGST_AVG(IJ,:,:)   = SLGST(IJ,:,:,IGST)
    SPOSGST_AVG(IJ,:,:) = SPOSGST(IJ,:,:,IGST)
    FLGST_AVG(IJ,:,:)   = FLGST(IJ,:,:,IGST)
  END DO
DO IGST=2,NGST
  DO IJ = KIJS,KIJL
    TAUWGST_AVG(IJ)      = TAUWGST_AVG(IJ)     + TAUWGST(IJ,IGST)
    TAUWDIRGST_AVG(IJ)   = TAUWDIRGST_AVG(IJ)  + TAUWDIRGST(IJ,IGST)
    TAUNWGST_AVG(IJ)     = TAUNWGST_AVG(IJ)    + TAUNWGST(IJ,IGST)
    USTARGST_AVG(IJ)     = USTARGST_AVG(IJ)    + USTARGST(IJ,IGST)
    UPROXYGST_AVG(IJ)    = UPROXYGST_AVG(IJ)   + UPROXYGST(IJ,IGST)
    SLGST_AVG(IJ,:,:)    = SLGST_AVG(IJ,:,:)   + SLGST(IJ,:,:,IGST)
    SPOSGST_AVG(IJ,:,:)  = SPOSGST_AVG(IJ,:,:) + SPOSGST(IJ,:,:,IGST)
    FLGST_AVG(IJ,:,:)    = FLGST_AVG(IJ,:,:)   + FLGST(IJ,:,:,IGST)
  ENDDO
END DO

DO IJ = KIJS,KIJL
  TAUW(IJ)     = AVG_GST*TAUWGST_AVG(IJ)
  TAUWDIR(IJ)  = AVG_GST*TAUWDIRGST_AVG(IJ)
  TAUNW(IJ)    = AVG_GST*TAUNWGST_AVG(IJ)
  UFRIC(IJ)    = AVG_GST*USTARGST_AVG(IJ)
  UPROXY(IJ)   = AVG_GST*UPROXYGST_AVG(IJ)
  SL(IJ,:,:)   = AVG_GST*SLGST_AVG(IJ,:,:)
  SPOS(IJ,:,:) = AVG_GST*SPOSGST_AVG(IJ,:,:)
  FLD(IJ,:,:)  = AVG_GST*FLGST_AVG(IJ,:,:)

! 10) --- Calculate roughness length and charnock ------------- /

  USTM1     = 1.0_JWRB/MAX(UFRIC(IJ),EPSUS)                       ! Protect the code
  USTM2     = 1.0_JWRB/MAX(UFRIC(IJ)**2,EPSUS)                    ! Protect the code
  KUOUST    = MIN(50._JWRB,XKAPPA*WSWAVE(IJ)*USTM1)                 ! Protect the code 
  Z0        = ZNLEV / ( EXP(KUOUST) - 1.0_JWRB )
  Z0        = MAX(Z0, 0.0000001_JWRB)
  Z0M(IJ)   = Z0                                 ! Update z0 
  CHNKOG(IJ)    = ( Z0 - ZRN*USTM1 ) * USTM2     ! Update charnock (where Z0=Z0CH+Z0VIS from airsea_zbry)
  ALPHAOGMAXU10 = MIN(ALPHAMAX,AMAX+BMAX*WSWAVE(IJ))*GM1 ! protective code taken from outbeta (incl /G)
  CHNKOG(IJ)    = MIN(CHNKOG(IJ),ALPHAOGMAXU10)        ! protective code taken from outbeta (incl /G)

  IF(LLCAPCHNK) THEN
    CHARNOCK_MIN = CHNKMIN(WSWAVE(IJ))
    CHNKOG(IJ)   = MAX(CHARNOCK_MIN*GM1,CHNKOG(IJ))
  ELSE
    CHNKOG(IJ)   = MAX(CHNKOG(IJ), 1E-5_JWRB)
  ENDIF

  CHRNCK(IJ)    = CHNKOG(IJ)*G

! 11) --- PHIWA calculation using non-directional
!         spectral density of the wind input  ---------------------- /

  SPOSDENSIG = SPOS(IJ,:,:)
  SNEGDENSIG = SL(IJ,:,:) - SPOS(IJ,:,:)
  PHIWA(IJ)  = CALCPHIWA(SPOSDENSIG,SNEGDENSIG,DSII)
END DO

! XLLWS based on SL (mask for neg. input)
DO M = 1,NFRE
  DO K = 1, NANG
    DO IJ=KIJS,KIJL
      IF (SL(IJ,K,M)>0.0_JWRB) THEN
          XLLWS(IJ,K,M)=1.0_JWRB
      ELSE
          XLLWS(IJ,K,M)=0.0_JWRB
      END IF
    END DO
  END DO
END DO
! ---------------------

! input source term!!!! (end)
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
!$loki inline
CALL FEMEANWS(KIJS, KIJL, FL1, XLLWS, FMEANWS)

! COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!$loki inline
CALL FRCUTINDEX(KIJS, KIJL, FMEAN, FMEANWS, UFRIC, CICOVER, MIJ, RHOWGDFTH)

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SINFLX',1,ZHOOK_HANDLE)

END SUBROUTINE SINFLX_ZBRY
