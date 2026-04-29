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
 &                 FMEAN, HALP, FMEANWS,      &
 &                 FLM,                       &
 &                 UFRIC, TAUW, TAUWDIR,      & 
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

!     *CALL* *SINFLX_ZBRY (ICALL, NCALL, NGST, KIJS, KIJL, LUPDTUS,*
!    &                    FL1, WAVNUM, CGROUP, CINV,
!    &                    WSWAVE, WDWAVE, AIRD,
!    &                    RAORW, WSTAR, CICOVER,
!    &                    FMEAN, HALP, FMEANWS,
!    &                    FLM,
!    &                    UFRIC, TAUW, TAUWDIR,
!    &                    Z0M, Z0B, CHRNCK, PHIWA,
!    &                    FLD, SL, SPOS,
!    &                    MIJ, RHOWGDFTH, XLLWS)
!        *ICALL*     - CALL NUMBER.
!        *NCALL*     - TOTAL NUMBER OF CALLS.
!        *NGST*      - NUMBER OF GUST COMPONENTS.
!        *KIJS*      - INDEX OF FIRST GRIDPOINT.
!        *KIJL*      - INDEX OF LAST GRIDPOINT.
!        *LUPDTUS*   - IF TRUE, UFRIC/Z0M ARE UPDATED VIA AIRSEA.
!        *FL1*       - WAVE SPECTRUM.
!        *WAVNUM*    - WAVE NUMBER.
!        *CGROUP*    - GROUP VELOCITY.
!        *CINV*      - INVERSE PHASE VELOCITY.
!        *WSWAVE*    - WIND SPEED IN M/S.
!        *WDWAVE*    - WIND DIRECTION (OCEANOGRAPHIC CONVENTION).
!        *AIRD*      - AIR DENSITY (KG/M**3).
!        *RAORW*     - AIR/WATER DENSITY RATIO.
!        *WSTAR*     - FREE CONVECTION VELOCITY SCALE (M/S).
!        *CICOVER*   - SEA ICE COVER.
!        *FMEAN*     - MEAN FREQUENCY.
!        *HALP*      - 1/2 PHILLIPS PARAMETER.
!        *FMEANWS*   - MEAN FREQUENCY OF WINDSEA.
!        *FLM*       - SPECTRAL DENSITY MINIMUM VALUE.
!        *UFRIC*     - FRICTION VELOCITY (M/S).
!        *TAUW*      - WAVE STRESS ((M/S)**2).
!        *TAUWDIR*   - WAVE STRESS DIRECTION.
!        *Z0M*       - ROUGHNESS LENGTH (M).
!        *Z0B*       - BACKGROUND ROUGHNESS LENGTH.
!        *CHRNCK*    - CHARNOCK COEFFICIENT.
!        *PHIWA*     - ENERGY FLUX FROM WIND INTO WAVES.
!        *FLD*       - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!        *SL*        - TOTAL SOURCE FUNCTION ARRAY.
!        *SPOS*      - POSITIVE INPUT SOURCE COMPONENT.
!        *MIJ*       - LAST FREQUENCY INDEX OF PROGNOSTIC RANGE.
!        *RHOWGDFTH* - WATER DENSITY * G * DF * DTHETA.
!        *XLLWS*     - WINDSEA MASK FROM INPUT SOURCE TERM.

!     METHOD.
!     -------

!       SEE REFERENCE.

!     EXTERNALS.
!     ----------
!     TAU_WAVE_ATMOS
!     LFACTOR

!     ORIGIN.
!     ----------
!     Adapted from Babanin Young Donelan & Banner (ZBRY) physics 
!     as implemented as ST6 in WAVEWATCH-III
!     Implementation into ECWAM DECEMBER 2021 by J. Kousal


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU    ,LLCAPCHNK

      USE YOWWNDG  , ONLY : ICODE    ,ICODE_CPL

      USE YOWFRED  , ONLY : FR       ,TH       ,DFIM     ,COSTH  ,SINTH, ZPIFR, DELTH, FRATIO, FRIC, &
     &                      SIG      ,DSII     ,SIGM1    ,DF     ,NFRE_EXT
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,GM1      ,EPSMIN, EPSUS, ZPI, ROWATER
      USE YOWPHYS  , ONLY : ZALP     ,TAUWSHELTER, XKAPPA, BETAMAXOXKAPPA2,    &
     &                      RNU      ,RNUM, &
     &                      SWELLF   ,SWELLF2  ,SWELLF3  ,SWELLF4  , SWELLF5, &
     &                      SWELLF6  ,SWELLF7  ,SWELLF7M1, Z0RAT   ,Z0TUBMAX , &
     &                      ABMIN  ,ABMAX, DTHRN_A  ,DTHRN_U, RNU_WATER, &
     &                      ZSIN6A0, FRQMAX, LLFACT
      USE YOWTEST  , ONLY : IU06
      USE YOWTABL  , ONLY : IAB      ,SWELLFT
      USE YOWSTAT  , ONLY : IPHYS2_AIRSEA, LLLOWWINDS, ZCDFAC

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

REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: FMEAN  !! MEAN FREQUENCY.
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: HALP   !! 1/2 PHILLIPS PARAMETER  
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(OUT) :: FMEANWS  !! MEAN FREQUENCY OF THE WINDSEA.
REAL(KIND=JWRB), DIMENSION(KIJL,NANG), INTENT(IN) :: FLM  !! SPECTAL DENSITY MINIMUM VALUE
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: UFRIC !! FRICTION VELOCITY IN M/S.
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

INTEGER(KIND=JWIM) :: IJ, K, M, IND, IGST

REAL(KIND=JWRB), DIMENSION(NANG,NFRE) :: SPOSDENSIG, SNEGDENSIG

REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE) :: A
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE)      :: ADENSIG, KMAX, ANAR, SQRTBN
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE) :: KK
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE,NGST) :: W1, W2, S, D
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE,NGST)      :: LFACT
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE,NGST) :: DINPOS, DINTOT


REAL(KIND=JWRB), DIMENSION(KIJL,NGST) :: TAUWX, TAUWY ! Component of the wave-supported stress
REAL(KIND=JWRB), DIMENSION(KIJL) :: COSU, SINU
REAL(KIND=JWRB), DIMENSION(KIJL,NGST) :: UPROXYGST

REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: XK, CGG_WAM, CM

! For USTAR, Z0, CHNK
REAL(KIND=JWRB), DIMENSION(KIJL,NGST) :: TAU
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
REAL(KIND=JWRB) :: SUMDIR

REAL(KIND=JWRB), DIMENSION(KIJL)      :: ROAIRN, CHNKOG

! For GUSTINESS
REAL(KIND=JWRB)                          :: AVG_GST
REAL(KIND=JWRB), DIMENSION(KIJL)      :: SIG_N, SIG_U10, TAUWGST_AVG, TAUWDIRGST_AVG, USTARGST_AVG
REAL(KIND=JWRB), DIMENSION(KIJL,NGST) :: TAUWGST, TAUWDIRGST, UABSGST, USTARGST
! REAL(KIND=JWRB), DIMENSION(KIJL,NGST) :: TAUNWGST
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE) :: SLGST_AVG, SPOSGST_AVG, FLGST_AVG
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE,NGST) :: SLGST, SPOSGST, FLGST
REAL(KIND=JWRB), DIMENSION(NANG,NFRE) :: SPOS_IJ, SNEG_IJ
REAL(KIND=JWRB), DIMENSION(NANG,NFRE) :: SDENSIG_IJ
LOGICAL, DIMENSION(KIJL) :: LREDUCE

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

IF(LUPDTUS) THEN
  !$loki inline
  CALL AIRSEA (KIJS, KIJL,                                  &
&              U10=WSWAVE, U10DIR=WDWAVE, TAUW=TAUW, TAUWDIR=TAUWDIR, &
&              US=UFRIC, Z0=Z0M, Z0B=Z0B, CHRNCK=CHRNCK, ICODE_WND=ICODE_WND, IUSFG=IUSFG) 

ENDIF


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! input source term!!!! (start)

!     Wind height
ZNLEV    = 10._JWRB

DO M=1,NFRE
  DO IJ=KIJS,KIJL
    CM(IJ,M)     = WAVNUM(IJ,M)*SIGM1(M)
  ENDDO
ENDDO


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

!/  --- Main loop over LOC ----------------------------------- /

!/ 0) --- set up a basic variables ----------------------------------- /
DO IJ = KIJS,KIJL
  COSU(IJ)  = COS(WDWAVE(IJ))
  SINU(IJ)  = SIN(WDWAVE(IJ))
END DO  
!
DO IGST=1,NGST
  DO IJ = KIJS,KIJL
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
  CASE(0,1,2)
    ! IPHYS2_AIRSEA=0,1,2 use USTAR
    DO IJ = KIJS,KIJL
        UPROXYGST(IJ,IGST) = 32.0_JWRB * ZCDFAC * USTARGST(IJ,IGST)
    END DO
  CASE(3)
    ! IPHYS2_AIRSEA=3 is based on wind directly
    DO IJ = KIJS,KIJL  
        UPROXYGST(IJ,IGST) = UABSGST(IJ,IGST) * ZCDFAC ! because FRIC=1/sqrt(CD), then this turns to purely a wind dependence (USTARGST cancels out)
    ENDDO
  END SELECT
END DO    

! ------------------------------------------------------------------- /
! START: using intrinsic frequency spectra

!
!
!/ 1) --- calculate 1d action density spectrum (A(sigma)) and
!/        zero-out values less than 1.0E-32 to avoid NaNs when
!/        computing directional narrowness in step 4). --------------- /

DO M = 1, NFRE
  DO K = 1, NANG
    DO IJ = KIJS,KIJL
      A(IJ,K,M) = FL1(IJ,K,M) * CGROUP(IJ,M) / ( ZPI * SIG(M) )  ! ACTION DENSITY SPECTRUM
      KK(IJ,K,M) = A(IJ,K,M)
    END DO
  END DO

  DO IJ = KIJS,KIJL
    ADENSIG(IJ,M) = 0.0_JWRB
    KMAX(IJ,M) = 0.0_JWRB
  END DO
  DO K = 1, NANG
    DO IJ = KIJS,KIJL
      ADENSIG(IJ,M) = ADENSIG(IJ,M) + KK(IJ,K,M)
      KMAX(IJ,M) = MAX(KMAX(IJ,M), KK(IJ,K,M))
    END DO
  END DO
  DO IJ = KIJS,KIJL
    ADENSIG(IJ,M) = ADENSIG(IJ,M) * SIG(M) * DELTH  ! Integrate over directions.
  END DO
!
!/ 2) --- calculate normalised directional spectrum K(theta,sigma) --- /

  DO K = 1, NANG
    DO IJ = KIJS,KIJL
      IF (KMAX(IJ,M).LT.1.0E-34_JWRB) THEN
        KK(IJ,K,M) = 1.0_JWRB
      ELSE
        KK(IJ,K,M) = KK(IJ,K,M)/KMAX(IJ,M)
      END IF
    END DO
  END DO
!
!/ 3) --- calculate normalised spectral saturation BN(M) ------------ /
  DO IJ = KIJS,KIJL
    ANAR(IJ,M) = 0.0_JWRB
  END DO
  DO K = 1, NANG
    DO IJ = KIJS,KIJL
      ANAR(IJ,M) = ANAR(IJ,M) + KK(IJ,K,M)
    END DO
  END DO
  DO IJ = KIJS,KIJL
    ANAR(IJ,M) = 1.0_JWRB/( ANAR(IJ,M) * DELTH )          ! directional narrowness
    SQRTBN(IJ,M)  = SQRT( ANAR(IJ,M) * ADENSIG(IJ,M) * WAVNUM(IJ,M)**3 )
  END DO

END DO
!
!/ 4) --- calculate growth rate GAMMA and S for all directions for
!/        following winds (U10/c - 1 is positive; W1) and in 7) for
!/        adverse winds (U10/c -1 is negative, W2). W1 and W2
!/        complement one another. ------------------------------------ /
DO IGST=1,NGST
  DO M = 1, NFRE
    DO K = 1, NANG
      DO IJ = KIJS,KIJL
        W1(IJ,K,M,IGST) = MAX(0.0_JWRB,                   &
  &                 UPROXYGST(IJ,IGST)*CINV(IJ,M)*(COSTH(K)*COSU(IJ) + SINTH(K)*SINU(IJ)) - 1.0_JWRB)**2
        D(IJ,K,M,IGST) = (RAORW(IJ)) * SIG(M) * &
                (2.8_JWRB-(1.0_JWRB+TANH(10.0_JWRB*SQRTBN(IJ,M)*W1(IJ,K,M,IGST)-11.0_JWRB)))* &
  &             SQRTBN(IJ,M)*W1(IJ,K,M,IGST)

        IF (LLLOWWINDS .AND. UABSGST(IJ,IGST)<=1.5_JWRB) THEN
          ! Reduce growth rates for low winds (following Muhammad Yasrab's work)
          D(IJ,K,M,IGST) = D(IJ,K,M,IGST) - (4._JWRB*(RNU_WATER)*(WAVNUM(IJ,M)**2))
        END IF
        S(IJ,K,M,IGST) = D(IJ,K,M,IGST) * A(IJ,K,M)
      END DO
    END DO
  END DO
END DO

IF (LLFACT) THEN ! TODO: how to make more efficient? is difficult...

  !/ 5) --- calculate reduction factor LFACT using non-directional
  !         spectral density of the wind input ------------------------- /

  DO IGST=1,NGST
    DO IJ = KIJS,KIJL
      DO M = 1, NFRE
        DO K = 1, NANG
          SDENSIG_IJ(K,M) = S(IJ,K,M,IGST)*SIG(M)/CGROUP(IJ,M)
        END DO
      END DO
      CALL LFACTOR(SDENSIG_IJ, CINV(IJ,:), UABSGST(IJ,IGST), USTARGST(IJ,IGST), UPROXYGST(IJ,IGST), WDWAVE(IJ),    &
  &                      ROAIRN(IJ), LFACT(IJ,:,IGST), LREDUCE(IJ), TAUWX(IJ,IGST), TAUWY(IJ,IGST), TAU(IJ,IGST))
    ENDDO

  !/ 6) --- apply reduction (LFACT) to the entire spectrum ------------- /

    DO M = 1, NFRE
      DO K = 1, NANG
        DO IJ = KIJS,KIJL
          IF (LREDUCE(IJ)) THEN
            D(IJ,K,M,IGST) = D(IJ,K,M,IGST) * LFACT(IJ,M,IGST)
            S(IJ,K,M,IGST) = D(IJ,K,M,IGST) * A(IJ,K,M)
          END IF
          DINPOS(IJ,K,M,IGST) = D(IJ,K,M,IGST)
        END DO
      END DO
    END DO

  ENDDO
END IF

!
!/ 7) --- compute negative wind input for adverse winds. negative
!/        growth is typically smaller by a factor of ~2.5 (=.28/.11)
!/        than those for the favourable winds [Donelan, 2006, Eq. (7)].
!/        the factor is adjustable with namelist parameter ZSIN6A0 ---- /
DO IGST=1,NGST
  IF (ZSIN6A0.GT.0.0_JWRB) THEN
    DO M = 1, NFRE
      DO K = 1, NANG
        DO IJ = KIJS,KIJL
          W2(IJ,K,M,IGST)  = MIN( 0.0_JWRB,UPROXYGST(IJ,IGST) * CINV(IJ,M) * &
  &                               (COSTH(K)*COSU(IJ) + SINTH(K)*SINU(IJ)) - 1.0_JWRB )**2
          D(IJ,K,M,IGST)   = D(IJ,K,M,IGST) - ( RAORW(IJ) * SIG(M) * ZSIN6A0 * &
                    (2.8_JWRB-(1.0_JWRB+TANH(10.0_JWRB*SQRTBN(IJ,M)*W2(IJ,K,M,IGST) - 11.0_JWRB)))&
  &                 *SQRTBN(IJ,M)*W2(IJ,K,M,IGST) )
          DINTOT(IJ,K,M,IGST) = D(IJ,K,M,IGST)
          S(IJ,K,M,IGST)      = D(IJ,K,M,IGST) * A(IJ,K,M)

! !     --- compute negative component of the wave supported stresses
! !         from negative part of the wind input  ---------------------- /
        END DO
      END DO
    END DO
  ELSE
    DO M = 1, NFRE
      DO K = 1, NANG
        DO IJ = KIJS,KIJL
          DINTOT(IJ,K,M,IGST) = DINPOS(IJ,K,M,IGST)
        END DO
      END DO
    END DO
  END IF
ENDDO
!
DO IGST=1,NGST
  DO IJ = KIJS,KIJL
    TAUWGST(IJ,IGST)     = SQRT(TAUWX(IJ,IGST)**2+TAUWY(IJ,IGST)**2) / ROAIRN(IJ) ! KINEMATIC TAUW
    TAUWDIRGST(IJ,IGST)  = ATAN2(TAUWX(IJ,IGST),TAUWY(IJ,IGST))
    ENDDO
  ENDDO

SELECT CASE (IPHYS2_AIRSEA)
  ! IPHYS2_AIRSEA=0      DO NOT USE STRESS BALANCE (LFAC) TO UPDATE USTAR
  ! IPHYS2_AIRSEA=1,2,3         USE STRESS BALANCE (LFAC) TO UPDATE USTAR
  CASE(1,2,3)
  DO IGST=1,NGST
    DO IJ = KIJS,KIJL
        USTARGST(IJ,IGST)    = SQRT(TAU(IJ,IGST) / ROAIRN(IJ) ) 
    ENDDO
  ENDDO
END SELECT

! END: using intrinsic frequency spectra
! ------------------------------------------------------------------- /
! NOTE: below this line is using frequency (i.e. not intrinsic frequency)

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
  DO M = 1,NFRE
    DO K = 1, NANG
      DO IJ = KIJS,KIJL
        FLGST(IJ,K,M,IGST) = DINTOT(IJ,K,M,IGST)
      END DO
    END DO
  END DO
END DO

! 9) --- Averaging over gust components ------------- /
DO IJ = KIJS,KIJL
  TAUWGST_AVG(IJ)     = 0.0_JWRB
  TAUWDIRGST_AVG(IJ)  = 0.0_JWRB
  USTARGST_AVG(IJ)    = 0.0_JWRB
END DO
DO M = 1,NFRE
  DO K = 1, NANG
    DO IJ = KIJS,KIJL
      SLGST_AVG(IJ,K,M)   = 0.0_JWRB
      SPOSGST_AVG(IJ,K,M) = 0.0_JWRB
      FLGST_AVG(IJ,K,M)   = 0.0_JWRB
    END DO
  END DO
END DO

DO IGST=1,NGST
  DO IJ = KIJS,KIJL
    TAUWGST_AVG(IJ)      = TAUWGST_AVG(IJ)     + TAUWGST(IJ,IGST)
    TAUWDIRGST_AVG(IJ)   = TAUWDIRGST_AVG(IJ)  + TAUWDIRGST(IJ,IGST)
    USTARGST_AVG(IJ)     = USTARGST_AVG(IJ)    + USTARGST(IJ,IGST)
  ENDDO
  DO M = 1,NFRE
    DO K = 1, NANG
      DO IJ = KIJS,KIJL
        SLGST_AVG(IJ,K,M)    = SLGST_AVG(IJ,K,M)   + SLGST(IJ,K,M,IGST)
        SPOSGST_AVG(IJ,K,M)  = SPOSGST_AVG(IJ,K,M) + SPOSGST(IJ,K,M,IGST)
        FLGST_AVG(IJ,K,M)    = FLGST_AVG(IJ,K,M)   + FLGST(IJ,K,M,IGST)
      END DO
    END DO
  END DO
END DO

DO IJ = KIJS,KIJL
  TAUW(IJ)     = AVG_GST*TAUWGST_AVG(IJ)
  TAUWDIR(IJ)  = AVG_GST*TAUWDIRGST_AVG(IJ)
  UFRIC(IJ)    = AVG_GST*USTARGST_AVG(IJ)
END DO

DO M = 1,NFRE
  DO K = 1, NANG
    DO IJ = KIJS,KIJL
      SL(IJ,K,M)   = AVG_GST*SLGST_AVG(IJ,K,M)
      SPOS(IJ,K,M) = AVG_GST*SPOSGST_AVG(IJ,K,M)
      FLD(IJ,K,M)  = AVG_GST*FLGST_AVG(IJ,K,M)
    END DO
  END DO
END DO

DO IJ = KIJS,KIJL

! 10) --- Calculate roughness length and charnock ------------- /

  USTM1     = 1.0_JWRB/MAX(UFRIC(IJ),EPSUS)                       ! Protect the code
  USTM2     = 1.0_JWRB/MAX(UFRIC(IJ)**2,EPSUS)                    ! Protect the code
  KUOUST    = MIN(50._JWRB,XKAPPA*WSWAVE(IJ)*USTM1)                 ! Protect the code 
  Z0        = ZNLEV / ( EXP(KUOUST) - 1.0_JWRB )
  Z0        = MAX(Z0, 0.0000001_JWRB)
  Z0M(IJ)   = Z0                                 ! Update z0 
  CHNKOG(IJ)    = ( Z0 - RNU_WATER*USTM1 ) * USTM2     ! Update charnock (where Z0=Z0CH+Z0VIS from airsea_zbry)
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
!              CALCPHIWA(SPOS        ,SNEG                     )               
  DO M = 1,NFRE
    DO K = 1, NANG
      SPOS_IJ(K,M) = SPOS(IJ,K,M)
      SNEG_IJ(K,M) = SL(IJ,K,M) - SPOS(IJ,K,M)
    END DO
  END DO
  PHIWA(IJ)  = CALCPHIWA(SPOS_IJ,SNEG_IJ)
END DO

! XLLWS based on SL (mask for pos. input)
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
