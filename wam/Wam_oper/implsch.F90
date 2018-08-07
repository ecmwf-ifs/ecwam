      SUBROUTINE IMPLSCH (FL3, IJS, IJL, IG,                            &
     &                    THWOLD, USOLD,                                &
     &                    TAUW, Z0OLD,                                  &
     &                    ROAIRO, WSTAROLD,                             &
     &                    CICVR, CIWA,                                  &
     &                    U10NEW, THWNEW, USNEW,                        &
     &                    Z0NEW, ROAIRN, WSTARNEW,                      &
     &                    MIJ, XLLWS)

! ----------------------------------------------------------------------

!**** *IMPLSCH* - IMPLICIT SCHEME FOR TIME INTEGRATION OF SOURCE
!****             FUNCTIONS.

!     S.D.HASSELMANN.  MPI
!     H. GUENTHER AND L. ZAMBRESKY  OPTIMIZATION PERFORMED.
!     H. GUENTHER      GKSS/ECMWF   OCTOBER 1989  NEW WIND FIELD
!                                                 INTERFACE AND
!                                                 TIME COUNTING
!     P.A.E.M. JANSSEN KNMI         AUGUST  1990  COUPLED MODEL
!     H. GUENTHER      GKSS/ECMWF   JUNE    1991  NEW SEPARATION OF
!                                                  DIAG- AND PROGNOSTIC
!                                                  PART OF SPECTRUM.
!     P.A.E.M. JANSSEN ECMWF        FEBRUARY 1995  ADD MINIMUM VALUE
!                                                  (FLMIN).
!     J. BIDLOT ECMWF               FEBRUARY 1996 MESSAGE PASSING
!     J. BIDLOT ECMWF               FEBRUARY 1997 MESSAGE PASSING
!     J. BIDLOT ECMWF               FEBRUARY 2001 MODIFY CALLING ORDER
!                                   BETWEEN SINPUT-STRESSO-AIRSEA
!     S. ABDALLA       ECMWF        OCTOBER 2001
!                                   INCLUSION OF AIR DENSITY AND Zi/L.

!*    PURPOSE.
!     --------

!       THE IMPLICIT SCHEME ENABLES THE USE OF A TIMESTEP WHICH IS
!       LARGE COMPARED WITH THE CHARACTERISTIC DYNAMIC TIME SCALE.
!       THE SCHEME IS REQUIRED FOR THE HIGH FREQUENCIES WHICH
!       RAPIDLY ADJUST TO A QUASI-EQUILIBRIUM.

!**   INTERFACE.
!     ----------

!       *CALL* *IMPLSCH (FL3, FL, IJS, IJL, IG,
!    1                    THWOLD,USOLD,TAUW,Z0OLD,
!    &                    ROAIRO, WSTAROLD, 
!    &                    CICVR, CIWA,
!    2                    U10NEW,THWNEW,USNEW,Z0NEW,ROAIRN,WSTARNEW,
!    &                    MIJ,  XLLWS)
!          *FL3*    - FREQUENCY SPECTRUM(INPUT AND OUTPUT).
!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT
!          *IG*     - BLOCK NUMBER
!      *U10NEW*    NEW WIND SPEED IN M/S.
!      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!      *THWOLD*    INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
!                  WIND VELOCITY.
!      *USNEW*     NEW FRICTION VELOCITY IN M/S.
!      *USOLD*     INTERMEDIATE STORAGE OF MODULUS OF FRICTION
!                  VELOCITY.
!      *Z0NEW*     ROUGHNESS LENGTH IN M.
!      *Z0OLD*     INTERMEDIATE STORAGE OF ROUGHNESS LENGTH IN
!                  M.
!      *TAUW*      WAVE STRESS IN (M/S)**2
!      *ROAIRN*    AIR DENSITY IN KG/M3.
!      *ROAIRO*    INTERMEDIATE STORAGE OF AIR DENSITY.
!      *WSTARNEW*  FREE CONVECTION VELOCITY SCALE (M/S)
!      *WSTAROLD*   INTERMEDIATE STORAGE OF WSTAR
!      *CICVR*     SEA ICE COVER.
!      *CIWA*      SEA ICE WAVE ATTENUATION.
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
!       BETWEEN TWO PROPAGATION CALLS - ARRAY FL3 CONTAINS THE
!       SPECTRUM AND FL IS USED AS AN INTERMEDIATE STORAGE FOR THE
!       DIAGONAL TERM OF THE FUNCTIONAL MATRIX.

!     EXTERNALS.
!     ---------

!       *AIRSEA*    - SURFACE LAYER STRESS AND ROUGHNESS LENGTH.
!       *FEMEAN*    - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT.
!       *INCDATE*   - UPDATE DATE TIME GROUP.
!SHALLOW
!       *SBOTTOM*   - COMPUTES BOTTOM DISSIPATION SOURCE TERM AND
!                     LINEAR CONTRIBUTION TO FUNCTIONAL MATRIX.
!SHALLOW
!       *SDISSIP*   - COMPUTATION OF DISSIPATION SOURCE FUNCTION
!                     AND LINEAR CONTRIBUTION OF DISSIPATION TO
!                     FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
!       *SEMEAN*    - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.
!       *SINPUT*    - COMPUTATION OF INPUT SOURCE FUNCTION, AND
!                     LINEAR CONTRIBUTION OF INPUT SOURCE FUNCTION
!                     TO FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
!       *SNONLIN*   - COMPUTATION OF NONLINEAR TRANSFER RATE AND
!                     DIAGONAL LINEAR CONTRIBUTION OF NONLINEAR SOURCE
!                     FUNCTION TO  FUNCTIONAL MATRIX.
!       *STRESSO*   - COMPUTATION NORMALISED WAVE STRESS.
!           !!!!!!! MAKE SURE THAT SINPUT IS CALLED FIRST, STRESSO
!           !!!!!!! NEXT, AND THEN THE REST OF THE SOURCE FUNCTIONS.
!       *FRCUTINDEX*
!       *IMPHFTAIL*

!     REFERENCE.
!     ----------

!       S. HASSELMANN AND K. HASSELMANN, "A GLOBAL WAVE MODEL",
!       30/6/85 (UNPUBLISHED NOTE)

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU    ,LWFLUX   , LWVFLX_SNL , LWNEMOCOU
      USE YOWCOUT  , ONLY : LWFLUXOUT 
      USE YOWFRED  , ONLY : FR       ,TH       ,DELTH       ,FRM5     , &
     &            COFRM4   ,FLMAX
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWICE   , ONLY : FLMIN    ,LCIWABR
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSHAL  , ONLY : DEPTH    ,INDEP    ,                        &
     &            IODP     ,IOBND    ,CINV     ,EMAXDPT
      USE YOWSTAT  , ONLY : IDELT    ,ISHALLO  ,CDTPRO   ,LBIWBK
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWUNPOOL ,ONLY : LLUNSTR
      USE YOWWNDG  , ONLY : ICODE    ,ICODE_CPL

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "imphftail.intfb.h"
#include "sdepthlim.intfb.h"
#include "airsea.intfb.h"
#include "ciwabr.intfb.h"
#include "femeanws.intfb.h"
#include "fkmean.intfb.h"
#include "frcutindex.intfb.h"
#include "sbottom.intfb.h"
#include "sdissip.intfb.h"
#include "sdiwbk.intfb.h"
#include "sinput.intfb.h"
#include "snonlin.intfb.h"
#include "stresso.intfb.h"
#include "wnfluxes.intfb.h"
#include "wrong_wnfluxes.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, IG
      INTEGER(KIND=JWIM), INTENT(OUT) :: MIJ(IJS:IJL)

      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(INOUT) :: THWOLD, USOLD, Z0OLD
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(INOUT) :: TAUW, ROAIRO, WSTAROLD
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: CICVR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(INOUT) :: U10NEW, USNEW 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: THWNEW
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: Z0NEW
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: ROAIRN, WSTARNEW
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE), INTENT(IN) :: CIWA
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL3
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(OUT) :: XLLWS

      INTEGER(KIND=JWIM) :: IJ, K, M
      INTEGER(KIND=JWIM) :: ILEV
      INTEGER(KIND=JWIM) :: ICODE_WND, ICODE_WAM

      REAL(KIND=JWRB) :: DELT, XIMP, DELT5
      REAL(KIND=JWRB) :: GTEMP1, GTEMP2, FLHAB
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: GAMOF
      REAL(KIND=JWRB) :: DELFL(NFRE)
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: EMEANALL, FMEANALL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: EMEANWS, FMEANWS, USFM, GADIAG 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: F1MEAN, AKMEAN, XKMEAN 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: PHIWA
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: DPTHREDUC
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG) :: FLM 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE) :: TEMP
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE) :: RHOWGDFTH
!     *FL*  DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!     *SL*  TOTAL SOURCE FUNCTION ARRAY.
!     *SPOS* : POSITIVE SINPUT ONLY
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: FL, SL, SPOS
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: CIREDUC 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: SSOURCE 

      LOGICAL :: LCFLX

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('IMPLSCH',0,ZHOOK_HANDLE)

!*    1. INITIALISATION.
!        ---------------

      DELT = IDELT
      XIMP = 1.0_JWRB
      DELT5 = XIMP*DELT

      LCFLX=LWFLUX.OR.LWFLUXOUT.OR.LWNEMOCOU

      IF (LWCOU) THEN
        ICODE_WND = ICODE_CPL
      ELSE
        ICODE_WND = ICODE
      ENDIF
      ICODE_WAM=3

! ----------------------------------------------------------------------

!*    2. COMPUTATION OF IMPLICIT INTEGRATION.
!        ------------------------------------

!         INTEGRATION IS DONE FROM CDATE UNTIL CDTPRO FOR A BLOCK
!         OF LATITUDES BETWEEN PROPAGATION CALLS.


!     REDUCE WAVE ENERGY IF LARGER THAN DEPTH LIMITED WAVE HEIGHT
      IF(ISHALLO.NE.1 .AND. LBIWBK) THEN
         CALL SDEPTHLIM(IJS,IJL,EMAXDPT(IJS),FL3)
      ENDIF

!*    2.2 COMPUTE MEAN PARAMETERS.
!        ------------------------


      CALL FKMEAN(FL3, IJS, IJL, EMEANALL, FMEANALL, F1MEAN, AKMEAN, XKMEAN)

      DO K=1,NANG
        DO IJ=IJS,IJL
          FLM(IJ,K)=FLMIN*MAX(0.0_JWRB,COS(TH(K)-THWNEW(IJ)))**2
        ENDDO
      ENDDO

!     COMPUTE DAMPING COEFFICIENT DUE TO FRICTION ON BOTTOM OF THE SEA ICE.
!!! testing sea ice attenuation (might need to restrict usage when needed)
      IF(LCIWABR) THEN
        CALL CIWABR(IJS, IJL, CICVR, FL3, CIREDUC) 
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=IJS,IJL
              CIREDUC(IJ,K,M)=CIWA(IJ,M)*CIREDUC(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=IJS,IJL
              CIREDUC(IJ,K,M)=CIWA(IJ,M)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    2.3 COMPUTATION OF SOURCE FUNCTIONS.
!         --------------------------------

!*      2.3.1 INITIALISE SOURCE FUNCTION AND DERIVATIVE ARRAY.
!             ------------------------------------------------

!!    FL AND SL ARE INITIALISED IN SINPUT

      ILEV=1
      CALL AIRSEA (U10NEW, TAUW, USNEW, Z0NEW, IJS, IJL, ILEV,ICODE_WND)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED BEFORE DO LOOP'
        CALL FLUSH (IU06)
      ENDIF

!*    2.3.2 ADD SOURCE FUNCTIONS AND WAVE STRESS.
!           -------------------------------------

      CALL SINPUT (FL3, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,             &
     &             ROAIRN, WSTARNEW, SL, SPOS, XLLWS)

!     diagnostic (not coded efficiently)
!!      DO IJ=IJS,IJL
!!        DO M=1,NFRE
!!          GAMOF=0.0_JWRB
!!          DO K=1,NANG
!!            GAMOF=GAMOF+FL(IJ,K,M)*DELTH
!!          ENDDO
!!          write(*,*) 'debile gamma ',CDTPRO," ",                     &
!!     &        USNEW(IJ)*CINV(INDEP(IJ),M),GAMOF/FR(M)
!!        ENDDO
!!      ENDDO

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: 1st SINPUT CALLED'
        CALL FLUSH (IU06)
      ENDIF

!     MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
      CALL FEMEANWS(FL3, IJS, IJL, EMEANWS, FMEANWS, XLLWS)

!     COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
      CALL FRCUTINDEX(IJS, IJL, FMEANALL, FMEANWS, USNEW, CICVR,        &
     &                MIJ, RHOWGDFTH)

      CALL STRESSO (FL3, SPOS, IJS, IJL,                                &
     &              MIJ, RHOWGDFTH,                                     &
     &              THWNEW, USNEW, Z0NEW, ROAIRN,                       &
     &              TAUW, PHIWA)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: STRESSO 1st CALLED'
        CALL FLUSH (IU06)
      ENDIF

      CALL AIRSEA (U10NEW, TAUW, USNEW, Z0NEW, IJS, IJL, ILEV,ICODE_WAM)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED'
        CALL FLUSH (IU06)
      ENDIF

!     RE-IMPOSE HIGH FREQUENCY TAIL
      CALL IMPHFTAIL (IJS, IJL, MIJ, FLM, FL3)


!*    REEVALUATE WIND INPUT SOURCE TERM
!     ---------------------------------

      CALL SINPUT (FL3, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,             &
     &             ROAIRN, WSTARNEW, SL, SPOS, XLLWS)

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: 2nd SINPUT CALL'
        CALL FLUSH (IU06)
      ENDIF

      IF(LCFLX) THEN
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=IJS,IJL
              SSOURCE(IJ,K,M) = SL(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      CALL STRESSO (FL3, SPOS, IJS, IJL,                                &
     &              MIJ, RHOWGDFTH,                                     &
     &              THWNEW, USNEW, Z0NEW, ROAIRN,                       &
     &              TAUW, PHIWA)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: STRESSO 2nd CALL'
        CALL FLUSH (IU06)
      ENDIF


!     2.3.3 ADD THE OTHER SOURCE TERMS.
!           ---------------------------

      CALL SDISSIP (FL3 ,FL, SL, IJS, IJL,                              &
     &              EMEANALL, F1MEAN, XKMEAN,                           &
     &              USNEW, THWNEW, ROAIRN)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SDISSIP CALLED'
        CALL FLUSH (IU06)
      ENDIF

      IF(LCFLX .AND. .NOT.LWVFLX_SNL) THEN

!        CALL WRONG_WNFLUXES (IJS, IJL,                                  &
!     &                       MIJ, RHOWGDFTH,                            &
!     &                       SSOURCE, SL,                               &
!     &                       PHIWA,                                     &
!     &                       EMEANALL, F1MEAN, U10NEW, THWNEW,          &
!     &                       USNEW, ROAIRN, .TRUE.)
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=IJS,IJL
!!!1          SSOURCE should only contain the positive contribution from sinput
              SSOURCE(IJ,K,M) = SL(IJ,K,M)-SSOURCE(IJ,K,M)+SPOS(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO
        CALL WNFLUXES (IJS, IJL,                                        &
     &                 MIJ, RHOWGDFTH,                                  &
     &                 SSOURCE, CICVR,                                  &
     &                 PHIWA,                                           &
     &                 EMEANALL, F1MEAN, U10NEW, THWNEW,                &
     &                 USNEW, ROAIRN, .TRUE.)
      ENDIF


      CALL SNONLIN (FL3, FL, IJS, IJL, IG, SL, AKMEAN)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SNONLIN CALLED'
        CALL FLUSH (IU06)
      ENDIF

      IF(LCFLX .AND. LWVFLX_SNL) THEN
!!!!!!  SL must only contain contributions contributed to fluxes into the oceans
!       MODULATE SL BY IMPLICIT FACTOR
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=IJS,IJL
!!!1          SSOURCE should only contain the positive contribution from sinput
              SSOURCE(IJ,K,M) = SL(IJ,K,M)-SSOURCE(IJ,K,M)+SPOS(IJ,K,M)
              GTEMP1 = MAX((1.0_JWRB-DELT5*FL(IJ,K,M)),1.0_JWRB)
              SSOURCE(IJ,K,M) = SSOURCE(IJ,K,M)/GTEMP1
            ENDDO
          ENDDO
        ENDDO
        CALL WNFLUXES (IJS, IJL,                                        &
     &                 MIJ, RHOWGDFTH,                                  &
     &                 SSOURCE, CICVR,                                  &
     &                 PHIWA,                                           &
     &                 EMEANALL, F1MEAN, U10NEW, THWNEW,                &
     &                 USNEW, ROAIRN, .TRUE.)
      ENDIF


!SHALLOW
      IF(ISHALLO.NE.1) THEN
        CALL SDIWBK(IJS, IJL, FL3 ,FL, SL, DEPTH(IJS,1), EMAXDPT(IJS), EMEANALL, F1MEAN)
        CALL SBOTTOM (IJS, IJL, DEPTH(IJS,1), FL3, FL, SL)
      ENDIF
!SHALLOW

! ----------------------------------------------------------------------

!*    2.4 COMPUTATION OF NEW SPECTRA.
!         ---------------------------

!     INCREASE OF SPECTRUM IN A TIME STEP IS LIMITED TO A FINITE
!     FRACTION OF A TYPICAL F**(-4) EQUILIBRIUM SPECTRUM.

      DO M=1,NFRE
        DELFL(M) = COFRM4(M)*DELT
      ENDDO
      DO IJ=IJS,IJL
        USFM(IJ) = USNEW(IJ)*MAX(FMEANWS(IJ),FMEANALL(IJ))
      ENDDO

      DO M=1,NFRE
        DO IJ=IJS,IJL
          TEMP(IJ,M) = USFM(IJ)*DELFL(M)
        ENDDO
      ENDDO

      IF (LLUNSTR) THEN
      DO K=1,NANG
        DO M=1,NFRE
          DO IJ=IJS,IJL
            GTEMP1 = MAX((1.0_JWRB-DELT5*FL(IJ,K,M)),1.0_JWRB)
            GTEMP2 = DELT*SL(IJ,K,M)/GTEMP1
            FLHAB = ABS(GTEMP2)
            FLHAB = MIN(FLHAB,TEMP(IJ,M))
            FL3(IJ,K,M) = FL3(IJ,K,M) + IOBND(IJ)*SIGN(FLHAB,GTEMP2)
            FL3(IJ,K,M) = MAX(IODP(IJ)*CIREDUC(IJ,K,M)*FL3(IJ,K,M),FLM(IJ,K))
            FL3(IJ,K,M) = MIN(FL3(IJ,K,M),FLMAX(M))
          ENDDO
        ENDDO
      ENDDO
      ELSE
      DO K=1,NANG
        DO M=1,NFRE
          DO IJ=IJS,IJL
            GTEMP1 = MAX((1.0_JWRB-DELT5*FL(IJ,K,M)),1.0_JWRB)
            GTEMP2 = DELT*SL(IJ,K,M)/GTEMP1
            FLHAB = ABS(GTEMP2)
            FLHAB = MIN(FLHAB,TEMP(IJ,M))
            FL3(IJ,K,M) = FL3(IJ,K,M) + SIGN(FLHAB,GTEMP2)
            FL3(IJ,K,M) = MAX(CIREDUC(IJ,K,M)*FL3(IJ,K,M),FLM(IJ,K))
            FL3(IJ,K,M) = MIN(FL3(IJ,K,M),FLMAX(M))
          ENDDO
        ENDDO
      ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    2.5 REPLACE DIAGNOSTIC PART OF SPECTRA BY A F**(-5) TAIL.
!         -----------------------------------------------------

      CALL FKMEAN(FL3, IJS, IJL, EMEANALL, FMEANALL, F1MEAN, AKMEAN, XKMEAN)

!     MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
      CALL FEMEANWS(FL3, IJS, IJL, EMEANWS, FMEANWS, XLLWS)

!     COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
      CALL FRCUTINDEX(IJS, IJL, FMEANALL, FMEANWS, USNEW, CICVR,        &
     &                MIJ, RHOWGDFTH)


      CALL IMPHFTAIL (IJS, IJL, MIJ, FLM, FL3)

! ----------------------------------------------------------------------

!*    2.6 SAVE WINDS INTO INTERMEDIATE STORAGE.
!         -------------------------------------

      DO IJ=IJS,IJL
        USOLD(IJ) = USNEW(IJ)
        Z0OLD(IJ) = Z0NEW(IJ)
        ROAIRO(IJ) = ROAIRN(IJ)
        WSTAROLD(IJ) = WSTARNEW(IJ)
      ENDDO

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('IMPLSCH',1,ZHOOK_HANDLE)

      END SUBROUTINE IMPLSCH
