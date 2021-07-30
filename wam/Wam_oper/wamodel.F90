SUBROUTINE WAMODEL (NADV, LDSTOP, LDWRRE)

! ----------------------------------------------------------------------

!**** *WAMODEL* - 3-G WAM MODEL - TIME INTEGRATION OF WAVE FIELDS.

!     S.D. HASSELMANN  MPI       1.12.85

!     G. KOMEN         KNMI         6.86  MODIFIED FOR SHALLOW WATER
!     P. JANSSEN                          ASPECTS.

!     S.D. HASSELMANN  MPI       15.2.87  MODIFIED FOR CYBER 205.

!     P. LIONELLO      ISDGM      6.3.87  MODIFIED TO OUTPUT SWELL.

!     S.D. HASSELMANN  MPI        1.6.87  ALL VERSIONS COMBINED INTO
!                                         ONE MODEL. DEEP AND SHALLOW
!                                         WATER , CRAY AND CYBER 205
!                                         VERSION.

!     CYCLE_2 MODICIFATIONS:
!     ----------------------

!     L. ZAMBRESKY     GKSS        10.87  OPTIMIZED FOR CRAY, CYBER 205
!     H. GUNTHER

!     A. SPEIDEL       MPI          4.88  VARIABLE DIMENSIONS, INTERNAL
!                                         CHECKS (CFL-CRITERION).

!     A. SPEIDEL       MPI         11.88  CHANGES FOR CRAY-2.

!     K. HUBBERT       POL          6.89  DEPTH AND CURRENT REFRACTION.
!                                         PRECALCULATION OF TERMS IN
!                                         *PROPDOT*.
!                                         SOLVE WAVE ACTION EQUATION
!                                         FOR CURRENT REFRACTION.

!     CYCLE_3 MODICIFATIONS:
!     ----------------------

!     R. PORTZ , S.D. HASSELMANN   MPI          1990

!      - RESTRUCTURE MODEL TO CALL THE ACTUAL INTEGRATION IN TIME
!        AS A SUBROUTINE: WAMODEL. A SHELL PROGRAM "WAMSHELL" READS
!        OUTPUT FROM PREPROC AND COMPUTES THE WIND ARRAYS FOR THE
!        INTEGRATION PERIOD FROM PREWIND, WHICH HAS BEEN INCORPORATED
!        AS A SUBROUTINE.
!      - ALL INTERMEDIATE AND RESTART I/O IS DONE IN THE SUBROUTINE
!        WAMODEL AND INPREST.
!      - THE YOWMON BLOCK IN THE PREPROCESSOR AND MODEL ARE MADE
!        COMPATIBLE.
!      - THE YOWPUTATION OF SEVERAL PARAMETERS HAS BEEN TRANSFERRED
!        FROM THE MODEL TO PREPROC.
!      - DEPTH AND CURRENT REFRACTION HAS BEEN INCORPORATED INTO THE
!        MODEL.
!      - OPEN BOUNDARIES ARE INCORPORATED IN THE MODEL.
!      - SEVERAL MINOR ERRORS HAVE BEEN REMOVED.
!      - THE BUFFERED I/O FOR THE CYBER 205 HAS BEEN CHANGED INTO A
!        BINARY READ AND WRITE.

!     CYCLE_4 MODICIFATIONS:
!     ----------------------

!     L. ZAMBRESKY   GKSS/ECMWF   6.89  ECMWF SUB VERSION
!                                       BASED ON CYCLE_2.

!     H. GUNTHER     GKSS/ECMWF 10.89  ECMWF SUB VERSION REORGANIZED.
!                                      - COMMON BLOCK STRUCTURE.
!                                      - BLOCKING STRUCTURE.
!                                      - TIME COUNTING.
!                                      - GRIDDED OUTPUT FIELDS.
!                                      - HEADERS ADDED TO OUTPUT FILES.
!                                      - ERRORS IN PROPAGATION CORRECTED

!     P.A.E.M. JANSSEN KNMI      1990  COUPLED MODEL.

!     H. GUNTHER     GKSS/ECMWF  8.91  LOGARITHMIC DEPTH TABLES.
!                                      MPI CYCLE_3 AND ECMWF VERSIONS
!                                      COMBINED INTO CYCLE_4.

!     J. BIDLOT ECMWF 1996   MESSAGE PASSING
!     J. DOYLE  ECMWF OCTOBER 1996   ATMOSPHERIC COUPLING
!     J. BIDLOT ECMWF FEBRUARY 1997   MODULE
!     J. BIDLOT ECMWF MARCH 1997  ADD SAVSTRESS AND SAVSPEC 
!     B. HANSEN ECMWF MARCH 1997  SIGNAL HANDLING.
!          LDSTOP* - SET .TRUE. IF STOP SIGNAL RECEIVED.
!          LDWRRE* - SET .TRUE. IF RESTART SIGNAL RECEIVED.
!     S. ABDALLA  ECMWF  OCTOBER 2000  INCLUDE AIR DENSITY & Zi/L
!                                      & CALL OUTBETA MOVED TO WAVEMDL
!

!SHALLOW
!          DIFFERENCES FOR SHALLOW WATER RUNS TO DEEP WATER RUNS
!          ARE ENCLOSED IN COMMENT LINES : 'CSHALLOW'.
!SHALLOW
!NEST
!          DIFFERENCES FOR NESTED GRID RUNS TO NORMAL RUNS
!          ARE ENCLOSED IN COMMENT LINES : 'CNEST'.
!NEST
!REFRA
!          DIFFERENCES FOR REFRACTION RUNS TO NORMAL RUNS
!          ARE ENCLOSED IN COMMENT LINES : 'CREFRA'.
!REFRA

!*    PURPOSE.
!     --------

!       COMPUTATION OF THE 2-D FREQUENCY-DIRECTION WAVE SPECTRUM AT ALL
!       GRID POINTS FOR A GIVEN INITIAL SPECTRUM AND FORCING SURFACE
!       STRESS FIELD.

!**   INTERFACE.
!     ----------

!     *CALL* *WAMODEL (NADV, LDSTOP, LDWRRE)*
!        *NADV*      INTEGER   NUMBER OF ADVECTION ITERATIONS
!                              PER CALL OF WAMODEL, OUTPUT PARAMETER.
!        *LDSTOP*    LOGICAL   SET .TRUE. IF STOP SIGNAL RECEIVED.
!        *LDWRRE*    LOGICAL   SET .TRUE. IF RESTART SIGNAL RECEIVED.

!     METHOD.
!     -------

!       GRID POINTS ARE LAT - LONG,VECTORIZATION IS ACHIEVED BY RUNNING
!       THROUGH THE GRID POINTS IN AN INNER LOOP ORGANIZED AS 1-D ARRAY
!       IN BLOCKS,-ALL COMPUTATIONS ARE CARRIED OUT FOR ONE BLOCK AT A
!       TIME (SEE "BLOCK STRUCTURE" BELOW)

!       ALL COMPONENTS OF THE SPECTRUM ARE YOWPUTED PROGNOSTICALLY FROM
!       THE SPECTRAL TRANSPORT EQUATION UP TO A VARIABLE CUT-OFF
!       FREQUENCY = MAX(4*FPM,2.5*FMEAN),WHERE FPM IS THE
!       PIERSON MOSKOVITZ FREQUENCY AND FMEAN IS THE MEAN FREQUENCY,
!       BEYOND THE PROGNOSTIC CUTOFF A DIAGNOSTIC F**-5 TAIL IS ATTACHED
!       CONTINUOUSLY FOR EACH DIRECTION,

!       SOURCE FUNCTIONS ARE TAKEN FROM KOMEN ET AL(1984)

!       THE NONLINEAR TRANSFER IS PARAMETERIZED BY THE DISCRETE INTER-
!       ACTION APPROXIMATION OF HASSELMANN ET AL (1985B)

!       THE SOURCE FUNCTION AND THE ADVECTION TERM ARE INTEGRATED ON TWO
!       DIFFERENT TIME STEP LEVELS AND WITH DIFFERENT METHODS,-THE
!       ADVECTION TIME STEP IS A MULTIPLE OF THE SOURCE FUNCTION
!       TIME STEP.

!       THE SOURCE FUNCTIONS ARE INTEGRATED IMPLICITLY ACCORDING TO
!       HASSELMANN AND HASSELMANN (1985A),-THE RELEVANT FUNCTIONAL
!       DERIVATIVES OF THE INDIVIDUAL SOURCE FUNCTIONS REQUIRED FOR THE
!       SOLUTION OF THE IMPLICIT EQUATION ARE YOWPUTED WITHIN THE SOURCE
!       FUNCTION SUBS,- THE TIME STEP IS TYPICALLY 20 MIN,

!       THE ADVECTION IS INTEGRATED BY A FIRST ORDER UPWIND SCHEME,ALSO
!       ACCORDING TO HASSELMANN AND HASSELMANN (1985A),-THE ADVECTIVE
!       TIMESTEP IS DEPENDENT ON THE FREQUENCY AND SPATIAL GRID IN
!       ACCORDANCE WITH CFL,

!       WINDS ARE READ IN EVERY WIND TIME STEP.IF THE WIND TIME STEP IS
!       GREATER THAN THE SOURCE TERM TIME STEP DELTWIND/DELTSOURCE STEPS
!       ARE INTEGRATED WITH CONSTANT WINDS,
!       WIND TIME STEP,PROPAGATION TIME STEP AND SOURCE TERM TIME STEP
!       SHOULD HAVE INTEGER RATIOS, THEY ARE GIVEN IN SECONDS AT
!       FULL MINUTES.

!NEST
!       ZERO ENERGY INFLUX IS ASSUMED AT COAST LINES. OPEN BOUNDARIES
!       ARE INCORPORATED IN THE MODEL, IF IT RUNS AS A NESTED GRID.
!NEST

!       BLOCK STRUCTURE (SEE PREPROC FOR DETAILS):
!       SEA POINTS ARE COLLECTED INTO A 1-DIMENSIONAL ARRAY.
!       BLOCKS OF MAXIMALLY NIBLO ELEMENTS.
!       SEA POINTS ARE COUNTED ALONG LINES OF LATITUDES FROM LEFT COAST
!       TO RIGHT COAST WORKING FROM SOUTH TO NORTH.
!       BLOCKS OVERLAP OVER TWO LATITUDE LINES,TO COMPUTE NORTH-SOUTH
!       ADVECTION TERMS, SEE ALSO COMMON GRIDPAR AND UBUF.

!       THE WIND FILES FOR THE BLOCKED WINDS CREATED BY PREWIND ARE
!       READ AND DELETED IN SUB IMPLSCH (IU17 AND IU18).


!       ALL PARAMETERS HAVE TO BE THE VALUES GIVEN AT THE END OF THE
!       PREPROC OUTPUT IN COLUMN 'REQUIRED'.

!     EXTERNALS.
!     ----------

!       *ABORT1*    - TERMINATES PROCESSING.
!       *AIRSEA*    - SURFACE LAYER STRESS.
!NEST
!       *BOUINPT*   - BOUNDARY VALUE INPUT.
!NEST
!REFRA
!       *DOTDC*     - READ COMMON REFDOT.
!REFRA
!       *FILLBL*    - ADD LATITUDES TO A BLOCK.
!       *GSFILE*    - ROUTINE TO DYNAMICALLY FETCH OR DISPOSE FILES.
!       *HEADBC*    - WRITE BOUNDARY OUTPUT FILE HEADER.
!       *IMPLSCH*   - IMPLICIT SCHEME FOR INTEGRATION OF SOURCE
!                     FUNCTIONS IN TIME AND INPUT OF WINDS.
!       *INCDATE*   - UPDATE DATE TIME GROUP.
!NEST
!       *INTSPEC*   - INTERPOLATION OF SPECTRA.
!NEST
!       *MAKEGRID*  - MAKES GRIDDED FIELDS.
!       *OUTBC*     - OUTPUT OF BOUNDARY VALUES.
!       *OUTBS*     - CONTROLS OUTPUT FROM BLOCKS.
!       *OUTGRID*   - SAVE BLOCKED PARAMETERS INTO GRID ARRAYS.
!       *OUTWINT*   - OUTPUT OF INTEGRATED PARAMETERS.
!       *PEAKFR*    - COMUTE PEAK FREQUENCY.
!       *PROPAG*    - PROPAGATION SCHEME.
!NEST
!       *ROTSPEC*   - ROTATE A SPECTRUM.
!NEST
!       *SAVSTRESS* - DISPOSE STRESS/WIND RESTART FILES.
!       *SAVSPEC    - DISPOSE SPECTRUM RESTART FILES.
!SHALLOW
!       *SBOTTOM*   - COMPUTES BOTTOM DISSIPATION SOURCE TERM AND
!                     LINEAR CONTRIBUTION TO FUNCTIONAL MATRIX.
!SHALLOW
!       *SDISSIP*   - COMPUTATION OF DISSIPATION SOURCE FUNCTION
!                     AND LINEAR CONTRIBUTION OF DISSIPATION TO
!                     FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
!       *SEMEAN*    - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.
!       *SEPWISW*   - COMPUTATION OF 2-DIMENSIONAL SWELL DISTRIBUTION
!                     TOTAL SWELL ENERGY, MEAN SWELL DIRECTION, AND
!                     MEAN SWELL FREQUENCY AT EACH GRID POINT.
!ICE
!       *SETICE*    - SET SPECTRA ON ICE EDGE TO ZERO.
!ICE
!       *SINFLX*    - COMPUTATION OF INPUT SOURCE FUNCTION, AND
!                     LINEAR CONTRIBUTION OF INPUT SOURCE FUNCTION
!                     TO FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
!                     and COMPUTATION OF WAVE STRESS.
!       *SNONLIN*   - COMPUTATION OF NONLINEAR TRANSFER RATE AND
!                     DIAGONAL LINEAR CONTRIBUTION OF NONLINEAR SOURCE
!                     FUNCTION TO FUNCTIONAL MATRIX.
!       *STHQ*      - COMPUTATION OF MEAN WAVE DIRECTION AT EACH
!                     GRID POINT.
!NEST
!       *STRSPEC*   - STRETCH A SPECTRUM.
!NEST
!       *OUTWNORM*  - COMPUTES A FEW NORMS OF GRIDDED FIELDS
!       *UNSETICE*  - FILL GAPS LEFT WHEN SEA ICE MASK CHANGES INITIALLY
!                     FROM THE PREVIOUS RUN.

!     REFERENCE.
!     ----------

!       SNYDER, R.L., F.W. DOBSON, J.A. ELLIOT, AND R.B. LONG:
!          ARRAY MEASUREMENTS OF ATMOSPHERIC PRESSURE FLUCTUATIONS
!          ABOVE SURFACE GRAVITY WAVES. J.FLUID MECH. 102, 1-59 ,1981.
!       G. KOMEN, S. HASSELMANN, K. HASSELMANN:
!          ON THE EXISTENCE OF A FULLY DEVELOPED WIND SEA SPECTRUM.
!          JPO,1984.
!       S. HASSELMANN, K. HASSELMANN, J.H. ALLENDER, T.P. BARNETT:
!          IMPROVED METHODS OF COMPUTING AND PARAMETERIZING THE
!          NONLINEAR ENERGY TRANSFER IN A GRAVITY WAVE SPECTRUM.
!          JPO, 1985.
!       S. HASSELMANN, K. HASSELMANN: A GLOBAL WAVE MODEL,
!          WAM REPORT,JUNE,30/1985.
!       P. JANSSEN, G. KOMEN: A SHALLOW WATER EXTENSION OF THE
!          3-G WAM-MODEL. WAM REPORT 1985.
!       THE WAMDI GROUP: THE WAM MODEL - A THIRD GENERATION OCEAN
!          WAVE PREDICTION MODEL. JPO, VOL. 18, NO. 12, 1988.
!       P.A.E.M JANSSEN: JPO, 1989 AND 1991.
!       K. HASSELMANN: TRANSPORT EQUATION OF FINITE DEPTH SURFACE
!          WAVE SPECTRUM IN TIME DEPENDANT CURRENT AND DEPTH FIELD USING
!          NONCANONICAL SPATIAL (SPHERICAL) AND WAVE NUMBER (FRQUENCY-
!          DIRECTION) COORDINATES. WAM REPORT 1988.

! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCPBO  , ONLY : IBOUNC   ,NBOUNC    ,GBOUNC  , IPOGBO  ,    &
     &            CBCPREF
      USE YOWCOUP  , ONLY : LWCOU    ,LLNORMWAMOUT,                     &
     &                      KCOUSTEP, LWNEMOCOU, NEMONTAU,              &
     &                      LWNEMOCOUSTK ,LWNEMOCOUSTRN,                &
     &                      NEMOWSTEP, NEMOFRCO     ,                   &
     &                      NEMOCSTEP, NEMONSTEP   ,                    &
     &                      NEMOSTRN, NEMOUSTOKES, NEMOVSTOKES
      USE YOWCOUT  , ONLY : COUTT    ,COUTS    ,FFLAG20  ,GFLAG20  ,    &
     &            JPPFLAG  ,FFLAG    ,GFLAG    ,                        &
     &            IRCD     ,IRU10    , IRALTHS  ,IRALTHSC  ,IRALTRC ,   &
     &            NGOUT    ,LLOUTERS ,                                  &
     &            NIPRMOUT ,                                            &
     &            LFDB     ,NOUTT    ,NOUTS    ,                        &
     &            CASS     ,NASS     ,LOUTINT  ,                        &
     &            LRSTPARALW, LRSTINFDAT,                               &
     &            LRSTST0  ,LWAMANOUT
      USE YOWCURR  , ONLY : CDTCUR
      USE YOWFPBO  , ONLY : IBOUNF
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIBHD, ONLY : LGRHDIFS 
      USE YOWGRID  , ONLY : IJS      ,IJL       ,                       &
     &            IJSLOC   ,IJLLOC   ,IJGLOBAL_OFFSET
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE ,CICOVER  ,CITHICK  ,    &
     &            CIWA
      USE YOWMEAN  , ONLY : WSEMEAN  ,WSFMEAN  ,USTOKES  ,VSTOKES  ,STRNMS
      USE YOWMESPAS, ONLY : LFDBIOOUT,LGRIBOUT ,LNOCDIN  ,LWAVEWIND 
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,KTAG 
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZMISS    ,DEG      ,EPSMIN
      USE YOWSHAl  , ONLY : DEPTH    ,EMAXDPT
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDTPRO   ,CDTRES   ,    &
     &            CDATER   ,CDATES   ,CDTINTT  ,IDELPRO  ,IDELT    ,    &
     &            IDELWI   ,IREST    ,IDELRES  ,IDELINT  ,              &
     &            ISHALLO  ,IASSI    ,                                  &
     &            CDTBC    ,IDELBC   ,                                  &
     &            IPROPAGS ,                                            &
     &            NENSFNB  ,NTOTENS  ,NSYSNB   ,                        &
     &            NMETNB   ,CDATEA   ,MARSTYPE ,YCLASS   ,YEXPVER  ,    &
     &            LLSOURCE ,                                            &
     &            LANAONLY ,LFRSTFLD ,NPROMA_WAM,IREFDATE
      USE YOWSPEC, ONLY   : NBLKS    ,NBLKE    ,                        &
     &            U10NEW   ,U10OLD   ,THWNEW   ,THWOLD   ,USNEW    ,    &
     &            USOLD    ,Z0NEW    ,Z0OLD    ,TAUW     ,TAUWDIR  ,    &
     &            Z0B      ,                                            &
     &            BETAOLD  ,ROAIRN   ,ROAIRO   ,ZIDLNEW  ,ZIDLOLD  ,    &
     &            FL1
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH    ,CWI      ,LRESTARTED
      USE YOWUNIT  , ONLY : IU02     ,IU04     ,              &
     &            IU19     ,IU20     ,IU30
      USE YOWWAMI  , ONLY : CBPLTDT  ,CEPLTDT  ,IANALPD  ,IFOREPD  ,    &
     &            IDELWIN  ,NFCST    ,ISTAT
      USE YOWWIND  , ONLY : CDATEWO
      USE YOWNEMOP , ONLY : NEMODP
      USE UNWAM, ONLY : EXCHANGE_FOR_FL1
      USE YOWUNPOOL ,ONLY : LLUNSTR, LLUNBINOUT
      USE YOWPD, ONLY : MNP => npa
      USE MPL_MODULE
      USE FDBSUBS_MOD, ONLY : IFLUSHFDBSUBS
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWNODEPOOL, ONLY : NP, IPLG, IPGL
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal
      
! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "mpcrtbl.intfb.h"
#include "outwint.intfb.h"
#include "outwpsp.intfb.h"
#include "abort1.intfb.h"
#include "bouinpt.intfb.h"
#include "chesig.intfb.h"
#include "closend.intfb.h"
#include "difdate.intfb.h"
#include "gsfile_new.intfb.h"
#include "headbc.intfb.h"
#include "implsch.intfb.h"
#include "incdate.intfb.h"
#include "newwind.intfb.h"
#include "outbc.intfb.h"
#include "outbs.intfb.h"
#include "outint.intfb.h"
#include "outunwamtest.intfb.h"
#include "outwnorm.intfb.h"
#include "outspec.intfb.h"
#include "propag_wam.intfb.h"
#include "savspec.intfb.h"
#include "savstress.intfb.h"
#include "setice.intfb.h"
#include "unsetice.intfb.h"
#include "updnemofields.intfb.h"
#include "updnemostress.intfb.h"
#include "wdfluxes.intfb.h"
#include "writsta.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: NADV
      INTEGER(KIND=JWIM) :: IJ, K, M, J, IRA, KADV, ICH, ICALL
      INTEGER(KIND=JWIM) :: IFIL, IC, ICL, ICR, II
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA
      INTEGER(KIND=JWIM) :: JSTPNEMO, IDATE, ITIME
      INTEGER(KIND=JWIM) :: IWAM_GET_UNIT
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL) :: MIJ

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: XLLWS

      CHARACTER(LEN= 2) :: MARSTYPEBAK
      CHARACTER(LEN=14) :: CDATEWH, CZERO
      CHARACTER(LEN=14) :: CDTINTTBAK
      CHARACTER(LEN=14) :: CDATE, CDTPRA, CDTIMP, CDTIMPNEXT, CDTRCF

      LOGICAL, INTENT(INOUT) :: LDSTOP, LDWRRE
      LOGICAL :: LLFLUSH
      LOGICAL :: LSV, LRST, LOUT
      LOGICAL :: LDREPROD
      LOGICAL :: NEWREAD
      LOGICAL, SAVE :: NEWFILE
      LOGICAL :: LLNONASSI
      LOGICAL :: LL2NDCALL
      LOGICAL :: FFLAGBAK(JPPFLAG), GFLAGBAK(JPPFLAG)
      LOGICAL,ALLOCATABLE,DIMENSION(:) :: LCFLFAIL

      DATA NEWFILE / .FALSE. /

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WAMODEL',0,ZHOOK_HANDLE)


      IF (ITEST .GE. 2) THEN
        WRITE(IU06,*)
        WRITE(IU06,*) '   SUB. WAMODEL: JUST ENTERED FOR CDTPRO= ',CDTPRO, ' AND CDATEA= ', CDATEA
        CALL FLUSH(IU06)
      ENDIF

      CZERO = ' '
      LLFLUSH = .FALSE.
      KTAG=200
      LRSTST0=.FALSE.
      NPROMA=NPROMA_WAM

!     TIME FOR THE NEXT SOURCE TERM INTEGRATION
      CDTIMPNEXT=CDTPRO
      CALL INCDATE(CDTIMPNEXT,IDELT)   
!     TIME FOR WIND INPUT UPDATE (SEE NEWWIND)
      CDTIMP=CDTPRO

!*    0. WRITE OUT NORMS FOR WAVE ENERGY, PEAK PERIOD AND WIND SPEED
!        -----------------------------------------------------------

      IF(CDTPRO.EQ.CDATEA .AND. LLSOURCE ) THEN
        IF (ITEST.GE.1) THEN
           WRITE(IU06,*) '   SEEDING OF ZERO ACTION'
           FLUSH(IU06)
        ENDIF
!       INSURE THERE IS SOME WAVE ENERGY FOR GRID POINTS THAT HAVE BEEN
!       FREED FROM SEA ICE (ONLY DONE INITIALLY AND IF THE MODEL IS NOT
!       RESTARTED).
!       IT ALSO RESETS THE MIMIMUM ENERGY LEVEL THAT MIGHT HAVE BEEN LOST
!       WHEN GETTING THE DATA FROM GRIB.
        CALL GSTATS(1236,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          CALL UNSETICE(IJS, IJL, KIJS, KIJL,                      &
     &                  DEPTH(KIJS), EMAXDPT(KIJS),                &
     &                  CICOVER(KIJS), U10OLD(KIJS), THWOLD(KIJS), &
     &                  FL1)
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1236,1)
      ENDIF

      IF(CDTPRO.EQ.CDATEA .OR. CDTPRO.EQ.CDATEF) THEN

        CDTINTTBAK=CDTINTT
        CDTINTT=CDTPRO

!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ) 
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          DO IJ=KIJS,KIJL
            U10NEW(IJ) = U10OLD(IJ)
            THWNEW(IJ) = THWOLD(IJ)
            USNEW(IJ) = USOLD(IJ)
            Z0NEW(IJ) = Z0OLD(IJ)
            ROAIRN(IJ) = ROAIRO(IJ)
            ZIDLNEW(IJ) = ZIDLOLD(IJ)
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO

        IF(LLSOURCE) THEN
!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
          DO JKGLO=IJS,IJL,NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL)
            CALL WDFLUXES (IJS, IJL, KIJS, KIJL,                        &
     &                     MIJ(KIJS),                                   &
     &                     FL1, XLLWS(KIJS:KIJL,:,:),    &
     &                     DEPTH(KIJS),                                 &
     &                     CICOVER(KIJS),                               &
     &                     U10NEW(KIJS), THWNEW(KIJS), USNEW(KIJS),     &
     &                     Z0NEW(KIJS), Z0B(KIJS),                      &
     &                     ROAIRN(KIJS), ZIDLNEW(KIJS),                 &
     &                     WSEMEAN(KIJS), WSFMEAN(KIJS),                &
     &                     USTOKES(KIJS), VSTOKES(KIJS), STRNMS(KIJS) )
          ENDDO
!$OMP     END PARALLEL DO
          IF (ITEST.GE.2) THEN
            WRITE(IU06,*) '   SUB. WAMODEL: WDFLUXES CALLED'
           CALL FLUSH (IU06)
          ENDIF

          IF (LLUNSTR) CALL EXCHANGE_FOR_FL1(FL1)

        ELSE
          MIJ(:) = NFRE
          XLLWS(:,:,:)=0._JWRB
        ENDIF

!       SET FL1 ON ICE POINTS TO ZERO
        IF (LICERUN .AND. LMASKICE .AND. LLSOURCE) THEN
          IF (ITEST.GE.1) THEN
            WRITE(IU06,*) '   SUB. WAMODEL: INITIAL SPECTRUM = 0 AT ICE POINTS'
          ENDIF
! Mod for OPENMP
          CALL GSTATS(1439,0)
!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
          DO JKGLO=IJS,IJL,NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL)
            CALL SETICE(IJS, IJL, KIJS, KIJL, FL1,                 &
     &                  CICOVER(KIJS), U10NEW(KIJS), THWNEW(KIJS))
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1439,1)
        ENDIF

!       OUTPUT POINT SPECTRA
        IF(NGOUT.GT.0 .OR. LLOUTERS) THEN
          CALL OUTWPSP (IJSLOC, IJLLOC, IJGLOBAL_OFFSET,                &
     &                  FL1(IJSLOC:IJLLOC,:,:),USNEW(IJSLOC))
        ENDIF

!*      0.1 SAVE INITIAL INTEGRATED FIELDS
!           ------------------------------
        IF((MARSTYPE.EQ.'cf' .OR. MARSTYPE.EQ.'pf' .OR.                 &
     &      MARSTYPE.EQ.'fc' .OR. MARSTYPE.EQ.'4v' .OR.                 &
     &      LANAONLY         .OR. LFRSTFLD             )                &
     &     .AND. LWAMANOUT) THEN

          IF(LFRSTFLD) THEN
            FFLAGBAK=FFLAG
            GFLAGBAK=GFLAG
            MARSTYPEBAK=MARSTYPE
            IF(MARSTYPE.EQ.'fg'.OR.MARSTYPE.EQ.'an') THEN
              IF(.NOT.LNOCDIN) THEN
                FFLAG(IRCD)=.FALSE.
                GFLAG(IRCD)=.FALSE.
              ENDIF
              IF(LWAVEWIND) THEN
                FFLAG(IRU10)=.FALSE.
                GFLAG(IRU10)=.FALSE.
              ENDIF
              GFLAG(IRALTHS)=.FALSE.
              GFLAG(IRALTHSC)=.FALSE.
              GFLAG(IRALTRC)=.FALSE.
              MARSTYPE='an'
            ENDIF
!           change the output flags need to reset mapping
            CALL MPCRTBL
          ENDIF
!         NEED TO TEMPORARLY RESET THE IFS FORECAST STEP
!         IF THE IFS GRIB HEADER IS USED SUCH THAT IT POINTS TO
!         THE START OF THE RUN.
          IF(LGRHDIFS) THEN
             LRSTST0=.TRUE.
          ENDIF

!         COMPUTE OUTPUT PARAMETERS
          IF(NIPRMOUT.GT.0) THEN
            CALL OUTBS (IJS, IJL, MIJ, FL1, XLLWS)
!           PRINT OUT NORMS
            !!!1 to do: decide if there are cases where we might want LDREPROD false
            LDREPROD=.TRUE.
            CALL OUTWNORM(LDREPROD)
          ENDIF

          IF( .NOT. LRESTARTED ) THEN
            IF(IREST.EQ.1 .AND. MARSTYPE.NE.'an' .AND. LGRIBOUT) THEN
              CALL OUTSPEC(IJS, IJL, FL1, CICOVER)
              LLFLUSH = .TRUE.
            ENDIF

            IF(NIPRMOUT.GT.0 ) THEN
              CALL OUTWINT
              LLFLUSH = .TRUE.
            ENDIF

          ENDIF

          IF (LFDB.AND.LLFLUSH) THEN
             CALL GSTATS(1976,0)
             CALL IFLUSHFDBSUBS()
             CALL GSTATS(1976,1)
             LLFLUSH=.FALSE.
          ENDIF

          IF(LFRSTFLD) THEN
            MARSTYPE=MARSTYPEBAK
            FFLAG=FFLAGBAK
            GFLAG=GFLAGBAK
            LFRSTFLD=.FALSE.
!           change the output flags need to reset mapping
            CALL MPCRTBL
          ENDIF

          CDTINTT=CDTINTTBAK

          IF(LGRHDIFS) THEN
            LRSTST0=.FALSE.
          ENDIF

          IF (ITEST .GE. 2) THEN
            WRITE(IU06,*)
            WRITE(IU06,*) '   SUB. WAMODEL: INITIAL FIELDS SAVED '
            WRITE(IU06,*) '                 FOR FORECAST RUN     '
            WRITE(IU06,*) '    '
            CALL FLUSH(IU06)
          ENDIF
      
        ENDIF

      ENDIF ! DO_UPDATE_WIND_FLUX




!*    1. ADVECTION/PHYSICS TIME LOOP.
!        ----------------------------

      ADVECTION : DO KADV = 1,NADV

!*    1.1 FIX END DATE OF THIS PROPAGATION STEP AND OUTPUT TIMES.
!         -------------------------------------------------------

        CDTPRA = CDTPRO
        CALL INCDATE(CDTPRO,IDELPRO)
        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '   SUB. WAMODEL: START OF PROPAGATION '
          WRITE(IU06,*) '     START DATE IS    CDTPRA = ',CDTPRA
          WRITE(IU06,*) '     END DATE WILL BE CDTPRO = ',CDTPRO
          CALL FLUSH (IU06)
        ENDIF

!         UPDATE OUTPUT TIMES.

        IF (NOUTT.GT.0) THEN
          CDTINTT = CZERO
          DO J=1,NOUTT
            IF (CDTPRO.EQ.COUTT(J)) THEN
              IF (FFLAG20.OR.GFLAG20) CDTINTT = COUTT(J)
            ENDIF
          ENDDO
        ELSE
          IF ((FFLAG20.OR.GFLAG20) .AND. CDTINTT.LT.CDTPRO) CALL INCDATE (CDTINTT,IDELINT)
        ENDIF

!       UPDATE SPECTRA OUTPUT DATE
        IF (NOUTS.GT.0) THEN
!         reset CDATES to insure that spectra output is only controlled
!         by list COUTS
          CDATES='000000000000'
          DO J=1,NOUTS
            IF (CDTPRO.EQ.COUTS(J)) THEN
              CDTRES=CDTPRO
              CDATER=CDTRES
              CDATES=CDTRES
              EXIT
            ENDIF
          ENDDO
        ELSE
          IF (CDTRES.LT.CDTPRO) CALL INCDATE(CDTRES,IDELRES)
        ENDIF

        IF ((IBOUNC.EQ.1.OR.IBOUNF.EQ.1).AND.CDTBC.LT.CDTPRO) CALL INCDATE(CDTBC,IDELBC)


!*    1.5.5 SET TIME COUNTER.
!           -----------------
          CDATE   = CDTPRA
          CDATEWH = CDATEWO    
          NEWREAD = .FALSE.

 1550     CONTINUE                       

!*    1.5.5.2 COMPUTATION OF PROPAGATION
!*            INTEGRATION OF SOURCE TERMS OVER SUB TIME STEPS BETWEEN
!*            PROPAGATION TIME STEPS.
!             -------------------------------------------------------

! IF CDATE CORRESPONDS TO A PROPAGATION TIME
          IF (CDATE.EQ.CDTPRA) THEN

            CALL PROPAG_WAM(IJS, IJL, FL1)

            IF (ITEST.GE.2) THEN
              WRITE(IU06,*) '   SUB. WAMODEL: PROPAGATION DONE'
              CALL FLUSH (IU06)
            ENDIF

            CDATE=CDTPRO   
          ENDIF

!*        RETRIEVING NEW FORCING FIELDS FROM TEMPORARY STORAGE IF NEEDED.
!         ---------------------------------------------------------------
          CALL NEWWIND(IJS,IJL,CDTIMP,CDATEWH,                          &
     &                 NEWREAD,NEWFILE,U10OLD,THWOLD,U10NEW,THWNEW,     &
     &                 USOLD, USNEW,                                    &
     &                 ROAIRO, ROAIRN, ZIDLOLD, ZIDLNEW,                &
     &                 CICOVER, CITHICK, CIWA,                          &
     &                 TAUW, BETAOLD)

          IF (ITEST.GE.2) THEN
            WRITE(IU06,*) '   SUB. WAMODEL: NEWWIND CALLED' 
            CALL FLUSH (IU06)
          ENDIF

!         IT IS TIME TO INTEGRATE THE SOURCE TERMS
!         ----------------------------------------
          IF(CDATE.GE.CDTIMPNEXT) THEN
!           COMPUTE UPDATE DUE TO SOURCE TERMS
            CALL GSTATS(1431,0)
            IF(LLSOURCE) THEN
!$OMP         PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
              DO JKGLO=IJS,IJL,NPROMA
                KIJS=JKGLO
                KIJL=MIN(KIJS+NPROMA-1,IJL)
                CALL IMPLSCH (FL1, IJS, IJL, KIJS, KIJL,                &
     &                        DEPTH(KIJS), EMAXDPT(KIJS),               &
     &                        THWOLD(KIJS), USOLD(KIJS),                &
     &                        TAUW(KIJS), TAUWDIR(KIJS),                &
     &                        Z0OLD(KIJS),                              &
     &                        ROAIRO(KIJS), ZIDLOLD(KIJS),              &
     &                        CICOVER(KIJS), CIWA(KIJS:KIJL,:),         &
     &                        U10NEW(KIJS), THWNEW(KIJS), USNEW(KIJS),  &
     &                        Z0NEW(KIJS), Z0B(KIJS),                   &
     &                        ROAIRN(KIJS), ZIDLNEW(KIJS),              &
     &                        WSEMEAN(KIJS), WSFMEAN(KIJS),             &
     &                        USTOKES(KIJS), VSTOKES(KIJS),STRNMS(KIJS),&
     &                        MIJ(KIJS), XLLWS(KIJS:KIJL,:,:))

              ENDDO
!$OMP       END PARALLEL DO

              IF (LWNEMOCOU) THEN
                NEMONTAU = NEMONTAU + 1
              ENDIF

              IF (ITEST.GE.2) THEN
                WRITE(IU06,*) '   SUB. WAMODEL: IMPLSCH CALLED'
                CALL FLUSH (IU06)
              ENDIF

            ELSE
!             NO SOURCE TERM CONTRIBUTION
!$OMP         PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ,K,M)  
              DO JKGLO=IJS,IJL,NPROMA
                KIJS=JKGLO
                KIJL=MIN(KIJS+NPROMA-1,IJL)
                DO IJ=KIJS,KIJL
                  MIJ(IJ) = NFRE
                  DO K=1,NANG
                    DO M=1,NFRE
                      FL1(IJ,K,M) = MAX(FL1(IJ,K,M),EPSMIN)
                      XLLWS(IJ,K,M) = 0.0_JWRB
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
!$OMP         END PARALLEL DO
            ENDIF
            CALL GSTATS(1431,1)


!*          COPY AND CLOSE FILES IF NEEDED.
!           -------------------------------

            CALL CLOSEND(IJS,IJL,CDTIMP,CDATEWH,                        &
     &                   NEWREAD,NEWFILE,U10OLD,THWOLD,ROAIRO,ZIDLOLD,  &
     &                   U10NEW,THWNEW,ROAIRN,ZIDLNEW)
            IF (ITEST.GE.2) THEN
              WRITE(IU06,*) '   SUB. WAMODEL: CLOSEND CALLED' 
              CALL FLUSH (IU06)
            ENDIF

            CDTIMP=CDTIMPNEXT
            CALL INCDATE(CDTIMPNEXT,IDELT)   

            IF (ITEST.GE.2) THEN
              WRITE(IU06,*) '   SUB. WAMODEL: SOURCE FUNCTIONS INTEGRATED:'
              CALL FLUSH(IU06)
            ENDIF
          ENDIF

          IF (LLUNSTR .AND. LLUNBINOUT) THEN
            CALL OUTUNWAMTEST (FL1(1:MNP,:,:), .FALSE.)
          ENDIF

!*    1.5.5.4 UPDATE TIME;IF TIME LEFT BRANCH BACK TO 1.5.5 
!             ---------------------------------------------
          IF (CDTIMPNEXT.LE.CDTPRO) GO TO 1550

!        END OF TIME LOOP ALL TIME STEPS DONE

          IF (ITEST.GE.2) THEN
            WRITE(IU06,*) '   SUB. WAMODEL: TIME STEPS DONE.'
            CALL FLUSH(IU06)
          ENDIF

!NEST
!*    1.5.7 INPUT OF BOUNDARY VALUES.
!           -------------------------

          CALL GSTATS(1909,0)
          IF (IBOUNF.EQ.1) THEN
            CALL BOUINPT (IU02, FL1, IJS, IJL, NBLKS, NBLKE)
          ENDIF
!NEST


!*    1.5.9 OUTPUT OF BOUNDARY POINTS.
!           --------------------------

          IF (IBOUNC.EQ.1) THEN
            CALL OUTBC (FL1, IJS, IJL, IU19)
          ENDIF
!NEST

!*    1.5.10 MODEL OUTPUT INTEGRATED DATA ARE SAVED
!            --------------------------------------

          IF (ITEST.GE.2) THEN
              WRITE(IU06,*) '   SUB. WAMODEL: MODEL ', 'OUTPUT CDTINTT=',CDTINTT
          ENDIF

#ifdef ECMWF
          IF(.NOT.LWCOU .AND. .NOT. LDSTOP) THEN
!!!!      the call to CHESIG is a signal handeling facility which is
!!!!      specific to running WAM at ECMWF, it can be ignored when
!!!!      WAM is not run at ECMWF.
            CALL CHESIG (IU06, ITEST, IRANK, NPROC, LDSTOP, LDWRRE)
          ENDIF
#endif

          LRST=(LDWRRE .AND. KADV.EQ.NADV )
          IF(LRST) THEN
            WRITE(IU06,*) ' '
            WRITE(IU06,*) '  ******************************************'
            IF(LDSTOP) THEN
            WRITE(IU06,*) '  AN INTERRUPT SIGNAL HAS BEEN RECEIVED '
            ENDIF
            WRITE(IU06,*) '  THE NECESSARY BINARY RESTART FILES WILL BE'
            WRITE(IU06,*) '  GENERATED.'
            WRITE(IU06,*) '  ******************************************'
            WRITE(IU06,*) ' '
            CALL FLUSH (IU06)
          ENDIF
          CALL GSTATS(1909,1)

          IF ( CDTINTT.EQ.CDTPRO .OR. LRST ) THEN

!           OUTPUT POINT SPECTRA
            IF(NGOUT.GT.0 .OR. LLOUTERS) THEN
              CALL OUTWPSP (IJSLOC, IJLLOC, IJGLOBAL_OFFSET,            &
     &                      FL1(IJSLOC:IJLLOC,:,:),USNEW(IJSLOC))
            ENDIF

!           COMPUTE OUTPUT PARAMETERS
            IF(NIPRMOUT.GT.0) THEN
              CALL OUTBS (IJS, IJL, MIJ, FL1, XLLWS)

!!!1 to do: decide if there are cases where we might want LDREPROD false
              LDREPROD=.TRUE.
              IF(LLNORMWAMOUT) CALL OUTWNORM(LDREPROD)

            ENDIF

            IF (ITEST.GE.2) THEN
              WRITE(IU06,*) '   SUB. WAMODEL: MODEL OUTPUT PREPARED IN OUTBS'
              CALL FLUSH (IU06)
            ENDIF

          ENDIF



!*    1.7 ONE PROPAGATION TIMESTEP DONE FOR ALL BLOCKS.
!         ---------------------------------------------

        WRITE(IU06,*) ' !!!!!!!!!!!!!! WAVE FIELDS INTEGRATED FOR DATE : ', CDTPRO
        IF (ITEST.GE.1) CALL FLUSH (IU06)

        LLNONASSI=.TRUE.
        IF (IASSI.EQ.1) THEN
!         IS THIS AN ANALYSIS TIME FROM THE INPUT LIST ?
          DO J=1,NASS
            IF(CDTPRO.EQ.CASS(J)) THEN
              LLNONASSI=.FALSE.
              EXIT
            ENDIF
          ENDDO
        ENDIF


!*      SAVE 2D SPECTRA AND/OR RESTART FIELDS.
!       --------------------------------------
        LSV=(CDTRES.EQ.CDTPRO.OR.CDATEE.EQ.CDTPRO.OR.CDTPRO.EQ.CDATER)

        IF (LSV.OR.LRST) THEN

!         THIS WILL HAPPEN WHEN IT IS NOT IN DATA ASSIMILATION MODE AND
!         IT IS EITHER A DETERMINED OUTPUT TIME
!         OR THE INTERUPT SIGNAL HAS BEEN TRIGGERED and it will wait
!         until the end of the advection loop.
!         OTHERWISE THE OUTPUT WILL OCCUR IN WAMASSI.

          LOUT= ((IREST.EQ.1) .AND. (CDTPRO.EQ.CDATER .OR. CDTPRO.LE.CDATES)) .AND. LSV .AND. LWAMANOUT

          IF ( LOUT .OR. LRST ) THEN

!           SAVE SPECTRUM IN GRIB
            IF(LOUT.AND.LGRIBOUT) THEN
!             we have insured that the spectra will be written to FDB
!             even when the restart option is triggered and it is an
!             output step for the spectra.

!             IF THE OUTPUT TIME IS NOT AN ANALYSIS TIME THEN TYPE FG or 4V
!             BECOMES TYPE AN (i.e. pseudo analysis)
              MARSTYPEBAK=MARSTYPE
              IF((MARSTYPE.EQ.'fg' .AND. KADV.LT.NADV) .OR.             &
     &           (MARSTYPE.EQ.'4v' .AND. LLNONASSI) ) THEN
                MARSTYPE='an'
              ENDIF

              CALL OUTSPEC(IJS, IJL, FL1, CICOVER)
              LLFLUSH = .TRUE.

              MARSTYPE=MARSTYPEBAK

              WRITE(IU06,*) ' '
              WRITE(IU06,*) '  GRIB WAVE SPECTRA DISPOSED AT........ CDTPRO  = ', CDTPRO
              WRITE(IU06,*) ' '
              IF (ITEST.GE.1) CALL FLUSH (IU06)
            ENDIF  

!           SAVE RESTART FILES IN PURE BINARY FORM
            IF ( .NOT.LGRIBOUT .OR. LDWRRE ) THEN
              CALL SAVSTRESS(U10OLD, THWOLD, USOLD, TAUW, TAUWDIR,     &
     &                       Z0OLD, ROAIRO, ZIDLOLD, CICOVER, CITHICK, &
     &                       NBLKS, NBLKE, CDTPRO, CDATEF)
              WRITE(IU06,*) ' '
              WRITE(IU06,*) '  BINARY STRESS FILE DISPOSED AT........',  &
     &         ' CDTPRO  = ', CDTPRO
              WRITE(IU06,*) ' '
              CALL SAVSPEC(IJS, IJL, FL1, NBLKS, NBLKE, CDTPRO, CDATEF, CDATER)
              WRITE(IU06,*) '  BINARY WAVE SPECTRA DISPOSED AT........', &
     &         ' CDTPRO  = ', CDTPRO
              WRITE(IU06,*) ' '
              CALL FLUSH(IU06)
            ENDIF


!*          UPDATE, WRITE AND SAVE WAMINFO FILE.
            IF (LRST .AND. IRANK.EQ.1) THEN
              ICH = 7 
              CALL DIFDATE (CDATEF,CDATEE,IFOREPD)
              IF (CDTPRO.LE.CDATEF) THEN
                CALL DIFDATE (CDTPRO,CDATEF,IANALPD)
                CBPLTDT = CDTPRO
                NFCST = 1
              ELSE
                NFCST = 0
                IANALPD = 0
                CBPLTDT = CDATEF
                CALL DIFDATE (CDTPRO,CDATEE,IFOREPD)
              ENDIF
              ISTAT(:) = 0
              IF (CDATE.EQ.CDATEE) ISTAT(1) = 1
              IDELWIN = IDELWI

              CEPLTDT = CDATEF

              IU04 = IWAM_GET_UNIT (IU06,CWI(1:ICPLEN+8) , 'w', 'f', 0)

              CALL WRITSTA (IU04, CDTPRO, CDATEE, IANALPD, IFOREPD,     &
     &                      IDELWIN, CDATER, CDATES, CBPLTDT, CEPLTDT,  &
     &                      IASSI, NFCST, ISTAT, CDTCUR,                &
     &                      LRSTPARALW, NPROC)
 
              CLOSE (IU04)
              WRITE(IU06,*) ' WAMINFO FILE WRITTEN FOR RESTART...',     &
     &         ' CDTPRO  = ', CDTPRO
              WRITE(IU06,*) '                                    ',     &
     &         ' CDATEF  = ', CDATEF
              WRITE(IU06,*) ' TO ', CWI(1:ICPLEN+8)
              CALL FLUSH(IU06)

              IF (LRSTINFDAT) THEN
                CDTRCF=CDTPRO
                CALL INCDATE(CDTRCF,IDELPRO)
                IU04 =  IWAM_GET_UNIT (IU06,CWI(1:ICPLEN+8)//'.'//      &
     &                              CDTRCF(1:8)//'_'//CDTRCF(9:14),     &
     &                              'w', 'f', 0)
                CALL WRITSTA (IU04, CDTPRO, CDATEE, IANALPD, IFOREPD,   &
     &                        IDELWIN, CDATER, CDATES, CBPLTDT, CEPLTDT,&
     &                        IASSI, NFCST, ISTAT, CDTCUR,              &
     &                        LRSTPARALW, NPROC)
                CLOSE (IU04)
              ENDIF

            ENDIF

          ENDIF  ! END SAVE RESTART FIELDS

        ENDIF  ! END SAVE


!       WRITE INTEGRATED DATA TO FDB OR TO FILE
!       ---------------------------------------
        IF (LWAMANOUT .AND. CDTINTT.EQ.CDTPRO .AND. NIPRMOUT.GT.0 ) THEN
!         IF THE OUTPUT TIME IS NOT AN ANALYSIS TIME THEN TYPE FG or 4V
!         BECOMES TYPE AN (i.e. speudo analysis)
          MARSTYPEBAK=MARSTYPE
          IF((MARSTYPE.EQ.'fg' .AND. KADV.LT.NADV) .OR.                 &
     &       (MARSTYPE.EQ.'4v' .AND. LLNONASSI) ) THEN
            MARSTYPE='an'
          ENDIF

          CALL OUTWINT
          LLFLUSH = .TRUE.
          IF (ITEST.GE.1) CALL FLUSH (IU06)

          MARSTYPE=MARSTYPEBAK

          CALL GSTATS(753,0)
          CALL MPL_BARRIER(CDSTRING='WAMODEL:')
          CALL GSTATS(753,1)
        ENDIF

!*      FLUSH FDB IF IT HAS BEEN USED AND IT IS NOT AN ANALYSIS
!       -------------------------------------------------------

        IF( LFDB .AND. LLFLUSH .AND. (IASSI.NE.1 .OR. CDTPRO.GT.CDATEF) ) THEN
          CALL GSTATS(1976,0)
          CALL IFLUSHFDBSUBS()
          CALL GSTATS(1976,1)
          WRITE(IU06,*) ' ' 
          WRITE(IU06,*) '  FDB FLUSHED AT ',  CDTPRO, ' FROM WAMODEL. '
          CALL FLUSH (IU06)
          LLFLUSH=.FALSE.
        ENDIF


!*      OUTPUT FILES AND RECOVERY FILES ARE DISPOSED WHEN
!       TIME REACHES THE DISPOSE DATE OR WHEN THE MODEL
!       HAS BEEN SIGNALLED TO DO SO.
!       -------------------------------------------------

        IF (CDATEE.EQ.CDTPRO .AND. LOUTINT ) THEN
!       MOVE INTEGRATED PARAMETERS OF ENTIRE GRID TO PERMANENT FILES

!         GRIB DATA IF FDB IS NOT USED: 
          IF (GFLAG20 .AND. .NOT. LFDB .AND. IRANK.EQ.1) THEN
            CALL GSFILE (IU06, IU30, 0, CDTPRO, CDATEF, 'MPP', 'S')
          ENDIF

!         PURE BINARY DATA:
          IF (FFLAG20 .AND. IRANK.EQ.1 ) THEN
            CALL GSFILE (IU06, IU20, 0, CDTPRO, CDATEF, 'MAP', 'S')
          ENDIF
        ENDIF

!       SAVE BOUNDARY VALUE FILE.
        IF(CDTBC.EQ.CDTPRO) THEN
          IF (IBOUNC.EQ.1 .AND. IRANK.EQ.1 ) THEN
            DO II=1,GBOUNC
            CALL GSFILE(IU06, IU19(II), 0, CDTBC, CDTBC,                &
     &        CBCPREF(II), 'S')
            IF (CDTBC.LT.CDATEE)                                        &
     &        CALL HEADBC (IPOGBO(II)-IPOGBO(II-1), IDELPRO,            &
     &                     TH(1), FR(1), IU19(II), IU06) 
            ENDDO
          ENDIF
        ENDIF
 

!*      WAM-NEMO COUPLING (WHEN NO atmospheric model !!!!!!)
!       ----------------------------------------------------
        IF (LWNEMOCOU.AND.(.NOT.LWCOU)) THEN
          NEMOWSTEP=NEMOWSTEP+1

          IF (MOD(NEMOWSTEP,NEMOFRCO)==0) THEN
            CALL GSTATS(1432,0)
!$OMP       PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ) 
            DO JKGLO=IJS,IJL,NPROMA
              KIJS=JKGLO
              KIJL=MIN(KIJS+NPROMA-1,IJL)
              IF(LWNEMOCOUSTK) THEN
                NEMOUSTOKES(KIJS:KIJL) = USTOKES(KIJS:KIJL)
                NEMOVSTOKES(KIJS:KIJL) = VSTOKES(KIJS:KIJL)
              ELSE
                NEMOUSTOKES(KIJS:KIJL) = 0.0_NEMODP
                NEMOVSTOKES(KIJS:KIJL) = 0.0_NEMODP
              ENDIF
              IF(LWNEMOCOUSTRN)  NEMOSTRN(KIJS:KIJL) = STRNMS(KIJS:KIJL)
            ENDDO
!$OMP       END PARALLEL DO
            CALL GSTATS(1432,1)

            CALL UPDNEMOFIELDS
            CALL UPDNEMOSTRESS

            DO JSTPNEMO=NEMOCSTEP,NEMOCSTEP+NEMONSTEP-1
               ! Advance the NEMO model 1 time step
#ifdef WITH_NEMO
               CALL NEMOGCMCOUP_STEP( JSTPNEMO, IDATE, ITIME )
#endif
               WRITE(IU06,*)'NEMO TIME IS : ',JSTPNEMO, IDATE, ITIME
            ENDDO
            NEMOCSTEP=NEMOCSTEP+NEMONSTEP

          ENDIF

        ENDIF

!*    BRANCHING BACK TO 1.0 FOR NEXT PROPAGATION STEP.
      ENDDO ADVECTION

      IF (LHOOK) CALL DR_HOOK('WAMODEL',1,ZHOOK_HANDLE)

END SUBROUTINE WAMODEL
