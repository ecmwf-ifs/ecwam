      SUBROUTINE WAMODEL (NADV, LDSTOP, LDWRRE, L1STCALL)

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
!     J. BIDLOT ECMWF MAY 2002 INTRODUCTION OF NWAM_BLKS TO SPLIT THE 
!                              LENGTH OF THE IJ LOOP OR TO BE THE
!                              NUMBER OF THREADS IF OPEMMP IS USED.
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

!     *CALL* *WAMODEL (NADV, LDSTOP, LDWRRE, L1STCALL)*
!        *NADV*      INTEGER   NUMBER OF ADVECTION ITERATIONS
!                              PER CALL OF WAMODEL, OUTPUT PARAMETER.
!        *LDSTOP*    LOGICAL   SET .TRUE. IF STOP SIGNAL RECEIVED.
!        *LDWRRE*    LOGICAL   SET .TRUE. IF RESTART SIGNAL RECEIVED.
!        *L1STCALL*  LOGICAL   SET .FALSE. AFTER FIRST CALL.

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
!       READ AND DELETED IN SUB IMPLSCH (IU17 AND IU18). THE FILE
!       NAMES ARE CREATED IN SUB CREWFN AND AN IMPLICIT OPEN IS USED.


!       ALL PARAMETERS HAVE TO BE THE VALUES GIVEN AT THE END OF THE
!       PREPROC OUTPUT IN COLUMN 'REQUIRED'.

!     EXTERNALS.
!     ----------

!       *ABORT1*    - TERMINATES PROCESSING.
!       *AIRSEA*    - SURFACE LAYER STRESS.
!NEST
!       *BOUINPT*   - BOUNDARY VALUE INPUT.
!NEST
!       *CREWFN*    - CREATES A WIND FILE NAME.
!REFRA
!       *DOTDC*     - READ COMMON REFDOT.
!REFRA
!       *FEMEAN*    - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT.
!       *FILLBL*    - ADD LATITUDES TO A BLOCK.
!       *GSFILE*    - ROUTINE TO DYNAMICALLY FETCH OR DISPOSE FILES.
!       *HEADBC*    - WRITE BOUNDARY OUTPUT FILE HEADER.
!       *IMPLSCH*   - IMPLICIT SCHEME FOR INTEGRATION OF SOURCE
!                     FUNCTIONS IN TIME AND INPUT OF WINDS.
!       *INCDATE*   - UPDATE DATE TIME GROUP.
!REFRA
!       *INTPOL*    - MAP SPECTRUM FROM SIGMA TO OMEGA SPACE.
!REFRA
!NEST
!       *INTSPEC*   - INTERPOLATION OF SPECTRA.
!NEST
!       *MAKEGRID*  - MAKES GRIDDED FIELDS.
!       *MPEXCHNG*  - PERFORMED MESSAGE PASSING OF THE SPECTRUM
!       *OUTBC*     - OUTPUT OF BOUNDARY VALUES.
!       *OUTBS*     - CONTROLS OUTPUT FROM BLOCKS.
!       *OUTGRID*   - SAVE BLOCKED PARAMETERS INTO GRID ARRAYS.
!       *OUTINT*    - OUTPUT OF INTEGRATED PARAMETERS.
!       *OUTPP*     - WRITES GRIDDED FIELDS TO PRINTER.
!       *OUTSPP*    - OUTPUT OF SPECTRA AT SELECTED POINTS.
!       *PEAKFR*    - COMUTE PEAK FREQUENCY.
!       *PROPAG*    - PROPAGATION SCHEME.
!       *PRSPP*     - PRINT A SPECTRUM.
!       *PRSPPS*    - PRINT A SWELL SPECTRUM.
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
!       *SINPUT*    - COMPUTATION OF INPUT SOURCE FUNCTION, AND
!                     LINEAR CONTRIBUTION OF INPUT SOURCE FUNCTION
!                     TO FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
!       *SNONLIN*   - COMPUTATION OF NONLINEAR TRANSFER RATE AND
!                     DIAGONAL LINEAR CONTRIBUTION OF NONLINEAR SOURCE
!                     FUNCTION TO FUNCTIONAL MATRIX.
!       *SPLITBL*   - STORING OVERLAPPING LATITUDES OF BLOCKS.
!       *STHQ*      - COMPUTATION OF MEAN WAVE DIRECTION AT EACH
!                     GRID POINT.
!       *STRESSO*   - COMPUTATION OF WAVE STRESS.
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

      USE YOWCPBO  , ONLY : IBOUNC   ,NBOUNC    ,GBOUNC  , IPOGBO  ,
     &            CBCPREF
      USE YOWCOUP  , ONLY : LWCOU    ,KCOUSTEP
      USE YOWCURR  , ONLY : LLCHKCFL ,LLCHKCFLA
      USE YOWCOUT  , ONLY : COUTT    ,COUTS    ,FFLAG20  ,FFLAG21  ,
     &            FFLAG25  ,FFLAG26  ,PFLAG20  ,PFLAG21  ,PFLAG25  ,
     &            PFLAG26  ,GFLAG20  ,GFLAG21  ,GFLAG25  ,GFLAG26  ,
     &            JPPFLAG  ,FFLAG    ,PFLAG    ,GFLAG    ,
     &            LFDB     ,NOUTT    ,NOUTS    ,
     &            CASS     ,NASS     ,LOUTINT  ,LOUTSPP  ,
     &            LGRIBAPI_I,LGRIBAPI_O,LRSTST0
      USE YOWCURR  , ONLY : CDTCUR
      USE YOWFPBO  , ONLY : IBOUNF
      USE YOWREFD  , ONLY : THDD     ,THDC     ,SDOT
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3   ,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK0    ,KSEK1    ,
     &            KSEK2    ,KSEK3    ,KSEK4    ,PSEK2    ,PSEK3    ,
     &            PPMISS   ,PPEPS    ,PPREC    ,NTENCODE ,NGRBRESS ,
     &            HOPERS   ,PPRESOL  ,IFS_KSEC1,LGRHDIFS ,LNEWLVTP
      USE YOWGRID  , ONLY : NLONRGG  ,IGL      ,IJS      ,IJL2     ,
     &            IJL      ,IJLT
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE ,CICOVER  ,CIWA
      USE YOWINTP  , ONLY : WHGTTG   ,WDIRTG   ,WPKFTG   ,WMNFTG   ,
     &            USTARG   ,UDIRG    ,TAUWG    ,CDG      ,SMEANG   ,
     &            U10G     ,WHGTAG   ,CWHTAG   ,RANGAG   ,MWP1G    ,
     &            MWP2G    ,WSPRDG   ,C4G      ,BFG      ,QPG      ,
     &            DEPTHG   ,HMAXG    ,TMAXG    ,USTOKESG ,VSTOKESG ,
     &            UCURG    ,VCURG    ,PHIEPSG  ,PHIAWG   ,TAUOCG
      USE YOWINTS  , ONLY : WHGTSG   ,WDIRSG   ,WMNFSG   ,WHGTWG   ,
     &            WDIRWG   ,WMNFWG
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IRGG     ,AMOWEP   ,
     &            AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO
      USE YOWMESPAS, ONLY : LMESSPASS,LFDBIOOUT,LGRIBOUT ,LNOCDIN  ,
     &            LWAVEWIND
      USE YOWMLTSK , ONLY : NCHUNK
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,
     &            KTAG
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,
     &            NBLO     ,NIBLO    ,CLDOMAIN
      USE YOWPCONS , ONLY : ZMISS    ,DEG      ,EPSMIN
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTH    ,DEPTHA   ,DEPTHD   ,
     &            INDEP
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDTPRO   ,CDTRES   ,
     &            CDATER   ,CDATES   ,CDTINTT  ,CDTINTS  ,CDTSPT   ,
     &            CDTSPS   ,CFDB2DSP ,NWFDBREF ,IDELPRO  ,IDELT    ,
     &            IDELWI   ,IREST    ,IDELRES  ,IDELINT  ,IDELINS  ,
     &            IDELSPT  ,IDELSPS  ,ISHALLO  ,IREFRA   ,IASSI    ,
     &            CDTBC    ,IDELBC   ,
     &            IPROPAGS ,
     &            NTASKS   ,NSIZE    ,NENSFNB  ,NTOTENS  ,NSYSNB   ,
     &            NMETNB   ,CDATEA   ,MARSTYPE ,YCLASS   ,YEXPVER  ,
     &            LLSOURCE ,
     &            LANAONLY ,LFRSTFLD ,NWAM_BLKS,IREFDATE ,LFDBOPEN
      USE YOWSPEC, ONLY   : NSTART   ,NEND     ,
     &            U10NEW   ,U10OLD   ,THWNEW   ,THWOLD   ,USNEW    ,
     &            USOLD    ,Z0NEW    ,Z0OLD    ,TAUW     ,FCONST   ,
     &            ROAIRN   ,ROAIRO   ,ZIDLNEW  ,ZIDLOLD  ,
     &            FL1      ,FL2      ,FL3      ,SL
!!                  SL IS used as a FL array when intpol is called !
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH    ,CWI      ,LRESTARTED
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,SUMWN    ,
     &            WLATN    ,WLONN    ,WCORN    ,WKPMN    ,WMPMN    ,
     &            LLWLATN  ,LLWLONN  ,LLWCORN  ,LLWKPMN  ,LLWMPMN  ,
     &            LUPDTWGHT 
      USE YOWUNIT  , ONLY : IU02     ,IU08     ,IU11     ,IU12     ,
     &            IU13     ,IU14     ,IU19     ,IU20     ,IU21     ,
     &            IU25     ,IU26     ,IU04     ,IU30
      USE YOWWAMI  , ONLY : CBPLTDT  ,CEPLTDT  ,IANALPD  ,IFOREPD  ,
     &            IDELWIN  ,NFCST    ,ISTAT
      USE YOWWIND  , ONLY : CDATEWO
      USE MPL_MODULE
      USE FDBSUBS_MOD
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      INTEGER :: IFS_KSEC1_BAK(NKSEK1)
      CHARACTER*2 MARSTYPEBAK
      CHARACTER*14 CDATEWH, ZERO
      CHARACTER*14 CDTINTTBAK,CDTINTSBAK,CDTSPTBAK,CDTSPSBAK 
      CHARACTER*14 CDATE, CDTPRA, CDTIMP, CDTIMPNEXT
      LOGICAL LDSTOP, LDWRRE, LLFLUSH
      LOGICAL :: FFLAGBAK(JPPFLAG), PFLAGBAK(JPPFLAG), GFLAGBAK(JPPFLAG)
      LOGICAL :: LSV, LRST, LOUT,L1STCALL,LLINTPOL,L4VOUT 
      LOGICAL,ALLOCATABLE,DIMENSION(:) :: LCFLFAIL

!     DEFINE LOCAL ALLOCATABLE ARRAYS

      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: FL4
      REAL(KIND=JPRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

!*    *ADDITIONAL PARAMETERS FOR MULTITASKING

      INTEGER ITSKCA(9,NTASKS)
      EXTERNAL PHYS

! Mod for OPENMP
      INTEGER JKGLO,KIJS,KIJL,NPROMA
! End Mod for OPENMP

      LOGICAL NEWREAD, NEWFILE
      DATA NEWFILE / .FALSE. /
      ZERO = ' '
      IF (LHOOK) CALL DR_HOOK('WAMODEL',0,ZHOOK_HANDLE)

! ----------------------------------------------------------------------

!*    0. WRITE OUT NORMS FOR WAVE ENERGY, PEAK PERIOD AND WIND SPEED
!        -----------------------------------------------------------
      IF (ITEST .GE. 2) THEN
        WRITE(IU06,*)
        WRITE(IU06,*) '   SUB. WAMODEL: JUST ENTERED FOR CDTPRO= ',
     &   CDTPRO, ' AND CDATEA= ', CDATEA
        CALL FLUSH(IU06)
      ENDIF
      LLFLUSH = .TRUE.

      IG=1

      IF(CDTPRO.EQ.CDATEA .AND. LLSOURCE .AND. LICERUN ) THEN
!       INSURE THERE IS SOME WAVE ENERGY FOR GRID POINTS THAT HAVE BEEN
!       FREED FROM SEA ICE (ONLY DONE INITIALLY AND IF THE MODEL IS NOT
!       RESTARTED).
!       IT ALSO RESETS THE MIMIMUM ENERGY LEVEL THAT MIGHT HAVE BEEN LOST
!       WHEN GETTING THE DATA FROM GRIB.
        CALL GSTATS(1236,0)
        NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=IJS(IG),IJL(IG),NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
          CALL UNSETICE(FL1(KIJS:KIJL,1:NANG,1:NFRE),KIJS,KIJL,
     &                  U10OLD(KIJS:KIJL,1),THWOLD(KIJS:KIJL,1))
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1236,1)
      ENDIF

      KTAG=200

      LRSTST0=.FALSE.

!     TIME FOR THE NEXT SOURCE TERM INTEGRATION
      CDTIMPNEXT=CDTPRO
      CALL INCDATE(CDTIMPNEXT,IDELT)   
!     TIME FOR WIND INPUT UPDATE (SEE NEWWIND)
      CDTIMP=CDTPRO

      IF(CDTPRO.EQ.CDATEA .OR. CDTPRO.EQ.CDATEF) THEN

        CDTINTTBAK=CDTINTT
        CDTINTSBAK=CDTINTS
        CDTSPTBAK=CDTSPT
        CDTSPSBAK=CDTSPS
        CDTINTT=CDTPRO
        CDTINTS=CDTPRO
        CDTSPT=CDTPRO
        CDTSPS=CDTPRO


        DO IJ=IJS(IG),IJL(IG)
          U10NEW(IJ) = U10OLD(IJ,IG)
          THWNEW(IJ) = THWOLD(IJ,IG)
          USNEW(IJ) = USOLD(IJ,IG)
          ROAIRN(IJ) = ROAIRO(IJ,IG)
          ZIDLNEW(IJ) = ZIDLOLD(IJ,IG)
        ENDDO

        IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
!         convert output spectra to absolute frame of reference
          CALL GSTATS(1236,0)
          IRA=1
          NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1
!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
          DO JKGLO=IJS(IG),IJL(IG),NPROMA
             KIJS=JKGLO
             KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
             CALL INTPOL (FL1, SL, KIJS, KIJL, IG, IRA)
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1236,1)
          LLINTPOL=.TRUE.

          CALL OUTBS (SL, FL3, IJS(IG), IJL(IG), IG, IGL,
     &                IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,
     &                TAUW, NSTART,NEND)
        ELSE
          NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1
!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1)
!$OMP+    PRIVATE(JKGLO,KIJS,KIJL,K,M,IJ)
          DO JKGLO=IJS(IG),IJL(IG),NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
            DO M=1,NFRE
              DO K=1,NANG
                DO IJ=KIJS,KIJL
                  SL(IJ,K,M) = FL1(IJ,K,M)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO

          CALL OUTBS (SL, FL3, IJS(IG), IJL(IG), IG, IGL,
     &              IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,
     &              TAUW, NSTART,NEND)
        ENDIF

!       SET FL1 ON ICE POINTS TO ZERO

          IF (LICERUN .AND. LMASKICE) THEN
            IF (ITEST.GE.1) THEN
              WRITE(IU06,*) '   SUB. WAMODEL: INITIAL SPECTRUM = 0 AT',
     &         ' ICEPOINTS'
            ENDIF
! Mod for OPENMP
! If NWAM_BLKS=1 then the code will be as before
            NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1
            CALL GSTATS(1439,0)
!$OMP       PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
            DO JKGLO=IJS(IG),IJL(IG),NPROMA
              KIJS=JKGLO
              KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
              CALL SETICE(FL1(KIJS:KIJL,1:NANG,1:NFRE),KIJS,KIJL,IG)
            ENDDO
!$OMP       END PARALLEL DO
            CALL GSTATS(1439,1)
          ENDIF

!       PRINT OUT NORMS

        CALL OUTWNORM
        CALL FLUSH(IU06)

!*    0.1 SAVE INITIAL INTEGRATED FIELDS
!         ------------------------------

        IF(MARSTYPE.EQ.'cf'.OR.MARSTYPE.EQ.'pf'.OR.MARSTYPE.EQ.'fc'.OR.
     &     MARSTYPE.EQ.'4v'.OR.
     &     LANAONLY.OR.LFRSTFLD) THEN

          IF(LFRSTFLD) THEN
            FFLAGBAK=FFLAG
            PFLAGBAK=PFLAG
            GFLAGBAK=GFLAG
            MARSTYPEBAK=MARSTYPE
            IF(MARSTYPE.EQ.'fg'.OR.MARSTYPE.EQ.'an') THEN
              IF(.NOT.LNOCDIN) THEN
                FFLAG(7)=.FALSE.
                PFLAG(7)=.FALSE.
                GFLAG(7)=.FALSE.
              ENDIF
              IF(LWAVEWIND) THEN
                FFLAG(10)=.FALSE.
                PFLAG(10)=.FALSE.
                GFLAG(10)=.FALSE.
              ENDIF
              GFLAG(22)=.FALSE.
              GFLAG(23)=.FALSE.
              GFLAG(24)=.FALSE.
              MARSTYPE='an'
            ENDIF
          ENDIF

!         NEED TO TEMPORARLY RESET THE IFS FORECAST STEP
!         IF THE IFS GRIB HEADER IS USED SUCH THAT IT POINTS TO
!         THE START OF THE RUN.
          IF(LGRHDIFS) THEN
             LRSTST0=.TRUE.

             IF(.NOT.LGRIBAPI_I) THEN 
               IFS_KSEC1_BAK(16)=IFS_KSEC1(16)  ! value to reset after output
               IF(IFS_KSEC1(16).GT.0) THEN
                 IF(IFS_KSEC1(15).EQ.0) THEN
                   IFS_KSEC1(16)=IFS_KSEC1(16)-(KCOUSTEP/60) 
                 ELSEIF(IFS_KSEC1(15).EQ.1) THEN
                   IFS_KSEC1(16)=IFS_KSEC1(16)-(KCOUSTEP/3600) 
                 ELSEIF(IFS_KSEC1(15).EQ.10) THEN
                   IFS_KSEC1(16)=IFS_KSEC1(16)-(KCOUSTEP/(3*3600)) 
                 ELSEIF(IFS_KSEC1(15).EQ.11) THEN
                   IFS_KSEC1(16)=IFS_KSEC1(16)-(KCOUSTEP/(6*3600)) 
                 ELSEIF(IFS_KSEC1(15).EQ.12) THEN
                   IFS_KSEC1(16)=IFS_KSEC1(16)-(KCOUSTEP/(12*3600)) 
                 ELSEIF(IFS_KSEC1(15).EQ.2) THEN
                   IFS_KSEC1(16)=IFS_KSEC1(16)-(KCOUSTEP/(24*3600)) 
                 ELSEIF(IFS_KSEC1(15).EQ.254) THEN
                   IFS_KSEC1(16)=IFS_KSEC1(16)-KCOUSTEP 
                 ENDIF
               ENDIF
!              IN VAREPS THE FC STEP 0 IS ALWAYS AN OVERLAP FOR LEG 2
!              AND ABOVE
               IFS_KSEC1_BAK(40)=IFS_KSEC1(40)  ! value to reset after output
               IINDAT=IFS_KSEC1(46)
               IHHMM=IFS_KSEC1(47)
               NLEG=IFS_KSEC1(48)
               IF(NLEG.GE.2 .AND. IINDAT.NE.0 .AND. IHHMM.NE.0 ) THEN
                 IFS_KSEC1(40)=1034
               ENDIF
             ENDIF 

          ENDIF 

          IF( .NOT. LRESTARTED ) THEN
            CALL OUTINT (IU06,IU20,IU21,IU30,ITEST,NWFDBREF,LFDBOPEN)
          ENDIF

          IF (NWFDBREF.NE.-5.AND.LFDB) THEN
             CALL GSTATS(1976,0)
             ISTAT = IFLUSHFDBSUBS (NWFDBREF)
!J           CALL MPL_BARRIER(CDSTRING='WAMODEL:')
             CALL GSTATS(1976,1)
          ELSE
             CALL GSTATS(753,0)
!J           CALL MPL_BARRIER(CDSTRING='WAMODEL:')
             CALL GSTATS(753,1)
          ENDIF

          IF(LFRSTFLD) THEN
            MARSTYPE=MARSTYPEBAK
            FFLAG=FFLAGBAK
            PFLAG=PFLAGBAK
            GFLAG=GFLAGBAK
            LFRSTFLD=.FALSE.
          ENDIF

          CDTINTT=CDTINTTBAK
          CDTINTS=CDTINTSBAK

          IF(LGRHDIFS) THEN
            LRSTST0=.FALSE.
            IF(.NOT.LGRIBAPI_I) THEN 
              IFS_KSEC1(16)=IFS_KSEC1_BAK(16)
              IFS_KSEC1(40)=IFS_KSEC1_BAK(40)
            ENDIF
          ENDIF

          IF (ITEST .GE. 2) THEN
            WRITE(IU06,*)
            WRITE(IU06,*) '   SUB. WAMODEL: INITIAL FIELDS SAVED '
            WRITE(IU06,*) '                 FOR FORECAST RUN     '
            WRITE(IU06,*) '    '
            CALL FLUSH(IU06)
          ENDIF
      
        ENDIF

        IF(.NOT.LANAONLY.AND.ALLOCATED(WHGTTG)) DEALLOCATE(WHGTTG)
        IF(ALLOCATED(WDIRTG)) DEALLOCATE(WDIRTG)
        IF(ALLOCATED(WPKFTG)) DEALLOCATE(WPKFTG)
        IF(ALLOCATED(WMNFTG)) DEALLOCATE(WMNFTG)
        IF(ALLOCATED(USTARG)) DEALLOCATE(USTARG)
        IF(ALLOCATED(UDIRG)) DEALLOCATE(UDIRG)
        IF(ALLOCATED(TAUWG)) DEALLOCATE(TAUWG)
        IF(ALLOCATED(CDG)) DEALLOCATE(CDG)
        IF(ALLOCATED(SMEANG)) DEALLOCATE(SMEANG)
        IF(ALLOCATED(U10G)) DEALLOCATE(U10G)
        IF(ALLOCATED(MWP1G)) DEALLOCATE(MWP1G)
        IF(ALLOCATED(MWP2G)) DEALLOCATE(MWP2G)
        IF(ALLOCATED(WSPRDG)) DEALLOCATE(WSPRDG)
        IF(ALLOCATED(WHGTAG)) DEALLOCATE(WHGTAG)
        IF(ALLOCATED(CWHTAG)) DEALLOCATE(CWHTAG)
        IF(ALLOCATED(RANGAG)) DEALLOCATE(RANGAG)
        IF(ALLOCATED(C4G)) DEALLOCATE(C4G)
        IF(ALLOCATED(BFG)) DEALLOCATE(BFG)
        IF(ALLOCATED(QPG)) DEALLOCATE(QPG)
        IF(ALLOCATED(DEPTHG)) DEALLOCATE(DEPTHG)
        IF(ALLOCATED(HMAXG)) DEALLOCATE(HMAXG)
        IF(ALLOCATED(TMAXG)) DEALLOCATE(TMAXG)
        IF(ALLOCATED(USTOKESG)) DEALLOCATE(USTOKESG)
        IF(ALLOCATED(VSTOKESG)) DEALLOCATE(VSTOKESG)
        IF(ALLOCATED(UCURG)) DEALLOCATE(UCURG)
        IF(ALLOCATED(VCURG)) DEALLOCATE(VCURG)
        IF(ALLOCATED(PHIEPSG)) DEALLOCATE(PHIEPSG)
        IF(ALLOCATED(PHIAWG)) DEALLOCATE(PHIAWG)
        IF(ALLOCATED(TAUOCG)) DEALLOCATE(TAUOCG)

        IF(ALLOCATED(WHGTSG)) DEALLOCATE(WHGTSG)
        IF(ALLOCATED(WDIRSG)) DEALLOCATE(WDIRSG)
        IF(ALLOCATED(WMNFSG)) DEALLOCATE(WMNFSG)
        IF(ALLOCATED(WHGTWG)) DEALLOCATE(WHGTWG)
        IF(ALLOCATED(WDIRWG)) DEALLOCATE(WDIRWG)
        IF(ALLOCATED(WMNFWG)) DEALLOCATE(WMNFWG)

      ENDIF

      IF(IGL.GT.1) THEN
         WRITE (IU06,*) '**************************************'
         WRITE (IU06,*) '*                                    *'
         WRITE (IU06,*) '*    FATAL ERROR                     *'
         WRITE (IU06,*) '*    ===========                     *'
         WRITE (IU06,*) '* A MULTI BLOCK VERSION WAS REQUESTED*'
         WRITE (IU06,*) '* THIS OPTION IS NOT IMPLEMENTED FOR *'
         WRITE (IU06,*) '*                                    *'
         WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.  *'
         WRITE (IU06,*) '*                                    *'
         WRITE (IU06,*) '**************************************'
         CALL ABORT1
      ENDIF

!*    1. ADVECTION TIME LOOP.
!        --------------------

      ADVECTION : DO KADV = 1,NADV

        IF(.NOT.ALLOCATED(FL4).AND..NOT.LMESSPASS)
     &   ALLOCATE(FL4(NINF-1:NSUP,NANG,NFRE))

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
          CDTINTT = ZERO
          CDTINTS = ZERO
          CDTSPT  = ZERO
          CDTSPS  = ZERO
          DO J=1,NOUTT
            IF (CDTPRO.EQ.COUTT(J)) THEN
              IF (PFLAG20.OR.FFLAG20.OR.GFLAG20) CDTINTT = COUTT(J)
              IF (PFLAG21.OR.FFLAG21.OR.GFLAG21) CDTINTS = COUTT(J)
              IF (PFLAG25.OR.FFLAG25.OR.GFLAG25) CDTSPT  = COUTT(J)
              IF (PFLAG26.OR.FFLAG26.OR.GFLAG26) CDTSPS  = COUTT(J)
            ENDIF
          ENDDO
        ELSE
          IF ((PFLAG20.OR.FFLAG20.OR.GFLAG20) .AND.
     &         CDTINTT.LT.CDTPRO)
     &     CALL INCDATE (CDTINTT,IDELINT)
          IF ((PFLAG21.OR.FFLAG21.OR.GFLAG21) .AND.
     &         CDTINTS.LT.CDTPRO)
     &     CALL INCDATE (CDTINTS,IDELINS)
          IF ((PFLAG25.OR.FFLAG25.OR.GFLAG25) .AND.
     &         CDTSPT.LT.CDTPRO)
     &     CALL INCDATE (CDTSPT,IDELSPT)
          IF ((PFLAG26.OR.FFLAG26.OR.GFLAG26) .AND.
     &         CDTSPS.LT.CDTPRO)
     &     CALL INCDATE (CDTSPS,IDELSPS)
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

        IF ((IBOUNC.EQ.1.OR.IBOUNF.EQ.1).AND.CDTBC.LT.CDTPRO)
     &  CALL INCDATE(CDTBC,IDELBC)


!*    1.2 RESET FILES.
!         ------------
        IF(.NOT.LMESSPASS) REWIND IU12

        IF (IGL.NE.1) THEN
          REWIND IU08(IPROPAGS)
          REWIND IU11
          REWIND IU13
          REWIND IU14

!*    1.3 CHANGE WORKING FILE UNITS.
!         --------------------------

          IFIL=IU11
          IU11=IU12
          IU12=IFIL
          IFIL=IU13
          IU13=IU14
          IU14=IFIL

!*    1.4 READ FIRST BLOCK OF SPECTRA, IF MULTI BLOCK VERSION.
!         ----------------------------------------------------

          READ(IU11) FL1
          IF (ITEST.GE.2)
     &     WRITE(IU06,*) '   SUB. WAMODEL: START OF READING',
     &     ' FIRST BLOCK OF SPECTRA'
        ENDIF

!*    1.5  LOOP FOR BLOCKS OF LATITUDES.
!          -----------------------------

        BLOCK : DO IG=1,IGL


!*    1.5.1 EVEN OR ODD BLOCK NUMBER
!          -------------------------

          IGH1=(IG/2)*2


!*    1.5.2 INPUT NEIGHBOURING GRID POINT INDICES (COMMON BLOCK UBUF).
!           ----------------------------------------------------------

          IF (IGL.NE.1) THEN
            READ (IU08(IPROPAGS)) KLAT
            READ (IU08(IPROPAGS)) KLON
!SHALLOW

!*    1.5.3 COMPUTE SHALLOW WATER TABLE INDICES.
!           ------------------------------------

            IF (ISHALLO.NE.1) THEN
              DO IJ=NINF,IJLT(IG)
                XD = LOG(DEPTH(IJ,IG)/DEPTHA)/LOG(DEPTHD)+1.
                ID = NINT(XD)
                ID = MAX(ID,1)
                INDEP(IJ) = MIN(ID,NDEPTH)
              ENDDO
            ENDIF
!SHALLOW
            IF (ITEST.GE.2) THEN
              IF (ITESTB.GE.IG) THEN
                WRITE(IU06,*) '   SUB. WAMODEL: COMMON UBUF ',
     &           'READ AND INITIALIZED FOR BLOCK IG = ', IG
              ENDIF
            ENDIF

!*    1.5.4 COUPLING WITH NEIGHBOURING BLOCKS IG +- 1 AND START
!*          INPUT OF SPECTRA FOR BLOCK IG+1.
!           ----------------------------------------------------

            IF(IGH1.EQ.IG) THEN
              CALL FILLBL(FL2, IG, IU13)
              IF (IG.NE.IGL)
     &         READ(IU11) FL1
            ELSE
              CALL FILLBL(FL1, IG, IU13)
              IF (IG.NE.IGL)
     &         READ(IU11) FL2
            ENDIF
            IF (ITEST.GE.2) THEN
              IF (ITESTB.GE.IG) THEN
                WRITE(IU06,*) '   SUB. WAMODEL: LATITUDES ADDED',
     &           ' TO BLOCK AND START OF READING',
     &           ' FOR NEXT BLOCK'
              ENDIF
            ENDIF
          ENDIF


!*    1.5.5 SET TIME COUNTER.
!           -----------------
          CDATE   = CDTPRA
          CDATEWH = CDATEWO    
          NEWREAD = .FALSE.

 1550     CONTINUE                       

!*    1.5.5.2 COMPUTATION OF PROPAGATION &&& 
!*            INTEGRATION OF SOURCE TERMS OVER SUB TIME STEPS BETWEEN
!*            PROPAGATION TIME STEPS.
!             -------------------------------------------------------

! IF CDATE CORRESPONDS TO A PROPAGATION TIME, THEN CALL MPEXCHNG
! TO OBTAIN INFORMATION AT NEIGHBORING GRID POINTS LOCATED ON ADJACENT
! PROCESSES.
          IF(CDATE.EQ.CDTPRA) THEN
            CALL MPEXCHNG(FL1,NANG,NFRE)

            IF (ITEST.GE.2) THEN
              WRITE(IU06,*) '   SUB. WAMODEL: MPEXCHNG CALLED' 
              CALL FLUSH (IU06)
            ENDIF
          ENDIF

! Mod for OPENMP
! If NWAM_BLKS=1 then the code will be as before

          NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1

          IF (CDATE.EQ.CDTPRA) THEN

!           IF CTU SCHEME IS USED, COMPUTE THE WEIGHTS
            IF(IPROPAGS.EQ.2. AND. LUPDTWGHT) THEN
              ALLOCATE(LCFLFAIL(IJS(IG):IJL(IG)))

              IF(.NOT. ALLOCATED(SUMWN))
     &           ALLOCATE(SUMWN(IJS(IG):IJL(IG),NANG,NFRE))
              IF(.NOT. ALLOCATED(WLATN))
     &           ALLOCATE(WLATN(IJS(IG):IJL(IG),NANG,NFRE,2,2))
              IF(.NOT. ALLOCATED(LLWLATN))
     &           ALLOCATE(LLWLATN(NANG,NFRE,2,2))

              IF(.NOT. ALLOCATED(WLONN))
     &           ALLOCATE(WLONN(IJS(IG):IJL(IG),NANG,NFRE,2))
              IF(.NOT. ALLOCATED(LLWLONN))
     &           ALLOCATE(LLWLONN(NANG,NFRE,2))

              IF(.NOT. ALLOCATED(WCORN))
     &           ALLOCATE(WCORN(IJS(IG):IJL(IG),NANG,NFRE,4,2))
              IF(.NOT. ALLOCATED(LLWCORN))
     &           ALLOCATE(LLWCORN(NANG,NFRE,4,2))

              IF(.NOT. ALLOCATED(WKPMN))
     &           ALLOCATE(WKPMN(IJS(IG):IJL(IG),NANG,NFRE,-1:1))
              IF(.NOT. ALLOCATED(LLWKPMN))
     &           ALLOCATE(LLWKPMN(NANG,NFRE,-1:1))

              IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
                IF(.NOT. ALLOCATED(WMPMN))
     &             ALLOCATE(WMPMN(IJS(IG):IJL(IG),NANG,NFRE,-1:1))
                IF(.NOT. ALLOCATED(LLWMPMN))
     &             ALLOCATE(LLWMPMN(NANG,NFRE,-1:1))
              ENDIF


              CALL GSTATS(1430,0)
!$OMP         PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
              DO JKGLO=IJS(IG),IJL(IG),NPROMA
                KIJS=JKGLO
                KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
                CALL CTUW(KIJS,KIJL,IG,LCFLFAIL(JKGLO))
              ENDDO
!$OMP         END PARALLEL DO
              CALL GSTATS(1430,1)
              IF (ITEST.GE.2) THEN
                WRITE(IU06,*) '   SUB. WAMODEL: CTUW CALLED'
                CALL FLUSH (IU06)
              ENDIF

              DO IJ=IJS(IG),IJL(IG)
                IF(LCFLFAIL(IJ)) THEN
                  WRITE(0,*) '!!! ********************************* !!'
                  WRITE(0,*) '!!! WAVE MODEL HAS ABORTED !!!'
                  WRITE(0,*) '!!! FOLLOWING CFL CRITERION VIOLATION !!'
                  WRITE(0,*) '!!! ON PE ',IRANK
                  WRITE(0,*) '!!! ********************************* !!'
                  CALL ABORT1
                ENDIF
              ENDDO

!             FIND THE LOGICAL FLAGS THAT WILL LIMIT THE EXTEND OF THE CALCULATION
!             IN PROPAGS2

              DO IC=1,2
                DO ICL=1,2
                  DO K=1,NANG
                    DO M=1,NFRE
                      LLWLATN(K,M,IC,ICL)=.FALSE.
                      DO IJ=IJS(IG),IJL(IG)
                        IF(WLATN(IJ,K,M,IC,ICL).GT.0.) THEN
                          LLWLATN(K,M,IC,ICL)=.TRUE.
                          EXIT
                        ENDIF
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

              DO IC=1,2
                DO M=1,NFRE
                  DO K=1,NANG
                    LLWLONN(K,M,IC)=.FALSE.
                    DO IJ=IJS(IG),IJL(IG)
                      IF(WLONN(IJ,K,M,IC).GT.0.) THEN
                        LLWLONN(K,M,IC)=.TRUE.
                        EXIT
                      ENDIF
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

              DO ICL=1,2
                DO ICR=1,4
                  DO M=1,NFRE
                    DO K=1,NANG
                      LLWCORN(K,M,ICR,ICL)=.FALSE.
                      DO IJ=IJS(IG),IJL(IG)
                        IF(WCORN(IJ,K,M,ICR,ICL).GT.0.) THEN
                          LLWCORN(K,M,ICR,ICL)=.TRUE.
                          EXIT
                        ENDIF
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

              DO IC=-1,1
                DO M=1,NFRE
                  DO K=1,NANG
                   LLWKPMN(K,M,IC)=.FALSE.
                    DO IJ=IJS(IG),IJL(IG)
                      IF(WKPMN(IJ,K,M,IC).GT.0.) THEN
                        LLWKPMN(K,M,IC)=.TRUE.
                        EXIT
                      ENDIF
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

              IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
                DO IC=-1,1
                  DO M=1,NFRE
                    DO K=1,NANG
                     LLWMPMN(K,M,IC)=.FALSE.
                      DO IJ=IJS(IG),IJL(IG)
                        IF(WMPMN(IJ,K,M,IC).GT.0.) THEN
                          LLWMPMN(K,M,IC)=.TRUE.
                          EXIT
                        ENDIF
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF


              LUPDTWGHT=.FALSE.
              DEALLOCATE(LCFLFAIL)
              IF(ALLOCATED(THDD)) DEALLOCATE(THDD)
              IF(ALLOCATED(THDC)) DEALLOCATE(THDC)
              IF(ALLOCATED(SDOT)) DEALLOCATE(SDOT)
            ENDIF

            IF(L1STCALL .OR. LLCHKCFLA) LLCHKCFL=.TRUE.


            IF(NIBLO.GT.1) THEN
              IF(IPROPAGS.EQ.1) THEN
                CALL GSTATS(1430,0)
!$OMP           PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
                DO JKGLO=IJS(IG),IJL(IG),NPROMA
                  KIJS=JKGLO
                  KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
                  CALL PROPAGS1(FL1,FL3,KIJS,KIJL,IG,FCONST,L1STCALL)
                ENDDO
!$OMP           END PARALLEL DO
                CALL GSTATS(1430,1)
                IF (ITEST.GE.2) THEN
                  WRITE(IU06,*) '   SUB. WAMODEL: PROPAGS1 CALLED'
                  CALL FLUSH (IU06)
                ENDIF
              ELSEIF(IPROPAGS.EQ.2) THEN
                CALL GSTATS(1430,0)
!$OMP           PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
                DO JKGLO=IJS(IG),IJL(IG),NPROMA
                  KIJS=JKGLO
                  KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
                  CALL PROPAGS2(FL1,FL3,KIJS,KIJL,IG)
                ENDDO
!$OMP           END PARALLEL DO
                CALL GSTATS(1430,1)
                IF (ITEST.GE.2) THEN
                  WRITE(IU06,*) '   SUB. WAMODEL: PROPAGS2 CALLED'
                  CALL FLUSH (IU06)
                ENDIF
              ELSE
                CALL GSTATS(1430,0)
!$OMP           PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
                DO JKGLO=IJS(IG),IJL(IG),NPROMA
                  KIJS=JKGLO
                  KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
                  CALL PROPAGS(FL1,FL3,KIJS,KIJL,IG,FCONST,L1STCALL)
                ENDDO
!$OMP           END PARALLEL DO
                CALL GSTATS(1430,1)
                IF (ITEST.GE.2) THEN
                  WRITE(IU06,*) '   SUB. WAMODEL: PROPAGS CALLED'
                  CALL FLUSH (IU06)
                ENDIF
              ENDIF

            ELSE
              FL3 = FL1
            ENDIF

            L1STCALL=.FALSE.
            LLCHKCFL=.FALSE.

            CDATE=CDTPRO   
          ENDIF

!         IT IS TIME TO INTEGRATE THE SOURCE TERMS
          IF(CDATE.GE.CDTIMPNEXT) THEN

!*          READING OF NEW WINDS IF NEEDED.
!           -------------------------------

            CALL NEWWIND(IJS(IG),IJL(IG),IG,IGL,CDTIMP,CDATEWH,
     &                   NEWREAD,NEWFILE,U10OLD,THWOLD,U10NEW,THWNEW,
     &                   ROAIRO, ROAIRN, ZIDLOLD,ZIDLNEW,CICOVER)
            IF (ITEST.GE.2) THEN
              IF (ITESTB.GE.IG) THEN
                WRITE(IU06,*) '   SUB. WAMODEL: NEWWIND CALLED' 
              ENDIF
            ENDIF


            CALL GSTATS(1431,0)
            IF(LLSOURCE) THEN
!$OMP         PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
              DO JKGLO=IJS(IG),IJL(IG),NPROMA
                KIJS=JKGLO
                KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
                CALL IMPLSCH (FL3(KIJS:KIJL,:,:),
     &                      FL1(KIJS:KIJL,:,:),
     &                      KIJS, KIJL, IG,
     &                      THWOLD(KIJS,IG),USOLD(KIJS,IG),
     &                      TAUW(KIJS,IG),Z0OLD(KIJS,IG),
     &                      ROAIRO(KIJS,IG),ZIDLOLD(KIJS,IG),
     &                      CICOVER(KIJS,IG),CIWA(KIJS:KIJL,:,IG),
     &                      U10NEW(KIJS),THWNEW(KIJS),USNEW(KIJS),
     &                      Z0NEW(KIJS),ROAIRN(KIJS),ZIDLNEW(KIJS),
     &                      SL(KIJS:KIJL,:,:),
     &                      FCONST(KIJS:KIJL,:))
              ENDDO
!$OMP       END PARALLEL DO
              IF (ITEST.GE.2) THEN
                WRITE(IU06,*) '   SUB. WAMODEL: IMPLSCH CALLED'
                CALL FLUSH (IU06)
              ENDIF

            ELSE
!$OMP         PARALLEL DO SCHEDULE(STATIC)  
!$OMP+        PRIVATE(JKGLO,KIJS,KIJL,IJ,K,M)
              DO JKGLO=IJS(IG),IJL(IG),NPROMA
                KIJS=JKGLO
                KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
                DO IJ=KIJS,KIJL
                  DO K=1,NANG
                    DO M=1,NFRE
                      FL3(IJ,K,M) = MAX(FL3(IJ,K,M),EPSMIN)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
!$OMP         END PARALLEL DO
            ENDIF
            CALL GSTATS(1431,1)
! End Mod for OPENMP


!*          COPY AND CLOSE FILES IF NEEDED.
!           -------------------------------

            CALL CLOSEND(IJS(IG),IJL(IG),IG,IGL,CDTIMP,CDATEWH,
     &                   NEWREAD,NEWFILE,U10OLD,THWOLD,ROAIRO,ZIDLOLD,
     &                   U10NEW,THWNEW,ROAIRN,ZIDLNEW)
            IF (ITEST.GE.2) THEN
              IF (ITESTB.GE.IG) THEN
                WRITE(IU06,*) '   SUB. WAMODEL: CLOSEND CALLED' 
                CALL FLUSH (IU06)
              ENDIF
            ENDIF

            CDTIMP=CDTIMPNEXT
            CALL INCDATE(CDTIMPNEXT,IDELT)   
          ENDIF

!*    1.5.5.4 UPDATE TIME;IF TIME LEFT BRANCH BACK TO 1.5.5 
!             ---------------------------------------------

          IF (CDTIMPNEXT.LE.CDTPRO) GO TO 1550

!        END OF TIME LOOP ALL TIME STEPS DONE

          IF (ITEST.GE.2) THEN
            IF (ITESTB.GE.IG) THEN
              WRITE(IU06,*) '   SUB. WAMODEL: SOURCE FUNCTIONS',
     &         ' INTEGRATED:'
            ENDIF
          ENDIF
!ICE

!*    1.5.6 SET FL3 ON ICE POINTS TO ZERO
!           -----------------------------

          IF (LICERUN .AND. LMASKICE) THEN
            IF (ITEST.GE.1) THEN
              WRITE(IU06,*) '   SUB. WAMODEL: SPECTRUM = 0 AT',
     &         ' ICEPOINTS'
            ENDIF
            CALL GSTATS(1439,0)
!$OMP       PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
            DO JKGLO=IJS(IG),IJL(IG),NPROMA
              KIJS=JKGLO
              KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
              CALL SETICE(FL3(KIJS:KIJL,1:NANG,1:NFRE),KIJS,KIJL,IG)
            ENDDO
!$OMP       END PARALLEL DO
            CALL GSTATS(1439,1)
          ENDIF

!ICE
!NEST

!*    1.5.7 INPUT OF BOUNDARY VALUES.
!           -------------------------

          CALL GSTATS(1909,0)
          IF (IBOUNF.EQ.1) THEN
            CALL BOUINPT (FL3, IU02, NSTART, NEND)
            IF (ITEST.GE.2) THEN
              IF (ITESTB.GE.IG) WRITE(IU06,*)
     &         '   SUB. WAMODEL: BOUNDARY VALUES INSERTED'
            ENDIF
          ENDIF
!NEST

!*    1.5.8 OUTPUT OF SECOND LATITUDES FOR NEXT TIMESTEP
!*          START OF OUTPUT OF ALL SPECTRA IN BLOCK.
!           --------------------------------------------

          IF (IGL.NE.1) THEN
            IF (IG.NE.1) CALL SPLITBL(FL3, IJS(IG), IJL2(IG), IU14)
            WRITE(IU12) FL3
            IF (ITEST.GE.2) THEN
              IF (ITESTB.GE.IG) WRITE(IU06,*)
     &         '   SUB. WAMODEL: OVERLAPPING LATITUDES',
     &         ' SAVED AND START OF BLOCK OUTPUT'
            ENDIF
          ENDIF
!NEST

!*    1.5.9 OUTPUT OF BOUNDARY POINTS.
!           --------------------------

          IF (IBOUNC.EQ.1) THEN
            CALL OUTBC (FL3, IJS(IG), IJL(IG), IG, IU19, THWNEW, USNEW)
            IF (ITEST.GE.2) THEN
              IF (ITESTB.GE.IG) WRITE(IU06,*)
     &         '   SUB. WAMODEL: BOUNDARY OUTP',
     &         'UT (COARSE GRID) DONE IN SUB OUTBC'
            ENDIF
          ENDIF
!NEST

!*    1.5.10 MODEL OUTPUT INTEGRATED DATA ARE SAVED IN COMMON.
!*           SPECTRA ARE SAVED ON DISK AND PRINTED.
!            -------------------------------------------------

          IF (ITEST.GE.2) THEN
            IF (ITESTB.GE.IG) THEN
              WRITE(IU06,*) '   SUB. WAMODEL: MODEL ',
     &         'OUTPUT? CDTINTT=',CDTINTT,' CDTINTS=',CDTINTS
              WRITE(IU06,*) '   SUB. WAMODEL: MODEL ',
     &         'OUTPUT? CDTSPT=',CDTSPT,' CDTSPS=',CDTSPS
              WRITE(IU06,*) '   SUB. WAMODEL: MODEL ',
     &         'OUTPUT? IGH1=',IGH1,' IG=',IG
            ENDIF
          ENDIF

          IF(.NOT.LWCOU .AND. .NOT. LDSTOP) THEN
!!!!      the call to CHESIG is a signal handeling facility which is
!!!!      specific to running WAM at ECMWF, it can be ignored when
!!!!      WAM is not run at ECMWF.
            CALL CHESIG (IU06, ITEST, IRANK, NPROC, LDSTOP, LDWRRE)
          ENDIF

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

          LLINTPOL=.FALSE.
          IF (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO .OR.
     &     CDTSPT.EQ.CDTPRO .OR. CDTSPS.EQ.CDTPRO .OR.
     &     CDATEE.EQ.CDTPRO .OR. LRST .OR.
     &     (IASSI.EQ.1 .AND. KADV.EQ.NADV ))
     &    THEN

            IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
!             convert output spectra to absolute frame of reference
              CALL GSTATS(1236,0)
              IRA=1
              NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1
!$OMP         PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
              DO JKGLO=IJS(IG),IJL(IG),NPROMA
                KIJS=JKGLO
                KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
                CALL INTPOL (FL3, SL, KIJS, KIJL, IG, IRA)
              ENDDO
!$OMP         END PARALLEL DO
              CALL GSTATS(1236,1)
              LLINTPOL=.TRUE.

              CALL OUTBS (SL, FL1, IJS(IG), IJL(IG), IG, IGL,
     &                    IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,
     &                    TAUW, NSTART,NEND)
            ELSE

              NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1
!             CALL STATS(1432,0)
!$OMP         PARALLEL DO SCHEDULE(DYNAMIC,1)
!$OMP+        PRIVATE(JKGLO,KIJS,KIJL,K,M,IJ)
              DO JKGLO=IJS(IG),IJL(IG),NPROMA
                KIJS=JKGLO
                KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
                DO M=1,NFRE
                  DO K=1,NANG
                    DO IJ=KIJS,KIJL
                      SL(IJ,K,M) = FL3(IJ,K,M)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
!$OMP         END PARALLEL DO
!             CALL GSTATS(1432,1)

              CALL OUTBS (SL, FL1, IJS(IG), IJL(IG), IG, IGL,
     &                    IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,
     &                    TAUW, NSTART,NEND)
            ENDIF

            IF (ITEST.GE.2) THEN
              IF (ITESTB.GE.IG) WRITE(IU06,*)
     &         '   SUB. WAMODEL: MODEL OUTPUT DONE IN SUB OUTBS'
              CALL FLUSH (IU06)
            ENDIF
          ENDIF

!*    BRANCHING BACK TO 1.5 FOR NEXT BLOCK OF LATITUDES

        ENDDO BLOCK

!*    1.6 IF ONE BLOCK VERSION COPY RESULTS.
!         ----------------------------------
!!!
!!!   note in case of Alt and SAR assimilation both FL1 and FL3
!!!   are needed since FL1 will be first modified by the altimeter
!!!   data assimilation and FL3 is needed to keep the first guess
!!!   (unaltered by the alt data assimilation) for the calculation
!!!   of the difference between SAR data and model values
!!!   (see rearrangsar)
!!!
        IF (IGL.EQ.1) THEN
          IG=1
          NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1
          CALL GSTATS(1432,0)
!$OMP     PARALLEL DO SCHEDULE(STATIC) 
!$OMP+    PRIVATE(JKGLO,KIJS,KIJL,K,M,IJ)
          DO JKGLO=IJS(IG),IJL(IG),NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
            DO M=1,NFRE
              DO K=1,NANG
                DO IJ=KIJS,KIJL
                  FL1(IJ,K,M) = FL3(IJ,K,M)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1432,1)

!         SET THE DUMMY LAND POINTS TO 0.
          DO M=1,NFRE
            DO K=1,NANG
               FL1(NINF-1,K,M) = 0. 
            ENDDO
          ENDDO

          IF (ITEST.GE.2)
     &     WRITE(IU06,*) '   SUB. WAMODEL: BLOCK COPIED'
        ENDIF

!*    1.7 ONE PROPAGATION TIMESTEP DONE FOR ALL BLOCKS.
!         ---------------------------------------------


!         WRITE INTEGRATED DATA TO FILE AND/OR PRINTER
!         DATA WERE COLLECTED INSIDE THE BLOCK LOOP.

        IF (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO) THEN

!       PRINTS OUT NORMS

          CALL OUTWNORM

!         IF THE OUTPUT TIME IS NOT AN ANLYSIS TIME THEN TYPE FG
!         BECOMES TYPE AN (i.e. speudo analysis)
          MARSTYPEBAK=MARSTYPE
          IF((MARSTYPE.EQ.'fg' .OR. MARSTYPE.EQ.'4v') .AND.
     &        KADV.LT.NADV) THEN
            MARSTYPE='an'
          ENDIF
          CALL OUTINT (IU06, IU20, IU21, IU30, ITEST, NWFDBREF,LFDBOPEN)
          MARSTYPE=MARSTYPEBAK

          IF(LMESSPASS) THEN
            CALL GSTATS(753,0)
            CALL MPL_BARRIER(CDSTRING='WAMODEL:')
            CALL GSTATS(753,1)
          ENDIF

        ENDIF

!       DEALLOCATE ARRAYS WHICH WERE ALLOCATED IN OUTBS
        IF (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO .OR.
     &   CDTSPT.EQ.CDTPRO .OR. CDTSPS.EQ.CDTPRO .OR.
     &   CDATEE.EQ.CDTPRO .OR. LRST .OR.
     &   (IASSI.EQ.1 .AND. KADV.EQ.NADV ))
     &  THEN
          IF(IASSI.NE.1) THEN
            IF(ALLOCATED(WHGTTG)) DEALLOCATE(WHGTTG)
          ENDIF

          IF(ALLOCATED(USTARG)) DEALLOCATE(USTARG)
          IF(ALLOCATED(CDG)) DEALLOCATE(CDG)
          IF(ALLOCATED(WDIRTG)) DEALLOCATE(WDIRTG)
          IF(ALLOCATED(WPKFTG)) DEALLOCATE(WPKFTG)
          IF(ALLOCATED(WMNFTG)) DEALLOCATE(WMNFTG)
          IF(ALLOCATED(UDIRG)) DEALLOCATE(UDIRG)
          IF(ALLOCATED(TAUWG)) DEALLOCATE(TAUWG)
          IF(ALLOCATED(SMEANG)) DEALLOCATE(SMEANG)
          IF(ALLOCATED(U10G)) DEALLOCATE(U10G)
          IF(ALLOCATED(MWP1G)) DEALLOCATE(MWP1G)
          IF(ALLOCATED(MWP2G)) DEALLOCATE(MWP2G)
          IF(ALLOCATED(WSPRDG)) DEALLOCATE(WSPRDG)
          IF(ALLOCATED(WHGTAG)) DEALLOCATE(WHGTAG)
          IF(ALLOCATED(CWHTAG)) DEALLOCATE(CWHTAG)
          IF(ALLOCATED(RANGAG)) DEALLOCATE(RANGAG)
          IF(ALLOCATED(C4G)) DEALLOCATE(C4G)
          IF(ALLOCATED(BFG)) DEALLOCATE(BFG)
          IF(ALLOCATED(QPG)) DEALLOCATE(QPG)
          IF(ALLOCATED(DEPTHG)) DEALLOCATE(DEPTHG)
          IF(ALLOCATED(HMAXG)) DEALLOCATE(HMAXG)
          IF(ALLOCATED(TMAXG)) DEALLOCATE(TMAXG)
          IF(ALLOCATED(USTOKESG)) DEALLOCATE(USTOKESG)
          IF(ALLOCATED(VSTOKESG)) DEALLOCATE(VSTOKESG)
          IF(ALLOCATED(UCURG)) DEALLOCATE(UCURG)
          IF(ALLOCATED(VCURG)) DEALLOCATE(VCURG)

          IF(ALLOCATED(WHGTSG)) DEALLOCATE(WHGTSG)
          IF(ALLOCATED(WDIRSG)) DEALLOCATE(WDIRSG)
          IF(ALLOCATED(WMNFSG)) DEALLOCATE(WMNFSG)
          IF(ALLOCATED(WHGTWG)) DEALLOCATE(WHGTWG)
          IF(ALLOCATED(WDIRWG)) DEALLOCATE(WDIRWG)
          IF(ALLOCATED(WMNFWG)) DEALLOCATE(WMNFWG)

        ENDIF


!         PRINT TIME.

        WRITE(IU06,112) CDTPRO
        IF (ITEST .GE. 1 .OR. .NOT.LWCOU) CALL FLUSH (IU06)
  112   FORMAT(/,3X,'!!!!!!!!!!!!!! ',
     &   'WAVE FIELDS INTEGRATED.  DATE IS: ',A14,
     &   '  !!!!!!!!!!!! ')

!*    1.8 OUTPUT FILES AND RECOVERY FILES ARE DISPOSED WHEN
!         TIME REACHES THE DISPOSE DATE OR WHEN THE MODEL
!         HAS BEEN SIGNALLED TO DO SO.
!         -------------------------------------------------

        IF (ITEST.GE.1) THEN
          WRITE(IU06,*) " WAMODEL: 1.8  CDTRES=", CDTRES, " CDTPRO=",
     &     CDTPRO, " LDWRRE=", LDWRRE, "KADV=", KADV, " NADV=", NADV
          WRITE(IU06,*)
     &     " CDATEF=", CDATEF, "CDATEE=", CDATEE, " CDATER=", CDATER,
     &     " CDATES=", CDATES, "IREST=", IREST
        ENDIF

        IF (CDATEE.EQ.CDTPRO.OR.LDWRRE) THEN

!*    1.8.2 STORE 2-D SPECTRA AT SELECTED GRID POINTS,
!           INTEGRATED PARAMETERS OF ENTIRE GRID.
!           ------------------------------------------

!           SAVE INTEGRATED DATA FILE.

          IF (GFLAG20 .AND. .NOT. LFDB .AND. LOUTINT) THEN
            IF((LMESSPASS.AND.(IRANK.EQ.1))
     &       .OR.(.NOT.LMESSPASS)) THEN
              CALL GSFILE (IU06, IU30, 0, CDTPRO, CDATEF, 'MPP', 'S')
            ENDIF
          ENDIF
          IF (FFLAG20 .AND. LOUTINT) THEN
            IF((LMESSPASS.AND.(IRANK.EQ.NPROC)).OR.(.NOT.LMESSPASS)) 
     &       THEN
              CALL GSFILE (IU06, IU20, 0, CDTPRO, CDATEF, 'MAP', 'S')
            ENDIF
          ENDIF

!           SAVE INTEGRATED SWELL AND SEA DATA FILE.
!           ----------------------------------------

          IF (FFLAG21 .AND. LOUTINT) THEN
            IF((LMESSPASS.AND.(IRANK.EQ.MIN(2,NPROC)))
     &         .OR.(.NOT.LMESSPASS)) THEN
              CALL GSFILE (IU06, IU21, 0, CDTPRO, CDATEF, 'SWE', 'S')
            ENDIF
          ENDIF

!           SAVE SPECTRA FILE.

          IF (FFLAG25 .AND. LOUTSPP) THEN
            IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
              CALL GSFILE (IU06, IU25, 0, CDTPRO, CDATEF, 'OUT', 'S')
            ENDIF
          ENDIF

!           SAVE SWELL SPECTRA FILE.

          IF (FFLAG26 .AND. LOUTSPP) THEN
            IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
              CALL GSFILE (IU06, IU26, 0, CDTPRO, CDATEF, 'SWS', 'S')
            ENDIF
          ENDIF

 
        END IF

!NEST
!       SAVE BOUNDARY VALUE FILE.

        IF(CDTBC.EQ.CDTPRO) THEN
          IF (IBOUNC.EQ.1) THEN
            IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
              DO II=1,GBOUNC
              CALL GSFILE(IU06, IU19(II), 0, CDTBC, CDTBC,
     &          CBCPREF(II), 'S')
              IF (CDTBC.LT.CDATEE)
     &          CALL HEADBC (IPOGBO(II)-IPOGBO(II-1), IDELPRO,
     &                       TH(1), FR(1), IU19(II), IU06) 
              ENDDO
            ENDIF
          ENDIF
        ENDIF

!NEST
 
        LSV=(CDTRES.EQ.CDTPRO.OR.CDATEE.EQ.CDTPRO.OR.CDTPRO.EQ.CDATER)

        IF (LSV.OR.LRST) THEN


!*        1.8.3 SAVE RESTART FIELDS.
!               --------------------
!         THIS WILL HAPPEN WHEN IT IS NOT IN DATA ASSIMILATION MODE AND
!         IT IS EITHER A DETERMINED OUTPUT TIME
!         OR THE INTERUPT SIGNAL HAS BEEN TRIGGERED and it will wait
!         until the end of the advection loop.
!         OTHERWISE THE OUTPUT WILL OCCUR IN WAMASSI.

          LOUT=((IREST.EQ.1)
     &            .AND.
     &          (CDTPRO.EQ.CDATER .OR. CDTPRO.LE.CDATES))
     &            .AND.
     &             LSV


          IF ( LOUT .OR. LRST ) THEN


!*          STRESS FIELD IS DUMPED FOR RESTART.
!           -----------------------------------
!           NOTE :
!           STRESS FIELD IS DUMPED FOR RESTART
!           ONLY IF CURRENT TIME CORRESPONDS TO CDATER
!           EXCEPT IF CDATER SET TO 00000000000000

            IF((CDTPRO .EQ. CDATER ) .OR.
     &         (CDATER .EQ. '00000000000000' ) .OR.
     &         (LDWRRE)) THEN

              IF ( .NOT.LGRIBOUT .OR.LDWRRE ) THEN
                CALL SAVSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,ROAIRO,
     &            ZIDLOLD,CICOVER,NSTART,NEND,CDTPRO,CDATEF)

              ENDIF

            ENDIF

!           SAVE SPECTRUM
!           -------------

!           we have insured that the spectra will be written to FDB
!           even when the restart option is triggered and it is an
!           output step for the spectra.
            IF(LOUT.AND.LGRIBOUT) THEN

!             IF THE OUTPUT TIME IS NOT AN ANLYSIS TIME THEN TYPE FG
!             BECOMES TYPE AN (i.e. speudo analysis)
              MARSTYPEBAK=MARSTYPE
              IF((MARSTYPE.EQ.'fg' .OR. MARSTYPE.EQ.'4v') .AND.
     &            KADV.LT.NADV) THEN
                MARSTYPE='an'
              ENDIF
              IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
!               convert output spectra to absolute frame of reference
                IF (.NOT.LLINTPOL) THEN
                  CALL GSTATS(1236,0)
                  IRA=1
                  NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1
!$OMP          PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
                  DO JKGLO=IJS(IG),IJL(IG),NPROMA
                    KIJS=JKGLO
                    KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
                    CALL INTPOL (FL1, SL, KIJS, KIJL, IG, IRA)
                  ENDDO
!$OMP          END PARALLEL DO
                  CALL GSTATS(1236,1)
                ENDIF
                IF(LGRIBAPI_O) THEN
                  CALL OUTWSPEC(SL)
                ELSE
                  CALL SPEC2FDB (IU06, NPROC, IRANK, LWCOU, TH, FR,
     &                       SL, NINF, NSUP, NANG, NFRE, NGX, NGY,
     &                       NLONRGG, KXLT, IXLG, NIBLO, NBLO,
     &                       AMONOP, XDELLA, NSTART, NEND,
     &                       NKSEK0, NKSEK1, NKSEK2, NKSEK3, NKSEK4,
     &                       NPSEK2, NPSEK3, KSEK0, KSEK1, KSEK2, KSEK3,
     &                       KSEK4, PSEK2, PSEK3, NTENCODE, HOPERS,
     &                       LGRHDIFS, LNEWLVTP, IFS_KSEC1,
     &                       DEG, PPMISS, PPEPS, PPREC,NGRBRESS,PPRESOL,
     &                       CDTPRO, CDATEF, CDATEA,
     &                       CLDOMAIN, MARSTYPE, YCLASS, YEXPVER,
     &                       CFDB2DSP, NENSFNB, NTOTENS, NSYSNB, NMETNB,
     &                       IREFDATE, LFDBIOOUT,ICPLEN, CPATH, ITEST,
     &                       NWFDBREF, LFDBOPEN)
                ENDIF
              ELSE
                IF(LGRIBAPI_O) THEN
                  CALL OUTWSPEC(FL1)
                ELSE
                  CALL SPEC2FDB (IU06, NPROC, IRANK, LWCOU, TH, FR,
     &                       FL1, NINF, NSUP, NANG, NFRE, NGX, NGY,
     &                       NLONRGG, KXLT, IXLG, NIBLO, NBLO,
     &                       AMONOP, XDELLA, NSTART, NEND,
     &                       NKSEK0, NKSEK1, NKSEK2, NKSEK3, NKSEK4,
     &                       NPSEK2, NPSEK3, KSEK0, KSEK1, KSEK2, KSEK3,
     &                       KSEK4, PSEK2, PSEK3, NTENCODE, HOPERS,
     &                       LGRHDIFS, LNEWLVTP, IFS_KSEC1,
     &                       DEG, PPMISS, PPEPS, PPREC,NGRBRESS,PPRESOL,
     &                       CDTPRO, CDATEF, CDATEA,
     &                       CLDOMAIN, MARSTYPE, YCLASS, YEXPVER,
     &                       CFDB2DSP, NENSFNB, NTOTENS, NSYSNB, NMETNB,
     &                       IREFDATE, LFDBIOOUT,ICPLEN, CPATH, ITEST,
     &                       NWFDBREF, LFDBOPEN)
                  ENDIF
              ENDIF
              MARSTYPE=MARSTYPEBAK

              IF (LFDB) THEN
                CALL GSTATS(1976,0)
                ISTAT = IFLUSHFDBSUBS (NWFDBREF)
!J              CALL MPL_BARRIER(CDSTRING='WAMODEL:')
                CALL GSTATS(1976,1)
              ELSE
                CALL GSTATS(753,0)
!J              CALL MPL_BARRIER(CDSTRING='WAMODEL:')
                CALL GSTATS(753,1)
              ENDIF
              LLFLUSH = .FALSE.
            ENDIF  

            IF ( .NOT.LGRIBOUT .OR. LDWRRE ) THEN
              CALL SAVSPEC(FL3,NSTART,NEND,CDTPRO,CDATEF,CDATER)
            ENDIF

            WRITE(IU06,*) ' '
            WRITE(IU06,*) '  WAVE SPECTRA FILES DISPOSED AT........',
     &       ' CDTPRO  = ', CDTPRO
            WRITE(IU06,*) ' '


!*    1.8.4 UPDATE, WRITE AND SAVE WAMINFO FILE.
!           ------------------------------------

            IF (LRST .AND.
     &       (LMESSPASS .AND. IRANK .EQ. 1 .OR. .NOT. LMESSPASS )) THEN
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
              ISTAT = 0
              IF (CDATE.EQ.CDATEE) ISTAT(1) = 1
              IDELWIN = IDELWI

              CEPLTDT = CDATEF

              IU04 =  I_GET_UNIT (IU06,CWI(1:ICPLEN+8) , 'w', 'f', 0)

              CALL WRITSTA (IU04, CDTPRO, CDATEE, IANALPD, IFOREPD,
     &                      IDELWIN, CDATER, CDATES, CBPLTDT, CEPLTDT,
     &                      IASSI, NFCST, ISTAT, CDTCUR)
 
              CLOSE (IU04)
              WRITE(IU06,*) ' WAMINFO FILE WRITTEN FOR RESTART...',
     &         ' CDTPRO  = ', CDTPRO
              WRITE(IU06,*) '                                    ',
     &         ' CDATEF  = ', CDATEF
              WRITE(IU06,*) ' TO ', CWI(1:ICPLEN+8)
              CALL FLUSH(IU06)
            ENDIF

          ENDIF  ! END SAVE RESTART FIELDS

        ENDIF  ! END SAVE

        IF( LFDB .AND. LLFLUSH .AND. NWFDBREF.NE.-5 .AND. 
     &    (IASSI.NE.1 .OR. CDTPRO.GT.CDATEF) .AND.
     &    (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO) ) THEN
          CALL GSTATS(1976,0)
          ISTAT = IFLUSHFDBSUBS (NWFDBREF)
!J        CALL MPL_BARRIER(CDSTRING='WAMODEL:')
          CALL GSTATS(1976,1)
          WRITE(IU06,*) ' ' 
          WRITE(IU06,*) '  DB ', NWFDBREF , ' FLUSHED AT ',
     &    CDTPRO, ' FROM WAMODEL. '
          CALL FLUSH (IU06)
        ELSE
          CALL GSTATS(753,0)
!J        CALL MPL_BARRIER(CDSTRING='WAMODEL:')
          CALL GSTATS(753,1)
        ENDIF

!*    BRANCHING BACK TO 1.0 FOR NEXT PROPAGATION STEP.

      ENDDO ADVECTION

      IF(ALLOCATED(FL4)) DEALLOCATE(FL4)

      IF (LHOOK) CALL DR_HOOK('WAMODEL',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE WAMODEL
