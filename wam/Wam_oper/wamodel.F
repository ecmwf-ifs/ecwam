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
!          WAVE SPECTRUM IN TIME DPENDANT CURRENT AND DEPTH FIELD USING
!          NONCANONICAL SPACIAL (SPHERICAL) AND WAVE NUMBER (FRQUENCY-
!          DIRECTION) COORDINATES. WAM REPORT 1988.

! -------------------------------------------------------------------

      USE YOWCBOU  , ONLY : IBOUNC   ,NBOUNC
      USE YOWCOUP  , ONLY : LWCOU
      USE YOWCOUT  , ONLY : COUTT    ,COUTS    ,FFLAG20  ,FFLAG21  ,
     &            FFLAG25  ,FFLAG26  ,PFLAG20  ,PFLAG21  ,PFLAG25  ,
     &            PFLAG26  ,GFLAG20  ,GFLAG21  ,GFLAG25  ,GFLAG26  ,
     &            JPPFLAG  ,GFLAG    ,LFDB     ,NOUTT    ,NOUTS
      USE YOWFBOU  , ONLY : IBOUNF
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3   ,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK0    ,KSEK1    ,
     &            KSEK2    ,KSEK3    ,KSEK4    ,PSEK2    ,PSEK3    ,
     &            PPMISS   ,PPEPS    ,PPREC    ,NTENCODE ,NGRBRESS ,
     &            HOPERS
      USE YOWGRID  , ONLY : NLONRGG  ,IGL      ,IJS      ,IJL2     ,
     &            IJL      ,IJLT
      USE YOWICE   , ONLY : LICERUN
      USE YOWINTP  , ONLY : WHGTTG   ,WDIRTG   ,WPKFTG   ,WMNFTG   ,
     &            USTARG   ,UDIRG    ,TAUWG    ,CDG      ,SMEANG   ,
     &            U10G     ,MWP1G    ,MWP2G    ,WSPRDG
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
      USE YOWPCONS , ONLY : ZMISS    ,DEG
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTH    ,DEPTHA   ,DEPTHD   ,
     &            INDEP
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDTPRO   ,CDTRES   ,
     &            CDATER   ,CDATES   ,CDTINTT  ,CDTINTS  ,CDTSPT   ,
     &            CDTSPS   ,CFDB2DSP ,MFDBSF   ,IDELPRO  ,IDELT    ,
     &            IDELWI   ,IREST    ,IDELRES  ,IDELINT  ,IDELINS  ,
     &            IDELSPT  ,IDELSPS  ,ISHALLO  ,IREFRA   ,IASSI    ,
     &            NTASKS   ,NSIZE    ,NENSFNB  ,NTOTENS  ,NSYSNB   ,
     &            NMETNB   ,CDATEA   ,MARSTYPE ,YCLASS   ,YEXPVER  ,
     &            LANAONLY ,LALTCOR  ,LFRSTFLD ,NWAM_BLKS
      USE YOWSPEC, ONLY : NSTART   ,NEND     ,
     &            NTOPE  ,NTOPEMAX ,IJTOPE   ,NGBTOPE  ,NTOPELST   ,
     &            NFROMPE,NFROMPEMAX,NIJSTART,NGBFROMPE,NFROMPELST ,
     &            U10NEW   ,U10OLD   ,THWNEW   ,THWOLD   ,USNEW    ,
     &            USOLD    ,Z0NEW    ,Z0OLD    ,TAUW     ,FCONST   ,
     &            ROAIRN   ,ROAIRO   ,ZIDLNEW  ,ZIDLOLD  ,
     &            FL1      ,FL2      ,FL3      ,SL
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH    ,CWI
      USE YOWUBUF  , ONLY : KLAT     ,KLON
      USE YOWUNIT  , ONLY : IU02     ,IU08     ,IU11     ,IU12     ,
     &            IU13     ,IU14     ,IU19     ,IU20     ,IU21     ,
     &            IU25     ,IU26     ,IU04     ,IU30
      USE YOWWAMI  , ONLY : CBPLTDT  ,CEPLTDT  ,IANALPD  ,IFOREPD  ,
     &            IDELWIN  ,NFCST    ,ISTAT
      USE YOWWIND  , ONLY : CDATEWO
      USE MPL_MODULE

! ----------------------------------------------------------------------
      CHARACTER*2 MARSTYPEBAK
      CHARACTER*12 CDATE, CDATEWH, CDTPRA, ZERO
      CHARACTER*12 CDTINTTBAK,CDTINTSBAK,CDTSPTBAK,CDTSPSBAK 
      LOGICAL LDSTOP, LDWRRE, LLFLUSH
      LOGICAL GFLAGBAK(JPPFLAG)
!     DEFINE LOCAL ALLOCATABLE ARRAYS

      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: FL4

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

      KTAG=200

      IF(CDTPRO.EQ.CDATEA .OR. CDTPRO.EQ.CDATEF) THEN

        CDTINTTBAK=CDTINTT
        CDTINTSBAK=CDTINTS
        CDTSPTBAK=CDTSPT
        CDTSPSBAK=CDTSPS
        CDTINTT=CDTPRO
        CDTINTS=CDTPRO
        CDTSPT=CDTPRO
        CDTSPS=CDTPRO

        DO IG=1,IGL

          IGH1=(IG/2)*2

          DO IJ=IJS(IG),IJL(IG)
            U10NEW(IJ) = U10OLD(IJ,IG)
            THWNEW(IJ) = THWOLD(IJ,IG)
            USNEW(IJ) = USOLD(IJ,IG)
            ROAIRN(IJ) = ROAIRO(IJ,IG)
            ZIDLNEW(IJ) = ZIDLOLD(IJ,IG)
          ENDDO

          IF (IGH1.EQ.IG) THEN
!REFRA
            IF (IREFRA.EQ.2) THEN
              CALL INTPOL (FL1, SL, IJS(IG), IJL(IG), IG)
              CALL OUTBS (SL, FL2, IJS(IG), IJL(IG), IG, IGL,
     &                    IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,
     &                    TAUW, NSTART,NEND)
!REFRA
            ELSE
              CALL OUTBS (FL1, FL2, IJS(IG), IJL(IG), IG, IGL,
     &                    IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,
     &                    TAUW, NSTART,NEND)
            ENDIF
          ELSE
!REFRA
            IF (IREFRA.EQ.2) THEN
              CALL INTPOL (FL1, SL, IJS(IG), IJL(IG), IG)
              CALL OUTBS (SL, FL1, IJS(IG), IJL(IG), IG, IGL,
     &                    IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,
     &                    TAUW, NSTART,NEND)
!REFRA
            ELSE
              CALL OUTBS (FL1, FL3, IJS(IG), IJL(IG), IG, IGL,
     &                    IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,
     &                    TAUW, NSTART,NEND)
            ENDIF
          ENDIF

!         SET FL1 ON ICE POINTS TO ZERO

          IF (LICERUN) THEN
            IF (ITEST.GE.1) THEN
              WRITE(IU06,*) '   SUB. WAMODEL: SPECTRUM = 0 AT',
     &         ' ICEPOINTS'
            ENDIF
            CALL SETICE(FL1,IJS(IG),IJL(IG),IG)
          ENDIF


        ENDDO

!       PRINT OUT NORMS

        CALL OUTWNORM

!*    0.1 SAVE INITIAL INTEGRATED FIELDS
!         ------------------------------
        IF(MARSTYPE.EQ.'cf'.OR.MARSTYPE.EQ.'pf'.OR.MARSTYPE.EQ.'fc'.OR.
     &     MARSTYPE.EQ.'sf'.OR.MARSTYPE.EQ.'4v'.OR.
     &     LANAONLY.OR.LFRSTFLD) THEN

          IF(LFRSTFLD) THEN
            GFLAGBAK=GFLAG
            MARSTYPEBAK=MARSTYPE
            IF(MARSTYPE.EQ.'fg'.OR.MARSTYPE.EQ.'an') THEN
              IF(.NOT.LNOCDIN) GFLAG(7)=.FALSE.
              IF(LWAVEWIND) GFLAG(10)=.FALSE.
              GFLAG(22)=.FALSE.
              GFLAG(23)=.FALSE.
              GFLAG(24)=.FALSE.
              MARSTYPE='an'
            ENDIF
          ENDIF

          CALL OUTINT (IU06, IU20, IU21, IU30, ITEST, MFDBSF)

          IF (MFDBSF.NE.-5.AND.LFDB) ISTAT = IFLUSHFDB (MFDBSF)

          IF(LFRSTFLD) THEN
            MARSTYPE=MARSTYPEBAK
            GFLAG=GFLAGBAK
            LFRSTFLD=.FALSE.
          ENDIF

          CDTINTT=CDTINTTBAK
          CDTINTS=CDTINTSBAK

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

        IF(ALLOCATED(WHGTSG)) DEALLOCATE(WHGTSG)
        IF(ALLOCATED(WDIRSG)) DEALLOCATE(WDIRSG)
        IF(ALLOCATED(WMNFSG)) DEALLOCATE(WMNFSG)
        IF(ALLOCATED(WHGTWG)) DEALLOCATE(WHGTWG)
        IF(ALLOCATED(WDIRWG)) DEALLOCATE(WDIRWG)
        IF(ALLOCATED(WMNFWG)) DEALLOCATE(WMNFWG)

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
          IF ((PFLAG20.OR.FFLAG20.OR.GFLAG20).AND.CDTINTT.LT.CDTPRO)
     &     CALL INCDATE (CDTINTT,IDELINT)
          IF ((PFLAG21.OR.FFLAG21.OR.GFLAG21).AND.CDTINTS.LT.CDTPRO)
     &     CALL INCDATE (CDTINTS,IDELINS)
          IF ((PFLAG25.OR.FFLAG25.OR.GFLAG25).AND.CDTSPT.LT.CDTPRO)
     &     CALL INCDATE (CDTSPT,IDELSPT)
          IF ((PFLAG26.OR.FFLAG26.OR.GFLAG26).AND.CDTSPS.LT.CDTPRO)
     &     CALL INCDATE (CDTSPS,IDELSPS)
        ENDIF

!     UPDATE SPECTRA OUTPUT DATE
        IF (NOUTS.GT.0) THEN
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

!*    1.2 RESET FILES.
!         ------------
        IF(.NOT.LMESSPASS) REWIND IU12

        IF (IGL.NE.1) THEN
          REWIND IU08
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
            READ (IU08) KLAT, KLON
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


!*    1.5.5.1 READING OF NEW WINDS IF NEEDED.
!             -------------------------------

          CALL NEWWIND(IJS(IG),IJL(IG),IG,IGL,CDATE,CDATEWH,
     &                 NEWREAD,NEWFILE,U10OLD,THWOLD,U10NEW,THWNEW
     &               , ROAIRO, ROAIRN, ZIDLOLD,ZIDLNEW)
          IF (ITEST.GE.2) THEN
            IF (ITESTB.GE.IG) THEN
              WRITE(IU06,*) '   SUB. WAMODEL: NEWWIND CALLED' 
            ENDIF
          ENDIF

!*    1.5.5.2 COMPUTATION OF PROPAGATION &&& 
!*            INTEGRATION OF SOURCE TERMS OVER SUB TIME STEPS BETWEEN
!*            PROPAGATION TIME STEPS.
!             -------------------------------------------------------

! DETERMINE WHAT FORM OF PARALELLISM IS USED BASED ON THE TYPE
! OF MACHINE ARCHTECTURE

          IF(LMESSPASS) THEN

            IF(IGH1.NE.IG) THEN
! IF CDATE CORRESPONDS TO A PROPAGATION TIME , THEN CALL MPEXCHNG
! TO OBTAIN INFORMATION AT NEIGHBORING GRID POINTS LOCATED ON ADJACENT
! PROCESSES.
              IF(CDATE.EQ.CDTPRA) THEN
                CALL MPEXCHNG(FL1,NSTART,NEND,
     &                        NTOPE,NTOPEMAX,IJTOPE,NGBTOPE,NTOPELST,
     &                        NFROMPE,NFROMPEMAX,NIJSTART,NGBFROMPE,
     &                        NFROMPELST)

                IF (ITEST.GE.2) THEN
                  WRITE(IU06,*) '   SUB. WAMODEL: MPEXCHNG CALLED' 
                  CALL FLUSH (IU06)
                ENDIF
              ENDIF

! Mod for OPENMP
! If NWAM_BLKS=1 then the code will be as before

              NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1

              IF (CDATE.EQ.CDTPRA) THEN
!$OMP           PARALLEL PRIVATE(JKGLO,KIJS,KIJL)
!$OMP           DO
                DO JKGLO=IJS(IG),IJL(IG),NPROMA
                  KIJS=JKGLO
                  KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
                  IF(NIBLO.GT.1) THEN
                    CALL PROPAGS(FL1,FL3,KIJS,KIJL,IG,SL,FCONST)
                    IF (ITEST.GE.2) THEN
                      WRITE(IU06,*) '   SUB. WAMODEL: PROPAGS CALLED'
                      CALL FLUSH (IU06)
                    ENDIF
                  ELSE
                    FL3 = FL1
                  ENDIF
                ENDDO
!$OMP           END DO
!$OMP           END PARALLEL
              ENDIF

!$OMP         PARALLEL PRIVATE(JKGLO,KIJS,KIJL)
!$OMP         DO
              DO JKGLO=IJS(IG),IJL(IG),NPROMA
                KIJS=JKGLO
                KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
                CALL IMPLSCH (FL3, FL1, KIJS, KIJL, IG,
     &                        THWOLD,USOLD,TAUW,Z0OLD,   ROAIRO,ZIDLOLD,
     &                        U10NEW,THWNEW,USNEW,Z0NEW, ROAIRN,ZIDLNEW,
     &                        SL,FCONST)
                IF (ITEST.GE.2) THEN
                  WRITE(IU06,*) '   SUB. WAMODEL: IMPLSCH CALLED'
                  CALL FLUSH (IU06)
                ENDIF

              ENDDO
!$OMP         END DO
!$OMP         END PARALLEL
! End Mod for OPENMP

            ELSE
              WRITE (IU06,*) '**************************************'
              WRITE (IU06,*) '*                                    *'
              WRITE (IU06,*) '*    FATAL ERROR                     *'
              WRITE (IU06,*) '*    ===========                     *'
              WRITE (IU06,*) '* A MULTI BLOCK VERSION WAS REQUESTED*'
              WRITE (IU06,*) '* THIS OPTION IS NOT IMPLEMENTED FOR *'
              WRITE (IU06,*) '* MESSAGE PASSING PROTOCOL           *'
              WRITE (IU06,*) '*                                    *'
              WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.  *'
              WRITE (IU06,*) '*                                    *'
              WRITE (IU06,*) '**************************************'
              CALL ABORT1
            ENDIF

          ELSE
            NCHUNK = 0
            NLENGTH = (IJL(IG)-IJS(IG)+1)
            NWORK   = NLENGTH/NSIZE
            IF (MOD(NLENGTH,NSIZE).NE.0) NWORK = NWORK + 1

            NTSKCAP = 3
            NSTKSZ  = 1
            NSTKINC = 1
            DO JROC=1,NTASKS-1
              ITSKCA(1,JROC) = NTSKCAP
              ITSKCA(3,JROC) = JROC
              ITSKCA(4,JROC) = 1 ! !!!! JUST TO GET PAST COMP 'STACKSZW'
              ITSKCA(5,JROC) = NSTKSZ
              ITSKCA(6,JROC) = 1 ! !!!! JUST TO GET PAST COMP 'STACKINW'
              ITSKCA(7,JROC) = NSTKINC

              IF(IGH1.NE.IG) THEN
                CALL TSKSTART(ITSKCA(1,JROC),PHYS,FL3,FL1,FL4,IJS(IG),
     &                        IJL(IG),IG,NSIZE,NWORK,CDATE,CDTPRA)
              ELSE
                CALL TSKSTART(ITSKCA(1,JROC),PHYS,FL3,FL2,FL4,IJS(IG),
     &                        IJL(IG),IG,NSIZE,NWORK,CDATE,CDTPRA)
              ENDIF
            ENDDO
            IF (ITEST.GE.2) CALL FLUSH(IU06)

            IF(IGH1.NE.IG) THEN
              CALL PHYS (FL3,FL1,FL4,IJS(IG),IJL(IG),IG,NSIZE,NWORK,
     &                   CDATE, CDTPRA,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     &                   ROAIRO, ZIDLOLD,
     &                   U10NEW,THWNEW,USNEW,Z0NEW,
     &                   ROAIRN, ZIDLNEW, SL,FCONST)
            ELSE
              CALL PHYS (FL3,FL2,FL4,IJS(IG),IJL(IG),IG,NSIZE, NWORK,
     &                   CDATE, CDTPRA,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     &                   ROAIRO, ZIDLOLD,
     &                   U10NEW,THWNEW,USNEW,Z0NEW,
     &                   ROAIRN, ZIDLNEW, SL,FCONST)
            ENDIF

            DO JROC=1,NTASKS-1
              CALL TSKWAIT(ITSKCA(1,JROC))
            ENDDO

          ENDIF



!*    1.5.5.3 COPY AND CLOSE FILES IF NEEDED.
!             -------------------------------

          CALL CLOSEND(IJS(IG),IJL(IG),IG,IGL,CDATE,CDATEWH,
     &                 NEWREAD,NEWFILE,U10OLD,THWOLD,ROAIRO,ZIDLOLD,
     &                                 U10NEW,THWNEW,ROAIRN,ZIDLNEW)
          IF (ITEST.GE.2) THEN
            IF (ITESTB.GE.IG) THEN
              WRITE(IU06,*) '   SUB. WAMODEL: CLOSEND CALLED' 
              CALL FLUSH (IU06)
            ENDIF
          ENDIF

!*    1.5.5.4 UPDATE TIME;IF TIME LEFT BRANCH BACK TO 1.5.5 
!             ---------------------------------------------

          CALL INCDATE(CDATE,IDELT)   
          IF (CDATE.LT.CDTPRO) GO TO 1550

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

          IF (LICERUN) THEN
            IF (ITEST.GE.1) THEN
              WRITE(IU06,*) '   SUB. WAMODEL: SPECTRUM = 0 AT',
     &         ' ICEPOINTS'
            ENDIF
            CALL SETICE(FL3,IJS(IG),IJL(IG),IG)
          ENDIF

!ICE
!NEST

!*    1.5.7 INPUT OF BOUNDARY VALUES.
!           -------------------------

          IF (IBOUNF.EQ.1) THEN
            CALL BOUINPT (FL3, IG, IU02)
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

!!          IF(.NOT.LWCOU .AND. .NOT. LDSTOP) THEN
!!!!      the call to CHESIG is a signal handeling facility which is
!!!!      specific to running WAM at ECMWF, it can be ignored when
!!!!      WAM is not run at ECMWF.
!!            CALL CHESIG (IU06, ITEST, IRANK, NPROC, LDSTOP, LDWRRE)
!!          ENDIF

          IF(LDWRRE) THEN
            WRITE(IU06,*) ' '
            WRITE(IU06,*) '  ******************************************'
            WRITE(IU06,*) '  AN INTERRUPT SIGNAL HAS BEEN RECEIVED '
            WRITE(IU06,*) '  THE NECESSARY BINARY RESTART FILES WILL BE'
            WRITE(IU06,*) '  GENERATED BY SETTING LGRIBOUT AND '
            WRITE(IU06,*) '  LFDBIOOUT TO FALSE. '
            WRITE(IU06,*) '  ******************************************'
            WRITE(IU06,*) ' '
            CALL FLUSH (IU06)
            LGRIBOUT = .FALSE.
            LFDBIOOUT = .FALSE. 
          ENDIF

          IF (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO .OR.
     &     CDTSPT.EQ.CDTPRO .OR. CDTSPS.EQ.CDTPRO .OR.
     &     CDATEE.EQ.CDTPRO .OR. (LDWRRE .AND. KADV.EQ.NADV ) .OR.
     &     (IASSI.EQ.1 .AND. KADV.EQ.NADV ))
     &     THEN

            IF (IGH1.EQ.IG) THEN
!REFRA
              IF (IREFRA.EQ.2) THEN
                CALL INTPOL (FL3, SL, IJS(IG), IJL(IG), IG)
                CALL OUTBS (SL, FL2, IJS(IG), IJL(IG), IG, IGL,
     &                      IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,
     &                      TAUW, NSTART,NEND)

!REFRA
              ELSE
                CALL OUTBS (FL3, FL2, IJS(IG), IJL(IG), IG, IGL,
     &                      IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,
     &                      TAUW, NSTART,NEND)
              ENDIF
            ELSE
!REFRA
              IF (IREFRA.EQ.2) THEN
                CALL INTPOL (FL3, SL, IJS(IG), IJL(IG), IG)
                CALL OUTBS (SL, FL1, IJS(IG), IJL(IG), IG, IGL,
     &                      IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,
     &                      TAUW, NSTART,NEND)
!REFRA
              ELSE
                CALL OUTBS (FL3, FL1, IJS(IG), IJL(IG), IG, IGL,
     &                      IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,
     &                      TAUW, NSTART,NEND)
              ENDIF
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

        IF (IGL.EQ.1) THEN
          DO M=1,NFRE
            DO K=1,NANG
              DO IJ=NINF-1,NSUP
                FL1(IJ,K,M) = FL3(IJ,K,M)
!!!   note in case of Alt and SAR assimilation both FL1 and FL3
!!!   are needed since FL1 will be first modified by the altimeter
!!!   data assimilation and FL3 is needed to keep the first guess
!!!   (unaltered by the alt data assimilation) for the calculation
!!!   of the difference between SAR data and model values
!!!   (see rearrangsar)
              ENDDO
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

          CALL OUTINT (IU06, IU20, IU21, IU30, ITEST, MFDBSF)
          IF(LMESSPASS) THEN
            CALL MPL_BARRIER(CDSTRING='WAMODEL:')
          ENDIF

        ENDIF

!       DEALLOCATE ARRAYS WHICH WERE ALLOCATED IN OUTBS
        IF (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO .OR.
     &   CDTSPT.EQ.CDTPRO .OR. CDTSPS.EQ.CDTPRO .OR.
     &   CDATEE.EQ.CDTPRO .OR. (LDWRRE .AND. KADV.EQ.NADV ) .OR.
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

          IF(ALLOCATED(WHGTSG)) DEALLOCATE(WHGTSG)
          IF(ALLOCATED(WDIRSG)) DEALLOCATE(WDIRSG)
          IF(ALLOCATED(WMNFSG)) DEALLOCATE(WMNFSG)
          IF(ALLOCATED(WHGTWG)) DEALLOCATE(WHGTWG)
          IF(ALLOCATED(WDIRWG)) DEALLOCATE(WDIRWG)
          IF(ALLOCATED(WMNFWG)) DEALLOCATE(WMNFWG)

        ENDIF


!         PRINT TIME.

        IF (ITEST .GE. 1) WRITE(IU06,112) CDTPRO
        IF (ITEST .GE. 2) CALL FLUSH (IU06)
  112   FORMAT(/,3X,'!!!!!!!!!!!!!! ',
     &   'WAVE FIELDS INTEGRATED  DATE IS: ',A12,
     &   '  !!!!!!!!!!!!!! ')

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

        IF (CDATEE.EQ.CDTPRO.AND.LDWRRE) THEN

!*    1.8.2 STORE 2-D SPECTRA AT SELECTED GRID POINTS,
!           INTEGRATED PARAMETERS OF ENTIRE GRID.
!           ------------------------------------------

!           SAVE INTEGRATED DATA FILE.

          IF (GFLAG20 .AND. .NOT. LFDB) THEN
            IF((LMESSPASS.AND.(IRANK.EQ.1))
     &       .OR.(.NOT.LMESSPASS)) THEN
              CALL GSFILE (IU06, IU30, 0, CDTPRO, CDATEF, 'MPP', 'S')
            ENDIF
          ENDIF
          IF (FFLAG20) THEN
            IF((LMESSPASS.AND.(IRANK.EQ.NPROC)).OR.(.NOT.LMESSPASS)) 
     &       THEN
              CALL GSFILE (IU06, IU20, 0, CDTPRO, CDATEF, 'MAP', 'S')
            ENDIF
          ENDIF

!           SAVE INTEGRATED SWELL AND SEA DATA FILE.
!           ----------------------------------------

          IF (FFLAG21) THEN
            IF((LMESSPASS.AND.(IRANK.EQ.MIN(2,NPROC)))
     &         .OR.(.NOT.LMESSPASS)) THEN
              CALL GSFILE (IU06, IU21, 0, CDTPRO, CDATEF, 'SWE', 'S')
            ENDIF
          ENDIF

!           SAVE SPECTRA FILE.

          IF (FFLAG25) THEN
            IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
              CALL GSFILE (IU06, IU25, 0, CDTPRO, CDATEF, 'OUT', 'S')
            ENDIF
          ENDIF

!           SAVE SWELL SPECTRA FILE.

          IF (FFLAG26) THEN
            IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
              CALL GSFILE (IU06, IU26, 0, CDTPRO, CDATEF, 'SWS', 'S')
            ENDIF
          ENDIF

C
      END IF
C
      IF (CDTRES.EQ.CDTPRO.OR.CDATEE.EQ.CDTPRO.OR.CDTPRO.EQ.CDATER.OR.
     1    (LDWRRE .AND. KADV.EQ.NADV)) THEN

!NEST

!           SAVE BOUNDARY VALUE FILE.

          IF (IBOUNC.EQ.1) THEN
            CALL GSFILE (IU06, IU19, 0, CDTPRO, CDATEF, 'CBO', 'S')
            IF (CDTPRO.LT.CDATEE)
     &         CALL HEADBC (NBOUNC, IDELPRO, TH(1), FR(1), IU19, IU06)
          ENDIF
!NEST

!*    1.8.3 SAVE RESTART FIELDS.
!           --------------------
!     THIS WILL HAPPEN WHEN IT IS NOT IN DATA ASSIMILATION MODE AND
!     IT IS EITHER A DETERMINED OUTPUT TIME
!     OR THE INTERUPT SIGNAL HAS BEEN TRIGGERED and it will wait until
!     the end of the advection loop.
!     OTHERWISE THE OUTPUT WILL OCCUR IN WAMASSI.


          IF ( (IASSI.NE.1 .OR. CDTPRO.GT.CDATEF)
     &         .AND.
     &         ( ( (IREST.EQ.1)
     &             .AND.
     &             (CDTPRO.EQ.CDATER .OR. CDTPRO.LE.CDATES)
     &           )
     &           .OR.
     &           (LDWRRE .AND. KADV.EQ.NADV)  
     &         )
     &       ) THEN


!*       STRESS FIELD IS DUMPED FOR RESTART.
!        -----------------------------------
!        NOTE :
!        STRESS FIELD IS DUMPED FOR RESTART
!        ONLY IF CURRENT TIME CORRESPONDS TO CDATER
!        EXCEPT IF CDATER SET TO 000000000000

            IF((CDTPRO .EQ. CDATER ) .OR.
     &       (CDATER .EQ. '000000000000' ) .OR.
     &       (LDWRRE)) THEN

              IF ( .NOT.LGRIBOUT ) THEN
                CALL SAVSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,ROAIRO,
     &                         ZIDLOLD,NSTART,NEND,CDTPRO,CDATEF)
              ENDIF

            ENDIF

!         SAVE SPECTRUM
!         -------------

            IF ( LGRIBOUT ) THEN
              CALL SPEC2FDB (IU06, NPROC, IRANK, LWCOU, TH, FR,
     &                       FL1, NINF, NSUP, NANG, NFRE, NGX, NGY,
     &                       NLONRGG, KXLT, IXLG, NIBLO, NBLO,
     &                       AMONOP, XDELLA, NSTART, NEND,
     &                       NKSEK0, NKSEK1, NKSEK2, NKSEK3, NKSEK4,
     &                       NPSEK2, NPSEK3, KSEK0, KSEK1, KSEK2, KSEK3,
     &                       KSEK4, PSEK2, PSEK3, NTENCODE, HOPERS,
     &                       DEG, PPMISS, PPEPS, PPREC, NGRBRESS,
     &                       CDTPRO, CDATEF, CDATEA,
     &                       CLDOMAIN, MARSTYPE, YCLASS, YEXPVER,
     &                       CFDB2DSP, NENSFNB, NTOTENS, NSYSNB, NMETNB,
     &                       LFDBIOOUT,ICPLEN, CPATH, ITEST, MFDBSF)
              IF (LFDB) ISTAT = IFLUSHFDB (MFDBSF)
              LLFLUSH = .FALSE.
            ELSE
              CALL SAVSPEC(FL3,NSTART,NEND,CDTPRO,CDATEF,CDATER,
     &                     MARSTYPE,YCLASS, YEXPVER, CFDB2DSP)
            ENDIF
            WRITE(IU06,*) ' '
            WRITE(IU06,*) '  WAVE SPECTRA FILES DISPOSED AT........',
     &       ' CDTPRO  = ', CDTPRO
            WRITE(IU06,*) ' '


!*    1.8.4 UPDATE, WRITE AND SAVE WAMINFO FILE.
!           ------------------------------------

            IF (LDWRRE .AND. KADV.EQ.NADV .AND.
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

              OPEN (IU04,FILE=CWI(1:ICPLEN+8),
     &         STATUS='UNKNOWN',FORM='FORMATTED')
              CALL WRITSTA (IU04, CDTPRO, CDATEE, IANALPD, IFOREPD,
     &                      IDELWIN, CDATER, CDATES, CBPLTDT, CEPLTDT,
     &                      IASSI, NFCST, ISTAT)
 
              CLOSE (IU04)
              WRITE(IU06,*) ' WAMINFO FILE WRITTEN FOR RESTART...',
     &         ' CDTPRO  = ', CDTPRO
              WRITE(IU06,*) '                                    ',
     &         ' CDATEF  = ', CDATEF
              WRITE(IU06,*) ' TO ', CWI(1:ICPLEN+8)
            ENDIF
          ENDIF
        ENDIF
        IF ( LFDB .AND. LLFLUSH .AND. MFDBSF.NE.-5 .AND. 
     &   (IASSI.NE.1 .OR. CDTPRO.GT.CDATEF) .AND.
     &   (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO) ) THEN
          ISTAT = IFLUSHFDB (MFDBSF)
          WRITE(IU06,*) ' ' 
          WRITE(IU06,*) '  DB ', MFDBSF , ' FLUSHED AT ',
     &     CDTPRO, ' FROM WAMODEL. '
          CALL FLUSH (IU06)
        ENDIF

!*    BRANCHING BACK TO 1.0 FOR NEXT PROPAGATION STEP.

      ENDDO ADVECTION

      IF(ALLOCATED(FL4)) DEALLOCATE(FL4)

      RETURN
      END SUBROUTINE WAMODEL
