      SUBROUTINE USERIN (IU05, IBOUNC, IBOUNF, IFORCA)

! ----------------------------------------------------------------------

!**** *USERIN* - ROUTINE TO READ AND WRITE NAMELIST INPUT.

!     H. GUNTHER     GKSS/ECMWF   NOVEMBER 1989
!     J. BIDLOT      ECMWF        JUNE  1996
!     J. BIDLOT      ECMWF        JUNE  1997: coupled/uncoupled
!                                             restarted/normal
!     B. HANSEN      ECMWF        JULY  1997  COUPLED AND ANALYSIS
!                                             (IE. PREPAN)
!     B. HANSEN      ECMWF        JANUARY 1998 NAMELIST INPUT.
!     B. HANSEN      ECMWF        NOVEMBER 1998 WPOSNAM USED TO POSITION
!                                                 NAMLIST FOR READING.
!     S. ABDALLA     ECMWF        OCTOBER 2001  LGUST & LADEN ADDED TO
!                                               CONTROL GUSTINESS & AIR DENSITY 

!*    PURPOSE.
!     --------

!       READ USER INPUT CONCERNING PERIOD OF INTEREST,TIMESTEPS AND
!       MODEL OPTIONS TO INITIALIZE COMMON USERD. A CONSISTENCY CHECK
!       IS DONE TOO.

!**   INTERFACE.
!     ----------

!       *CALL* *USERIN (IU05, IBOUNC, IBOUNF, IFORCA)*
!          *IU05*   -  INPUT UNIT USER INPUT FILE.
!          *IBOUNC* -  COARSE GRID BOUNDARY FLAG.
!          *IBOUNF* -  FINE GRID BOUNDARY FLAG.
!          *IFORCA* -  FORCAST START OPTION.

!     METHOD.
!     -------

!        USER INFORMATION IS BEING READ WITH THE PRESUMPTIONS THAT:
!         1. EVERY LINE STARTING WITH 'C' IS A COMMENT LINE
!         2. VALUES ARE PUT IN BELOW POSITIONS INDICATED WITH '-'
!            (RIGHT-JUSTIFIED)

!     EXTERNALS.
!     ----------

!       *ABORT1*    - TERMINATES PROCESSING.
!       *DIFDATE*   - COMPUTES A TIME DIFFERENCE.
!       *IINITFDB*  - INITIALISE DATA BASE.
!       *WPOSNAM*   - POSITION NAMELIST FOR READING AND CONTROLLED
!                     TRAPPING OF EOF.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWALTAS , ONLY : NUMALT   ,IBUFRSAT  ,ALTSDTHRSH,ALTBGTHRSH,
     &            HSALTCUT, LALTGRDOUT,
     &            XKAPPA2  ,HSCOEFCOR,HSCONSCOR ,LALTCOR   ,LALTLRGR
      USE YOWCOUP  , ONLY : LWCOU    ,LWCOUNORMS
      USE YOWCOUT  , ONLY : COUTT    ,COUTS    ,CASS     ,FFLAG    ,
     &            FFLAG20  ,FFLAG21  ,FFLAG25  ,FFLAG26  ,PFLAG    ,
     &            PFLAG20  ,PFLAG21  ,PFLAG25  ,PFLAG26  ,GFLAG    ,
     &            GFLAG20  ,GFLAG21  ,GFLAG25  ,GFLAG26  ,NFLAG    ,
     &            NFLAGALL ,UFLAG    ,LFDB     ,NOUTT    ,NOUTS    ,
     &            NASS
      USE YOWGRIBHD, ONLY : IFS_KSEC1, LGRHDIFS ,LNEWLVTP
      USE YOWMESPAS, ONLY : LMESSPASS,LPBIOIN  ,LPBIOOUT ,LFDBIOIN ,
     &            LFDBIOOUT,LGRIBIN  ,LGRIBOUT ,LNOCDIN
      USE YOWPARAM , ONLY : SWAMPWIND,SWAMPWIND2,DTNEWWIND,LTURN90 ,
     &            LL1D
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDATER   ,CDATES   ,
     &            CFDBSF   ,CFDB2DSP ,IDELPRO  ,IDELT    ,IDELWI   ,
     &            IDELWO   ,IDELALT  ,IREST    ,IDELRES  ,IDELINT  ,
     &            IDELINS  ,IDELSPT  ,IDELSPS  ,ICASE    ,ISHALLO  ,
     &            IREFRA   ,IASSI    ,NTASKS   ,NSIZE    ,NENSFNB  ,
     &            NTOTENS  ,NSYSNB   ,NMETNB   ,CDATEA   ,YCLASS   ,
     &            YEXPVER  ,L4VTYPE  ,LFRSTFLD ,LALTAS   ,LSARAS   ,
     &            LSARINV  ,ISTREAM  ,NLOCGRB  ,NCONSENSUS,NDWD    ,
     &            NMFR     ,NNCEP    ,NUKM     ,IREFDATE ,
     &            LGUST    ,LADEN    ,NWAM_BLKS           
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWTEXT  , ONLY : LRESTARTED,ICPLEN   ,USERID   ,RUNID    ,
     &            PATH     ,CPATH    ,CWI
      USE YOWUNIT  , ONLY : IU20     ,IU21     ,IU23     ,IU25     ,
     &            IU26     ,IU27     ,IU28     ,IU04     ,IU30     ,
     &            IU31     ,IU32     ,IU33     ,IU35     ,IU36     ,
     &            IU37     ,IU38
      USE YOWWAMI  , ONLY : CBEGDT   ,CENDDT   ,CBPLTDT  ,CEPLTDT  ,
     &            CLSPDT   ,CRSTDT   ,IANALPD  ,IFOREPD  ,IDELWIN  ,
     &            IASSIM   ,NFCST    ,ISTAT
      USE FDBSUBS

! ----------------------------------------------------------------------

      LOGICAL LERROR, LMESSPIFS, LFDBOPIFS, LLEOF, LASTREAM
      LOGICAL LLNALTGO
      CHARACTER*2 :: MARSFCTYPE
      CHARACTER*4 :: CSTREAM
      CHARACTER (LEN=  12), PARAMETER :: ZERO = ' '
      CHARACTER (LEN=  70)            :: CLHEADER
      CHARACTER (LEN= 256)            :: CLFORM

! ----------------------------------------------------------------------

      CHARACTER*1 CLMTSU(4), CLOTSU(7)

!     caution : if you change NALINE, also change it in mpuserin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      NAMELIST /NALINE/ CLHEADER,
     &   CBPLTDT, CEPLTDT, CDATEF,
     &   IDELPRO, IDELT, IDELWO, IDELWI, CLMTSU, IDELALT,
     &   IDELINT, IDELINS, IDELAGI, IDELSPT, IDELSPS, IDELAGS, IDELRES,
     &   CLOTSU, CDATER, CDATES,
     &   PFLAG,  FFLAG,  GFLAG, NFLAG,
     &   LFDB, LPBIOIN, LGRIBIN, LPBIOOUT, LGRIBOUT, LFDBIOOUT,
     &   ICASE, ISHALLO, ITEST, ITESTB, IREST, IASSI, IFORCA,
     &   IBOUNC, IBOUNF,
     &   USERID, RUNID,  PATH, YCLASS, YEXPVER, CPATH,
     &   NENSFNB, NTOTENS, NSYSNB, NMETNB,
     &   LMESSPASS, LWCOU, NTASKS, NSIZE, CFDBSF, CFDB2DSP, LNOCDIN,
     &   LALTCOR, L4VTYPE, LFRSTFLD, LALTAS, LSARAS, LSARINV, XKAPPA2,
     &   IBUFRSAT, SWAMPWIND, SWAMPWIND2, DTNEWWIND, LTURN90,
     &   LALTLRGR, HSCOEFCOR, HSCONSCOR,ALTSDTHRSH,ALTBGTHRSH,HSALTCUT,
     &   ISTREAM, NLOCGRB, IREFDATE,
     &   NCONSENSUS, NDWD, NMFR, NNCEP, NUKM,
     &   LGUST, LADEN, LALTGRDOUT,
     &   NWAM_BLKS, LL1D, LWCOUNORMS, LGRHDIFS ,LNEWLVTP

      CHARACTER* 12 CLOUT
      NAMELIST /NAOT/ CLOUT

      CHARACTER* 12 CLSOUT
      NAMELIST /NAOS/ CLSOUT

      CHARACTER* 12 CLAOUT
      NAMELIST /NAAT/ CLAOUT

!     NAMELIST NALINE : 
!     ===============
!     CBPLTDT: USER INPUT START DATE OF RUN.
!     CEPLTDT: USER INPUT END DATE OF RUN.
!     CDATEF: BEGIN DATE OF FORECAST.
!     IDELPRO: PROPAGATION TIME STEP.
!     IDELT: SOURCE TERM INTEGRATION TIME STEP.
!     IDELWO: WIND OUTPUT TIME STEP.
!     IDELWI: WIND INPUT TIME STEP.
!     CLMTSU: STEP UNIT (S : seconds or H : hours).
!     IDELALT: ALTIMETER DATA TIME WINDOW (in seconds).
!     IDELINT: OUTPUT TIME STEP FOR INTEGRATED PARAMETER OF TOTAL SEA
!              (if output times are not specified by namelist NAOT).
!     IDELINS: OUTPUT TIME STEP INTEGRATED PARAMETER OF SWELL/WINDSEA
!              (if output times are not specified by namelist NAOT).
!     IDELAGI: TIME STEP NOT IN USE.
!     IDELSPT: OUTPUT TIME STEP for SPECTRUM OF TOTAL SEA AT SELECTED 
!              GRID POINTS (see preproc).
!     IDELSPS: OUTPUT TIME STEP for SPECTRUM OF SWELL AT SELECTED
!              GRID POINTS (see preproc).
!     IDELAGS: TIME STEP NOT IN USE.
!     IDELRES: OUTPUT TIME STEP for RESTART SPECTRA (if output times
!              are not specified by namelist NAOS).
!     CLOTSU: STEP UNIT (S : seconds or H : hours) FOR IDELINT, IDELINS,
!             IDELAGI, IDELSPT, IDELSPS, IDELAGS, IDELRES.
!     CDATER: ONE SPECIFIC OUTPUT TIME FOR RESTART FILES if output times
!             are not specified by namelist NAOS).
!     CDATES: LAST OUTPUT TIME ALLOWED BESIDE CDATER FOR RESTART SPECTRA
!             (if output times are not specified by namelist NAOS).
!     PFLAG: OUTPUT FLAG FOR OUTPUT TO FILE OF EACH OUTPUT TYPE.
!     FFLAG: OUTPUT FLAG FOR OUTPUT TO USER OUTPUT OF EACH OUTPUT TYPE.
!     GFLAG: OUTPUT FLAG FOR OUTPUT TO GRIB OF EACH OUTPUT TYPE.
!     NFLAG: OUTPUT FLAG FOR USER OUTPUT DISPLAY OF FIELD NORM FOR
!            EACH OUTPUT TYPE.
!     TYPE OF INTEGRATED PARAMETERS IN PFLAG FFLAG GFLAG (see OUTINT) :
!     1  : WAVE HEIGHT (M)
!     2  : MEAN WAVE DIRECTION (DEG.)
!     3  : WAVE MEAN PERIOD BASED ON THE 1/f INTEGRATION OF F (s.)
!     4  : FRICTION VELOCITY (M/S) (not in grib) 
!     5  : WIND DIRECTION (DEG.)
!     6  : PEAK WAVE PRERIOD (s)
!     7  : DRAG COEFFICIENT (-)
!     8  : NORMALISED WAVE STRESS (-) (not not grib)
!     9  : MEAN SQUARE SLOPE (-)
!     10 : WIND SPEED (M/S)
!     11 : WIND SEA WAVE HEIGHT (M)
!     12 : SWELL WAVE HEIGHT (M)
!     13 : MEAN WIND SEA DIRECTION (DEG.)
!     14 : MEAN SWELL DIRECTION (DEG.)
!     15 : MEAN WIND SEA PERIOD (s)
!     16 : MEAN SWELL PERIOD (s) 
!     17 :
!     18 :
!     19 :
!     20 :
!     21 :
!     22 : ALTIMETER WAVE HEIGHT (M)
!     23 : CORRECTED ALTIMETER WAVE HEIGHT (M)
!     24 : ALTIMETER RANGE RELATIVE CORRECTION (-)
!     25 : MEAN PERIOD BASED ON f * F INTEGRATION OF F (s.) 
!     26 : MEAN PERIOD BASED ON f**2 * F INTEGRATION OF F (s.) 
!     27 : MEAN DIRECTIONAL SPREAD (-) 
!     28 : WIND SEA MEAN PERIOD BASED ON f * F INTEGRATION OF F (s.) 
!     29 : SWELL MEAN PERIOD BASED ON f * F INTEGRATION OF F (s.) 
!     30 : WINDSEA MEAN PERIOD BASED ON f**2 * F INTEGRATION OF F (s.) 
!     31 : SWELL MEAN PERIOD BASED ON f**2 * F INTEGRATION OF F (s.) 
!     32 : WINDSEA MEAN DIRECTIONAL SPREAD (-) 
!     33 : SWELL MEAN DIRECTIONAL SPREAD (-) 
!     34 : KURTOSIS DERIVED FROM WAVE SPECTRA (-)
!     35 : BENJAMIN-FEIR INDEX (-) 
!     36 : PEAKEDNESS FACTOR OR GODA QUALITY FACTOR (-) 
!     37 : MODEL BATHYMETRY (m)

!     LFDB: TRUE IF INTEGRATED PARAMETERS ARE WRITTEN OUT TO FDB.
!     LPBIOIN: IF TRUE THE PBIO SOFTWARE IS USED FOR BINARY INPUT
!              (see LGRIBIN).
!     LGRIBIN : IF TRUE THE WAVE SPECTRA IS INPUT IN GRIB FORMAT, ELSE
!               THE BINARY RESTART FILES ARE USED.
!     LPBIOOUT : IF TRUE THE PBIO SOFTWARE WILL BE USED FOR BINARY
!                OUTPUT (see LGRIBOUT).
!     LGRIBOUT : IF TRUE THE WAVE SPECTRA IS OUTPUT IN GRIB FORMAT, ELSE
!                THE BINARY RESTART FILES ARE PRODUCED.
!     LFDBIOOUT : IF TRUE THE GRIB SPECTRA OUTPUT IS SENT TO THE FDB.
!     ICASE: 1 FOR SPHERICAL COORDINATES ELSE CARTESIAN COORDINATES.
!     ISHALLO: 1 FOR DEEP WATER MODE ELSE SHALLOW WATER MODEL.
!     ITEST: TEST OUTPUT LEVEL.
!     ITESTB: MAX BLOCK NUMBER FOR OUTPUT IN BLOCK LOOPS.
!     IREST: 1 FOR THE PRODUCTION OF RESTART FILE(S).
!     IASSI: 1 ASSIMILATION IS DONE IF ANALYSIS RUN.
!     IFORCA: 1 FOR START RUN WITH ANALYSIS FIELD (OBSOLETE). 
!     IBOUNC: 1 FOR RUN WHICH INCLUDES BOUNDARY POINTS (COARSE GRID).
!     IBOUNF: 1 FOR RUN WHICH INCLUDES BOUNDARY POINTS (FINE GRID).
!     USERID, RUNID, PATH : OUT OF DATE, SHOULD ONLY BE USED WHENEVER
!                           USER WANTS TO WRITE OUTPUT DIRECTLY TO ECFS
!                           WHICH IS HIGHLY NOT RECOMMENDATED ON THE VPP
!     YCLASS: DATA CLASS USED TO CODE DATA IN GRIB.
!     YEXPVER: EXPERIMENT VERSION USED TO CODE DATA IN GRIB.
!     CPATH: PATH FOR OUTPUT TO DISK (USED WHENEVER OUTPUT IS NOT
!            TO THE FIELD DATA BASE (SEE LFDBIOOUT)).
!     NTOTENS: TOTAL ENSEMBLE FORECAST MEMBERS (DEFAULT=0).
!     NENSFNB: ENSEMBLE FORECAST NUMBER (DEFAULT=0).
!     NSYSNB : SYSTEM NUMBER TO BE USED FOR GRIBBING OF SEASONAL DATA.
!              or MONTHLY FORECAST RUNS
!     NMETNB : METHOD NUMBER TO BE USED FOR GRIBBING OF SEASONAL DATA.
!              or MONTHLY FORECAST RUNS
!     LMESSPASS: TRUE FOR MESSAGE PASSING ARCHITECHTURE.
!     LWCOU: FALSE FOR UNCOUPLED RUN.
!     NTASKS: NUMBER OF PROCESSORS FOR PARALELL RUN ON CRAY (OBSOLETE).
!     NSIZE: WORK LOAD PER PROCESSOR FOR PARALELL RUN ON CRAY (OBSOLETE)
!     CFDBSF: FDB ROOT DIRECTORY FOR SCALAR GRIB FIELDS.
!     CFDB2DSP : FDB ROOT DIRECTORY FOR 2D SPECTRA. 
!     LNOCDIN: IF TRUE THEN GRIB INPUT OF A DRAG COEFFICIENT FIELD IS
!              NOT REQUIRED.  
!     LALTCOR: IF TRUE THEN THE ALTIMETER DATA WILL BE CORRECTED
!              SEE GRFIELD. 
!     L4VTYPE: IF TRUE THEN MARS TYPE 4V IS USED INSTEAD OF FG.
!     LFRSTFLD: IF TRUE THEN INITIAL INTEGRATED PARAMETER FIELDS WILL
!               OUTPUT PROVIDED THEY WERE NOT INPUT.
!     XKAPPA2: KAPPA2 PARAMETER USED IN THE FORMULA FOR THE
!              DETERNMINATION OF ALTIMETER WAVE HEIGHT.
!     LALTAS: CONTROLS WHETHER ALTIMETER DATA ARE ASSIMILATED
!     LSARAS: CONTROLS WHETHER SAR DATA ASSIMILATED
!     LSARINV: CONTROLS WHETHER SAR INVERSION IS DONE
!     SWAMPWIND : CONSTANT WIND SPEED USED BY SWAMP CASE.
!     SWAMPWIND2 : SECOND WIND SPEED VALUE TO BE USED (if different than 0)
!                  IN SWAMP CASE. IT WILL BE APPLIED DTNEWWIND HOURS
!                  AFTER STATING TIME.
!     DTNEWWIND : TIME IN HOURS AFTER WHICH SWAMPWIND2 WILL BE APPLIED
!     LTURN90 : IF TRUE THE NEW WIND IN SWAMP CASE WILL TURN BY 90 DEGREES.
!     LALTLRGR: IF TRUE THEN THE ALTIMETER DATA WILL BE CORRECTED
!               PRIOR TO THEIR ASSIMILATION USING A LINEAR REGRESSION AS
!               PROVIDED BY HSCOEFCOR AND HSCONSCOR (see below)
!     HSCOEFCOR: COEFFICIENT OF THE CORRECTIVE LINEAR REGRESSION FOR
!                ALTIMETER WAVE HEIGHTS.
!     HSCONSCOR: CONSTANT OF THE CORRECTIVE LINEAR REGRESSION FOR
!                ALTIMETER WAVE HEIGHTS.
!     ALTSDTHRSH:THRESHOLD FOR SUSPICIOUS DATA (SEE GRFIELD).
!     ALTBGTHRSH:THRESHOLD FOR BACKGROUND CHECK (SEE GRFIELD).
!     HSALTCUT: USER INPUT OF THE MINIMUM WAVE HEIGHT ALLOWED IN ALTAS
!              (SEE GRFIELD). 
!     ISTREAM: STREAM NUMBER USED WHEN GRIBBING THE DATA 
!                     IT MUST BE SPECIFIED EXCEPT IF LGRHDIFS IS TRUE ! 
!     LALTGRDOUT: FLAG USED TO CONTROL THE OUTPUT OF GRIDDED ALTIMETER
!                 PRODUCTS OF SPECIFIC INSTRUMENT (TRUE TO OUTPUT) 
!                (SEE GRFIELD). 
!     LGUST:   FLAG USED TO ACTIVATE COMPUTATIONS RELATED TO GUSTINESS 
!     LADEN:   FLAG USED TO ACTIVATE COMPUTATIONS RELATED TO AIR DENSITY
!     NLOCGRB: LOCAL GRIB TABLE NUMBER
!     IREFDATE: REFERENCE DATE FOR MONTHLY FORECAST HINDCAST.
!     NCONSENSUS: ORIGINE OF THE INITIAL CONDITONS FOR
!                 MULTI-ANALYSIS ENSEMBLE RUNS
!     NDWD: IF 1 DWD ANALYSES ARE USED IN THE MULTI-ANALYSIS ENSEMBLE
!           RUNS
!     NMFR: IF 1 METEO FRANCE ANALYSES ARE USED IN THE MULTI-ANALYSIS
!           ENSEMBLE RUNS
!     NNCEP: IF 1 NCEP ANALYSES ARE USED IN THE MULTI-ANALYSIS ENSEMBLE
!           RUNS
!     NUKM: IF 1 MET OFFICE ANALYSES ARE USED IN THE MULTI-ANALYSIS
!           ENSEMBLE RUNS
!     NWAM_BLKS: IF OPENMP IS USED THEN NWAM_BLKS WILL BE THE NUMBER
!                USED TO SPLIT THE GRID POINTS INTO NWAM_BLKS CHUNKS
!                THAT ARE DISTRIBUTED TO THE DIFFERENT THREADS. 
!                OTHERWISE IT IS STILL USED TO SPLIT THE INNER LOOPS
!                INTO CHUNKS (SEE WAMODEL).
!     LL1D: IF TRUE THEN THE DOMAIN IS ONLY DIVIDED INTO LATITUDINAL BANDS
!     LWCOUNORMS : FLAG THAT CONTROLS WHETHER OR NOT GLOBAL NORMS
!                  ARE PRODUCED FOR THE DIFFERENT FIELDS EXCHANGED 
!                  BETWEEN WAM AND THE COUPLED ATMOSPHERIC MODEL.
!                  IT IS FALSE BY DEFAULT.
!     LGRHDIFS : FLAGS CONTROLLING WHETHER OR NOT GRIB HEADER INFORMATION
!                IS COPIED FROM THE ATMOSPHERIC MODEL (ONLY USEFULL IF
!                COUPLED TO THE IFS).
!     LNEWLVTP : FLAGS CONTROLLING WHETHER OR NOT THE NEW GRIB LEVTYPE
!                DEFINITIONS ARE USED.
!     NAMELIST NAOT : 
!     ===============
!     CLOUT : LIST OF OUTPUT TIME FOR INTEGRATED PARAMETERS.
      
!     NAMELISS NAOT : 
!     ===============
!     CLSOUT : LIST OF OUTPUT TIME FOR SPECTRA. 

!     NAMELISS NAAT : 
!     ===============
!     CLAOUT : LIST OF TIME FOR ALTIMETER DATA ASSIMILATION. 

! ----------------------------------------------------------------------

!*    0. SET DEFAULT VALUES FOR THE NAMELIST ELEMENTS.
!        ---------------------------------------------

      CLMTSU    = 'S'
      CLOTSU    = 'H'
      CLHEADER  = ZERO
      CBPLTDT   =  ZERO
      CEPLTDT   =  ZERO
      CDATEF    =  ZERO
      IDELPRO   =  0
      IDELT     =  0
      IDELWO    =  0
      IDELWI    =  0
      IDELALT   =  21600 
      IDELINT   =  0
      IDELINS   =  0
      IDELAGI   =  0
      IDELSPT   =  0
      IDELSPS   =  0
      IDELAGS   =  0
      IDELRES   =  0
      CDATER    = ZERO
      CDATES    = ZERO
      PFLAG     = .FALSE. 
      FFLAG     = .FALSE. 
!     note: if IREST>0 then GFLAG(7) and GFLAG(10) will always be reset
!     true if grib restart files are requested (LGRIBOUT=.true.).
!     see below.
      GFLAG     = .FALSE. 
      NFLAG     = .FALSE. 
      NFLAG(1)  = .TRUE. 
      NFLAG(6)  = .TRUE. 
      NFLAG(7)  = .TRUE. 
      NFLAG(10) = .TRUE. 
!  SAVE THE VALUE OF LFDB IN LFDBOPIFS AS LFDB SHOULD HAVE BEEN
!  DETERMINED BY THE IFS WHEN COUPLED TO IT, SEE WVALLOC.
!  AFTER READING THE NAMELIST LFDB WILL BE RESET TO LFDBOPIFS IF
!  IT IS A COUPLED RUN.
      LFDBOPIFS = LFDB
      LFDB      = .TRUE. 
      LPBIOIN   = .TRUE. 
      LGRIBIN   = .TRUE. 
!     LFDBIOIN : obsolete, was used to read directly from fdb
!                it should not be used
      LFDBIOIN  = .FALSE. 
      LPBIOOUT  = .TRUE. 
      LGRIBOUT  = .TRUE. 
      LFDBIOOUT = .TRUE. 
      ICASE     = 1 
      ISHALLO   = 0 
      ITEST     = 0 
      ITESTB    = 0 
      IREST     = 0 
      IASSI     = 0 
      IFORCA    = 1 
      IBOUNC    = 0 
      IBOUNF    = 0 
      USERID    = ZERO
      RUNID     = ZERO
      PATH      = ZERO
      CPATH     = ZERO
      YCLASS    = ZERO
      YEXPVER   = ZERO
      NENSFNB   = 0
      NTOTENS   = 0
      NSYSNB    = -1
      NMETNB    = -1
!  SAVE THE VALUE OF LMESSPASS IN LMESSPIFS AS LMESSPASS SHOULD HAVE
!  BEEN DETERMINED BY THE IFS WHEN COUPLED TO IT, SEE WVALLOC.
!  AFTER READING THE NAMELIST LMESSPASS WILL BE RESET TO LMESSPIFS IF
!  IT IS A COUPLED RUN.
      LMESSPIFS = LMESSPASS
      LMESSPASS = .TRUE.
      NTASKS    = 0
      NSIZE     = 0

      NOUTT     = 0

      CFDBSF    = ""
      CFDB2DSP  = ""

      LNOCDIN   = .FALSE. 

      DO ISAT=1,NUMALT
        IBUFRSAT(ISAT) = 0
        LALTCOR(ISAT)   = .FALSE. 
        XKAPPA2(ISAT)   = 0.0700
        LALTLRGR(ISAT)  = .FALSE.
        HSCOEFCOR(ISAT) = 1.
        HSCONSCOR(ISAT) = 0.

        ALTBGTHRSH(ISAT) = 1.
!       if no value is provided in the namelist ALTSDTHRSH will
!       be set in grfield.
        ALTSDTHRSH(ISAT) = -1.

!       HSALTCUT is used in combination with the error estimate of
!       the altimeter data to determine the minimum Hs allowed for
!       altimeter data.
        HSALTCUT(ISAT) = 999999.
        LALTGRDOUT(ISAT)  = .FALSE.
      ENDDO

      L4VTYPE   = .FALSE. 

      LFRSTFLD  = .FALSE.

      LALTAS    = .TRUE. 
      LSARAS    = .FALSE. 
      LSARINV   = .FALSE. 


      SWAMPWIND = 18.45
      SWAMPWIND2 = 0.
      DTNEWWIND = 0.
      LTURN90 = .FALSE.

      ISTREAM   = 0

      NLOCGRB   = 1

      IREFDATE  = 0

      NCONSENSUS= 0
      NDWD      = 0
      NMFR      = 0
      NNCEP     = 0
      NUKM      = 0

      NOUTS     = 0

      NASS      = 0

      LGUST = .FALSE.
      LADEN = .FALSE.

      NWAM_BLKS = 1

      LL1D = .TRUE.

      LWCOU = .FALSE.

      LWCOUNORMS = .FALSE.

      LGRHDIFS = .FALSE.

      LNEWLVTP = .FALSE.
! ----------------------------------------------------------------------

!*    1. READ NAMELIST NALINE.
!        ---------------------

      CALL WPOSNAM (IU05, 'NALINE', LLEOF)
      IF (.NOT. LLEOF) THEN
        READ (IU05, NALINE)
      ELSE
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'+ SUBROUTINE USERIN :                      +'
        WRITE(IU06,*)'+ READ NAMELIST FAILED                     +' 
        WRITE(IU06,*)'+ PROGRAM WILL ABORT                       +'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        CALL ABORT1
      ENDIF
      IF(.NOT.LWCOU) LGRHDIFS = .FALSE.  ! by definition
      IF(LGRHDIFS) THEN
!       GET ISTREAM THAT CORRESPONDS TO IFS_STREAM
        IFS_STREAM=IFS_KSEC1(40)
        IF(.NOT.LNEWLVTP) THEN
          CALL WSTREAM_STRG(IFS_STREAM, CSTREAM, NENSFNB, NTOTENS,
     &                      MARSFCTYPE, ISTREAM, LASTREAM)
          IF(CSTREAM.EQ.'****') THEN
            WRITE(IU06,*) '*****************************************'
            WRITE(IU06,*) ''
            WRITE(IU06,*) ' ERROR IN USERIN !!!!'
            WRITE(IU06,*) ' IFS STREAM UNKNOWN '
            WRITE(IU06,*) ' IFS STREAM = ', IFS_STREAM
            WRITE(IU06,*) ' BUT NOT DEFINED IN WSTREAM_STRG !!!!'
            WRITE(IU06,*) ''
            WRITE(IU06,*) '*****************************************'
            CALL ABORT1
          ENDIF
        ELSE
          ISTREAM=IFS_STREAM
        ENDIF
      ELSEIF(ISTREAM.LE.0) THEN
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'+ SUBROUTINE USERIN :                      +'
        WRITE(IU06,*)'+ READ NAMELIST FAILED                     +' 
        WRITE(IU06,*)'+ ISTREAM MUST BE SPECIFIED > 0 !!!!       +' 
        WRITE(IU06,*)'+ PROGRAM WILL ABORT                       +'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        CALL ABORT1
      ENDIF

      IF(ISTREAM.LE.0) THEN
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'+ SUBROUTINE USERIN :                      +'
        WRITE(IU06,*)'+ READ NAMELIST FAILED                     +' 
        WRITE(IU06,*)'+ ISTREAM MUST BE SPECIFIED > 0 !!!!       +' 
        WRITE(IU06,*)'+ PROGRAM WILL ABORT                       +'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        CALL ABORT1
      ENDIF

      IF(LWCOU) LMESSPASS = LMESSPIFS
      IF(LWCOU) LFDB = LFDBOPIFS

      IF(.NOT.LWCOU) THEN
        CDATEA = CBPLTDT
        CDATEE = CEPLTDT
      ELSE
        CBPLTDT = CDATEA
        CEPLTDT = CDATEF
      ENDIF
      CALL u2l1cr( YCLASS )

!           **** MODEL TIME STEPS ****
      IF (CLMTSU(1) .EQ. 'H') IDELPRO = IDELPRO*3600
      IF (CLMTSU(2) .EQ. 'H') IDELT   = IDELT*3600
      IF (CLMTSU(3) .EQ. 'H') IDELWO  = IDELWO*3600
      IF (CLMTSU(4) .EQ. 'H') IDELWI  = IDELWI*3600
!           **** OUTPUT TIME IN FIXED INTERVALS ****
      IF (CLOTSU(1) .EQ. 'H') IDELINT = IDELINT*3600
      IF (CLOTSU(2) .EQ. 'H') IDELINS = IDELINS*3600
      IF (CLOTSU(3) .EQ. 'H') IDELAGI = IDELAGI*3600
      IF (CLOTSU(4) .EQ. 'H') IDELSPT = IDELSPT*3600
      IF (CLOTSU(5) .EQ. 'H') IDELSPS = IDELSPS*3600
      IF (CLOTSU(6) .EQ. 'H') IDELAGS = IDELAGS*3600
      IF (CLOTSU(7) .EQ. 'H') IDELRES = IDELRES*3600
!           **** CHECK LENGTH OF YEXPVER AND PUT IT RIGHT JUSTIFIED ****
      LYEP=LEN_TRIM(YEXPVER)
      YEXPVER(5-LYEP:4)=YEXPVER(1:LYEP)
      YEXPVER(1:4-LYEP)='0000'
!           **** OUTPUT TIME AT SPECIFIED TIMES ****
      REWIND(IU05)
      SPT0: DO
        CALL WPOSNAM (IU05, 'NAOT', LLEOF)
        IF (LLEOF) EXIT SPT0
        READ (IU05, NAOT, END=1900)
        NOUTT = NOUTT+1
      ENDDO SPT0
1900  CONTINUE
      REWIND(IU05)
      ALLOCATE(COUTT(NOUTT))
      SPT: DO IC=1,NOUTT
        CALL WPOSNAM (IU05, 'NAOT', LLEOF)
        IF (LLEOF) EXIT SPT
        READ (IU05, NAOT, END=1910)
        COUTT(IC) = CLOUT
      ENDDO SPT
1910  CONTINUE
!           **** SPECTRA OUTPUT TIME AT SPECIFIED TIMES ****
      REWIND(IU05)
      SPS0: DO
        CALL WPOSNAM (IU05, 'NAOS', LLEOF)
        IF (LLEOF) EXIT SPS0
        READ (IU05, NAOS, END=1901)
        NOUTS = NOUTS+1
      ENDDO SPS0
1901  CONTINUE
      REWIND(IU05)
      ALLOCATE(COUTS(NOUTS))
      SPS: DO IC=1,NOUTS
        CALL WPOSNAM (IU05, 'NAOS', LLEOF)
        IF (LLEOF) EXIT SPS
        READ (IU05, NAOS, END=1911)
        COUTS(IC) = CLSOUT
      ENDDO SPS
1911  CONTINUE
!           **** ASSIMILATION AT SPECIFIED TIMES ****
      REWIND(IU05)
      SPA0: DO
        CALL WPOSNAM (IU05, 'NAAT', LLEOF)
        IF (LLEOF) EXIT SPA0
        READ (IU05, NAAT, END=1902)
        NASS = NASS+1
      ENDDO SPA0
1902  CONTINUE
      REWIND(IU05)
      ALLOCATE(CASS(NASS))
      SPA: DO IC=1,NASS
        CALL WPOSNAM (IU05, 'NAAT', LLEOF)
        IF (LLEOF) EXIT SPA
        READ (IU05, NAAT, END=1912)
        CASS(IC) = CLAOUT
      ENDDO SPA
1912  CONTINUE

      IF(LWCOU) THEN
!     IF COUPLED MODEL
!     SET TEMPORARILY IDELWO and IDELWI TO IDELPRO
!     to AVOID ALL THE CHECKS. SEE IN INITMDL WHERE THEY ARE RESET
        IDELWO=IDELPRO
        IDELWI=IDELPRO
      ENDIF

!   1 FORMAT(A,A3)
!   2 FORMAT(A,A3,A,I10,A)

!* CHECK FLAG FOR GRIBING AS SOME OPTION ARE NOT IMPLEMENTED
!* IN OUTINT

      IF (GFLAG(4)) THEN
        GFLAG(4) = .FALSE.
        WRITE(IU06,*) ' ******************* NOTE ********************'
        WRITE(IU06,*) ' NO GRIB OUTPUT POSSIBLE FOR FRICTION VELOCITY'
        WRITE(IU06,*) ' NOT IMPLEMENTED IN OUTINT ' 
        WRITE(IU06,*) ' '
      ENDIF
      IF (GFLAG(8)) THEN
        GFLAG(8) = .FALSE.
        WRITE(IU06,*) ' *********** NOTE *************** ' 
        WRITE(IU06,*) ' NO GRIB OUTPUT POSSIBLE FOR' 
        WRITE(IU06,*) ' NORMALISED WAVE STRESS'
        WRITE(IU06,*) ' NOT IMPLEMENTED IN OUTINT ' 
        WRITE(IU06,*) ' '
      ENDIF
      IF (IREST.GT.0.AND.LGRIBOUT.AND..NOT.GFLAG(7)) THEN
        GFLAG(7) = .TRUE.
        WRITE(IU06,*) ' ******** NOTE *** NOTE **** NOTE *************' 
        WRITE(IU06,*) ' YOU HAVE REQUESTED RESTART FILE IN GRIB FORMAT' 
        WRITE(IU06,*) ' BUT FAILED TO ASK FOR OUTPUT OF THE DRAG COEF.'
        WRITE(IU06,*) ' THE MODEL WILL RESET THE OUTPUT OPTION TO GET ' 
        WRITE(IU06,*) ' CD OUT AT TIME PRESCRIBED BY THE OUTPUT OF ALL'
        WRITE(IU06,*) ' OTHER INTEGRATED PARAMETERS.'
        WRITE(IU06,*) ' ******** NOTE *** NOTE **** NOTE *************' 
        WRITE(IU06,*) ' '
      ENDIF
      IF (IREST.GT.0.AND.LGRIBOUT.AND..NOT.GFLAG(10)) THEN
        GFLAG(10) = .TRUE.
        WRITE(IU06,*) ' ******** NOTE *** NOTE **** NOTE *************' 
        WRITE(IU06,*) ' YOU HAVE REQUESTED RESTART FILE IN GRIB FORMAT' 
        WRITE(IU06,*) ' BUT FAILED TO ASK FOR OUTPUT OF U10wave'
        WRITE(IU06,*) ' THE MODEL WILL RESET THE OUTPUT OPTION TO GET ' 
        WRITE(IU06,*) ' u10 OUT AT TIME PRESCRIBED BY THE OUTPUT OF ALL'
        WRITE(IU06,*) ' OTHER INTEGRATED PARAMETERS.'
        WRITE(IU06,*) ' ******** NOTE *** NOTE **** NOTE *************' 
        WRITE(IU06,*) ' '
      ENDIF
      IF (GFLAG(22).AND.(IASSI.NE.1.OR..NOT.LALTAS)) THEN
        GFLAG(22) = .FALSE.
      ENDIF
      IF (GFLAG(23).AND.(IASSI.NE.1.OR..NOT.LALTAS)) THEN
        GFLAG(23) = .FALSE.
      ENDIF
      IF (GFLAG(24).AND.(IASSI.NE.1.OR..NOT.LALTAS)) THEN
        GFLAG(24) = .FALSE.
      ENDIF
      FFLAG20 = FFLAG( 1).OR.FFLAG( 2).OR.FFLAG( 3).OR.FFLAG( 4).OR.
     &          FFLAG( 5).OR.FFLAG( 6).OR.FFLAG( 7).OR.FFLAG( 8).OR.
     &          FFLAG( 9).OR.FFLAG(10).OR.FFLAG(22).OR.FFLAG(23).OR.
     &          FFLAG(24).OR.FFLAG(25).OR.FFLAG(26).OR.FFLAG(27).OR.
     &          FFLAG(34).OR.FFLAG(35).OR.FFLAG(36).OR.FFLAG(37)
      PFLAG20 = PFLAG( 1).OR.PFLAG( 2).OR.PFLAG( 3).OR.PFLAG( 4).OR.
     &          PFLAG( 5).OR.PFLAG( 6).OR.PFLAG( 7).OR.PFLAG( 8).OR.
     &          PFLAG( 9).OR.PFLAG(10).OR.PFLAG(22).OR.PFLAG(23).OR.
     &          PFLAG(24).OR.PFLAG(25).OR.PFLAG(26).OR.PFLAG(27).OR.
     &          PFLAG(34).OR.PFLAG(35).OR.PFLAG(36).OR.PFLAG(37)
      GFLAG20 = GFLAG( 1).OR.GFLAG( 2).OR.GFLAG( 3).OR.GFLAG( 4).OR.
     &          GFLAG( 5).OR.GFLAG( 6).OR.GFLAG( 7).OR.GFLAG( 8).OR.
     &          GFLAG( 9).OR.GFLAG(10).OR.GFLAG(22).OR.GFLAG(23).OR.
     &          GFLAG(24).OR.GFLAG(25).OR.GFLAG(26).OR.GFLAG(27).OR.
     &          GFLAG(34).OR.GFLAG(35).OR.GFLAG(36).OR.GFLAG(37)
      FFLAG21 = FFLAG(11).OR.FFLAG(12).OR.FFLAG(13).OR.FFLAG(14).OR.
     &          FFLAG(15).OR.FFLAG(16).OR.FFLAG(28).OR.FFLAG(29).OR.
     &          FFLAG(30).OR.FFLAG(31).OR.FFLAG(32).OR.FFLAG(33) 
      PFLAG21 = PFLAG(11).OR.PFLAG(12).OR.PFLAG(13).OR.PFLAG(14).OR.
     &          PFLAG(15).OR.PFLAG(16).OR.PFLAG(28).OR.PFLAG(29).OR.
     &          PFLAG(30).OR.PFLAG(31).OR.PFLAG(32).OR.PFLAG(33) 
      GFLAG21 = GFLAG(11).OR.GFLAG(12).OR.GFLAG(13).OR.GFLAG(14).OR.
     &          GFLAG(15).OR.GFLAG(16).OR.GFLAG(28).OR.GFLAG(29).OR.
     &          GFLAG(30).OR.GFLAG(31).OR.GFLAG(32).OR.GFLAG(33)
      FFLAG25 = FFLAG(17)
      PFLAG25 = PFLAG(17)
      GFLAG25 = GFLAG(17)
      FFLAG26 = FFLAG(18)
      PFLAG26 = PFLAG(18)
      GFLAG26 = GFLAG(18)
      NFLAGALL= NFLAG( 1).OR.NFLAG( 2).OR.NFLAG( 3).OR.NFLAG( 4).OR.
     &          NFLAG( 5).OR.NFLAG( 6).OR.NFLAG( 7).OR.NFLAG( 8).OR.
     &          NFLAG( 9).OR.NFLAG(10).OR.NFLAG(11).OR.NFLAG(12).OR.
     &          NFLAG(13).OR.NFLAG(14).OR.NFLAG(15).OR.NFLAG(16).OR.
     &          NFLAG(22).OR.NFLAG(23).OR.NFLAG(24).OR.NFLAG(25).OR.
     &          NFLAG(26).OR.NFLAG(27).OR.NFLAG(28).OR.NFLAG(29).OR.
     &          NFLAG(30).OR.NFLAG(31).OR.NFLAG(32).OR.NFLAG(33).OR.
     &          NFLAG(34).OR.NFLAG(35).OR.NFLAG(36).OR.NFLAG(37)

!     DEFINE GRIDDED OUTPUT ARRAYS WHICH ARE USED IN THE ASSIMILATION
      UFLAG=.FALSE.
      IF(IASSI.EQ.1) THEN
        UFLAG(1)=.TRUE.
      ENDIF

      ICPLEN=LEN_TRIM(CPATH)
      IF(ICPLEN.GT.0.AND.CPATH(ICPLEN:ICPLEN).EQ.'/') THEN
        CPATH=CPATH(1:ICPLEN-1)
        ICPLEN=ICPLEN-1
      ENDIF
      CWI=CPATH(1:ICPLEN)//'/waminfo'

      IF(LFDBIOOUT.AND..NOT.LGRIBOUT) THEN
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'+ LFDBIOOUT = TRUE AND LGRIBOUT = FALSE    +'
        WRITE(IU06,*)'+ IS AN OBSOLETE OPTION:                   +'
        WRITE(IU06,*)'+ LFDBIOOUT IS RESET TO FALSE              +'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
      ENDIF

!*    1.1  READ THE WAMINFO FILE AND OVERWRITE INPUT.
!          ------------------------------------------

      INQUIRE(FILE=CWI,EXIST=LRESTARTED)
      IF (LRESTARTED) THEN
        OPEN  (IU04,FILE=CWI,STATUS='OLD', FORM='FORMATTED')
        CALL READSTA(IU04, CBEGDT, CENDDT, IANALPD, IFOREPD, IDELWIN,
     &               CRSTDT, CLSPDT, CBPLTDT, CEPLTDT, IASSIM, NFCST,
     &               ISTAT)
        CLOSE (IU04,STATUS='KEEP')
        WRITE(IU06,*)'+++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*)'+                                           +'
        WRITE(IU06,*)'+    WARNING MODEL RESTART FROM RC FILE     +'
        WRITE(IU06,*)'+    ==================================     +'
        WRITE(IU06,*)'+                                           +'
        WRITE(IU06,*)'+ RC: ', CWI(1:ICPLEN+8)
        WRITE(IU06,*)'+                                           +'
        WRITE(IU06,*)'+ START DATE OF RUN         CDATEA  = ', CDATEA
        WRITE(IU06,*)'+ START DATE OF RESTART     CBEGDT  = ', CBEGDT
        WRITE(IU06,*)'+ ANALYSIS PERIOD (SECONDS) IANALPD = ', IANALPD
        WRITE(IU06,*)'+ FORECAST PERIOD (SECONDS) IFOREPD = ', IFOREPD
        WRITE(IU06,*)'+ END DATE FROM WAMINFO     CENDDT  = ', CENDDT
        WRITE(IU06,*)'+ END DATE FROM INPUT       CDATEE  = ', CDATEE
        WRITE(IU06,*)'+                                           +'
        WRITE(IU06,*)'+ MODEL CONTINUES USING RESTART DATE.       +'
        WRITE(IU06,*)'+                                           +'
        WRITE(IU06,*)'+++++++++++++++++++++++++++++++++++++++++++++'

        CDATEF = CEPLTDT
        IF (.NOT. LWCOU) THEN
          CDATEA = CBEGDT

          CDATEE = CDATEF
!         IF ( IANALPD.LT.0 ) IFOREPD=IFOREPD+IANALPD*(-1) 
          CALL INCDATE (CDATEE,IFOREPD)
        ENDIF

        CDATER = CRSTDT
        CDATES = CLSPDT
        IDELWI = IDELWIN
        IASSI = IASSIM
        IFORCA = NFCST

        IF (.NOT. LWCOU) THEN
          IF (CDATEE.NE.CENDDT) THEN
            WRITE(IU06,*)'+++++++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+    WARNING ERROR IN SUB. USERIN           +'
            WRITE(IU06,*)'+    ============================           +'
            WRITE(IU06,*)'+ END DATE OF RUN AND RUN LENGTH DO NOT     +'
            WRITE(IU06,*)'+ MATCH IN THE WAM INFO FILE.               +'
            WRITE(IU06,*)'+ START DATE OF RUN         CBEGDT =', CBEGDT
            WRITE(IU06,*)'+ ANALYSIS PERIOD (SECONDS) IANALPD=', IANALPD
            WRITE(IU06,*)'+ FORECAST PERIOD (SECONDS) IFOREPD=', IFOREPD
            WRITE(IU06,*)'+ END DATE FROM WAMINFO     CENDDT =', CENDDT
            WRITE(IU06,*)'+ END DATE COMPUTED         CDATEE =', CDATEE
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+ MODEL CONTINUES USING COMPUTED END DATE.  +'
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+++++++++++++++++++++++++++++++++++++++++++++'
            CENDDT = CDATEE
          ENDIF
        ELSE
          IF ( CDATEA .NE. CBEGDT ) THEN
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+++++++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+   F A T A L   E R R O R  IN SUB. USERIN   +'
            WRITE(IU06,*)'+   =====================================   +'
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+   COUPLED RESTART FAILED.                 +'
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+ START DATE IN WAMINFO DOES NOT CORRESPOND +'
            WRITE(IU06,*)'+ WITH CURRENT DATE OF THE ATMOSPHERE       +'
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+   ABORT SERVICE ROUTINE CALLED BY USERIN  +'
            WRITE(IU06,*)'+   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  +'
            CALL ABORT1
          ENDIF
        ENDIF

        IF(LFDBIOIN) THEN
          WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
          WRITE(IU06,*)'+                                          +'
          WRITE(IU06,*)'+ LFDBIOIN = TRUE IS AN OBSOLETE OPTION    +'
          WRITE(IU06,*)'+ PROGRAM WILL ABORT                       +'
          WRITE(IU06,*)'+                                          +'
          WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
          CALL ABORT1
        ENDIF

        IF(LFDBIOOUT.AND..NOT.LGRIBOUT) THEN
          WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
          WRITE(IU06,*)'+                                          +'
          WRITE(IU06,*)'+ LFDBIOOUT = TRUE AND LGRIBOUT = FALSE    +'
          WRITE(IU06,*)'+ IS AN OBSOLETE OPTION                    +'
          WRITE(IU06,*)'+ PROGRAM WILL ABORT                       +'
          WRITE(IU06,*)'+                                          +'
          WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
          CALL ABORT1
        ENDIF

      ENDIF

      IF(CBPLTDT.GT.CDATEA.AND..NOT.LRESTARTED) THEN
        WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*) '+                                           +'
        WRITE(IU06,*) '+    WARNING IN SUB. USERIN                 +'
        WRITE(IU06,*) '+    ======================                 +'
        WRITE(IU06,*) '+    CBPLTDT GREATER THAN CDATEA            +'
        WRITE(IU06,*) '+    CBPLTDT = ',CBPLTDT
        WRITE(IU06,*) '+    CDATEA  = ',CDATEA
        WRITE(IU06,*) '+    CBPLTDT WAS RESET TO CDATEA            +'
        WRITE(IU06,*) '+                                           +'
        WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++'
        CBPLTDT=CDATEA
      ENDIF

      IF(CDATER .NE. '000000000000') THEN
        IF(CDATER .LT. CDATEA) CDATER = CDATEE 
        IF(CDATER .GT. CDATEE) CDATER = CDATEE 
      ENDIF
      IF(CDATES .LT. CDATEA) CDATES = CDATEE 
      IF(CDATES .GT. CDATEE) CDATES = CDATEE 

!     1.3 CHECK IF IDELALT WAS SET ELSE SET IT TO IDELWO
!         ----------------------------------------------

      IF (IDELALT.EQ.0.AND.IASSI.EQ.1) THEN
        IDELALT=IDELWO
        WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*) '+                                             +'
        WRITE(IU06,*) '+         WARNING IN SUB. USERIN              +'
        WRITE(IU06,*) '+         ======================              +'
        WRITE(IU06,*) '+THE ALTIMETER DATA ASSIMILATION WINDOW LENGTH+'
        WRITE(IU06,*) '+WAS NOT PROVIDED IN WAMINPUT, IT WAS SET TO  +'
        WRITE(IU06,*) '+IDELWO = ', IDELWO
        WRITE(IU06,*) '+                                             +'
        WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
      ENDIF

!*    1.4  PRINT INITIAL CONDITIONS.
!          -------------------------
      WRITE(IU06,*) '  '
      IV = INDEX(CLHEADER,'4.0')
      IF (IV .GT. 30 ) THEN
        CLHEADER(IV:IV+2) = '4.8'
      ENDIF
      WRITE(IU06,*) CLHEADER
      WRITE(IU06,*) ' '
      WRITE(IU06,*) ' COUPLING WITH ATMOS. MODEL (LWCOU) : ',LWCOU
      IF(.NOT.LMESSPASS) THEN
        WRITE(IU06,*) ' USING ',NTASKS,' PROCESSORS'
        WRITE(IU06,*) ' WITH A WORKLOAD OF ',NSIZE,' PER PROCESSOR'
      ENDIF
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' STARTING DATE ........... : ',CDATEA
      WRITE(IU06,*) ' END DATE ................ : ',CDATEE
      IF(CDATEF.LT.CDATEE) THEN
        WRITE(IU06,*) ' FORECAST STARTING DATE    : ',CDATEF
      ENDIF
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' MODEL TIME STEPS:'
      WRITE(IU06,*) ' SOURCE TERM INTEGRATION TIME STEP : ',
     & IDELT,' SECS'
      WRITE(IU06,*) ' PROPAGATION TIME STEP ............: ',
     & IDELPRO,' SECS'
      IF(.NOT.LWCOU) THEN
        WRITE(IU06,*) ' MODEL WIND INPUT TIME STEP .......: ',
     &   IDELWI,' SECS'
        WRITE(IU06,*) ' MODEL WIND OUTPUT TIME STEP.......: ',
     &   IDELWO,' SECS'
      ENDIF
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' MODEL OPTIONS:'
      WRITE(IU06,*) '  '
      IF (ICASE.EQ.1) THEN
        WRITE(IU06,*) ' PROPAGATION GRID SPHERICAL LAT/LON COORDINATES'
      ELSE
        WRITE(IU06,*) ' PROPAGATION GRID CARTESIAN COORDINATES'
      ENDIF
      IF(LL1D) THEN
        WRITE(IU06,*) ' 1D DECOMPOSITION OF THE DOMAIN ' 
      ELSE
        WRITE(IU06,*) ' 2D DECOMPOSITION OF THE DOMAIN ' 
      ENDIF
      IF (ISHALLO.EQ.1) THEN
        WRITE(IU06,*) ' THIS IS A DEEP WATER RUN '
      ELSE
        WRITE(IU06,*) ' THIS IS A SHALLOW WATER RUN '
      ENDIF
      IF (IREFRA.EQ.0) THEN
        WRITE(IU06,*) ' MODEL RUNS WITHOUT REFRACTION'
      ELSEIF (IREFRA.EQ.1) THEN
        WRITE(IU06,*) ' MODEL RUNS WITH DEPTH REFRACTION'
      ELSEIF (IREFRA.EQ.2) THEN
        WRITE(IU06,*) ' MODEL RUNS WITH DEPTH AND CURRENT REFRACTION'
      ENDIF
      IF (IDELWO.GE.IDELWI) THEN
        WRITE(IU06,*) ' WIND FIELDS ARE NOT INTERPOLATED IN TIME'
      ELSE
        WRITE(IU06,*) ' WIND FIELDS ARE INTERPOLATED IN TIME'
      ENDIF
      IF (LGUST) THEN
        WRITE(IU06,*) ' GUSTINESS EFFECT IS INCLUDED.'
      ENDIF
      IF (LADEN) THEN
        WRITE(IU06,*) ' VARIABLE AIR DENSITY EFFECT IS INCLUDED.'
      ENDIF
      IF(.NOT.LWCOU) THEN
        IF (IBOUNC.EQ.1) THEN
          WRITE(IU06,*) ' MODEL RUNS WITH BOUNDARY POINTS (COARSE GRID)'
        ELSE
          WRITE(IU06,*) ' MODEL RUNS WITHOUT BOUNDARY POINTS',
     &     ' (COARSE GRID)'
        ENDIF
        IF (IBOUNF.EQ.1) THEN
          WRITE(IU06,*) ' MODEL RUNS WITH BOUNDARY POINTS (FINE GRID)'
        ELSE
          WRITE(IU06,*) ' MODEL RUNS WITHOUT BOUNDARY POINTS',
     &     ' (FINE GRID)'
        ENDIF
        WRITE(IU06,*) ' '
        IF (IFORCA.EQ.1) THEN
          WRITE(IU06,*) ' MODEL STARTS FROM ANALYSIS FIELDS'
        ELSE
          WRITE(IU06,*) ' MODEL STARTS FROM FORECAST FIELDS'
        ENDIF
      ENDIF
      IF (IASSI.EQ.1) THEN
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' WAVE DATA ASSIMILATION IS CARRIED OUT'
        IF(NASS.GT.0) THEN
          WRITE(IU06,*) ' AT DATE(S) '
          WRITE(IU06,'(2X,A12)') (CASS(I),I=1,NASS)
        ELSE
          IF(LWCOU) THEN
            WRITE(IU06,*) ' AT DATE ', CDATEF
          ELSE
            WRITE(IU06,*) ' UNTIL THE END OF THE ANALYSIS PERIOD'
            WRITE(IU06,*) ' AT DATE ', CDATEF
          ENDIF
          CALL FLUSH(IU06)
        ENDIF
        IF(LALTAS) THEN
          WRITE(IU06,*) ' WITH ALTIMETER DATA IN TIME WINDOW(S) OF ' 
          WRITE(IU06,*) ' IDELALT = ', IDELALT,' SECONDS' 
          WRITE(IU06,*) ' CENTERED AROUND THE ASSIMILATION TIME(S) '
          LLNALTGO = .TRUE.
          DO ISAT=1,NUMALT
             WRITE(IU06,*) '  '
             WRITE(IU06,*) ' THE ALTIMETER DATA FROM SATELLITE ',
     &                       IBUFRSAT(ISAT) 
            IF(LALTLRGR(ISAT)) THEN
              WRITE(IU06,*) ' THE DATA WILL BE CORRECTED ' 
              WRITE(IU06,*) ' ACCORDING TO FOLLOWING LINEAR REGRESSION'
              WRITE(IU06,*) ' Hsnew= ',HSCOEFCOR(ISAT),' Hs + ',
     &                                 HSCONSCOR(ISAT) 
            ENDIF
            IF(LALTCOR(ISAT)) THEN
              WRITE(IU06,*) '  '
              WRITE(IU06,*) ' THE DATA WILL BE CORRECTED ' 
              WRITE(IU06,*) ' ACCORDING TO THE MODEL SEA STATE.' 
            ENDIF
            WRITE(IU06,*) ' THE THRESHOLD FOR BACKGROUND CHECK IS ',
     &                      ALTBGTHRSH(ISAT)
            IF(HSALTCUT(ISAT).LT.999999.) THEN
              WRITE(IU06,*) ' THE INPUT MINIMUM WAVE HEIGHT IS ',
     &                        HSALTCUT(ISAT)
            ELSE
              WRITE(IU06,*) ' THE MINIMUM WAVE HEIGHT WILL BE',
     &                      ' THE OBSERVATION ERROR.'
            ENDIF
            IF(LALTGRDOUT(ISAT)) THEN
              WRITE(IU06,*) ' GRIDDED ALTIMETER FIELDS WILL BE' 
              WRITE(IU06,*) ' PRODUCED FOR THIS ALTIMETER.' 
              LLNALTGO = .FALSE.
            ENDIF
            WRITE(IU06,*) '  '
          ENDDO
          IF (LLNALTGO) THEN
            WRITE(IU06,*) '  '
            WRITE(IU06,*) ' WARNING   WARNING   WARNING   WARNING'
            WRITE(IU06,*) ' GRIDDED ALTIMETER FIELDS WILL NOT BE' 
            WRITE(IU06,*) ' PRODUCED FOR ANY INSTRUMENT.'
            WRITE(IU06,*) '  '
          ENDIF
        ENDIF
        IF(LSARAS) THEN
          WRITE(IU06,*) ' '
          WRITE(IU06,*) ' WITH SAR DATA IN TIME WINDOW(S) OF ' 
          WRITE(IU06,*) ' IDELALT = ', IDELALT,' SECONDS' 
          WRITE(IU06,*) ' CENTERED AROUND THE ASSIMILATION TIME(S) '
        ENDIF
      ELSE
        WRITE(IU06,*) ' DATA ASSIMILATION IS NOT CARRIED OUT'
      ENDIF
      IF(LSARINV) THEN
        WRITE(IU06,*) '  '
        WRITE(IU06,*) ' SAR INVERSION IS CARRIED OUT'
        WRITE(IU06,*) '  '
      ENDIF
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' MODEL OUTPUT SELECTION:'
      WRITE(IU06,*) '  '
      IF (NOUTT.GT.0) THEN
        WRITE(IU06,*) ' NUMBER OF OUTPUT TIMES IS NOUTT = ', NOUTT
        WRITE(IU06,*) ' OUTPUT WILL BE PROCESSED AT:'
        WRITE(IU06,'(6(2X,A12))') (COUTT(I),I=1,NOUTT)
        WRITE(IU06,*) '  '
      ENDIF
      IF(FFLAG20) CALL OUTXT(IU06,IU20,IU30,IU20,IU32,IDELINT,IASSI,
     & NOUTT,'IU20','IU30','IU20','IU32')

      WRITE(IU06,*) '    F = FALSE   T = TRUE '
      WRITE(IU06,*) '                           PRINTER   FILE    GRIB'
      WRITE(IU06,*) '1  SIGNIFICANT WAVE HEIGHT ...', PFLAG( 1),
     & '......', FFLAG( 1),'.....',GFLAG( 1)
      WRITE(IU06,*) '2  MEAN WAVE DIRECTION .......', PFLAG( 2),
     & '......', FFLAG( 2),'.....',GFLAG( 2)
      WRITE(IU06,*) '3  WAVE MEAN PERIOD (-1) .....', PFLAG( 3),
     & '......', FFLAG( 3),'.....',GFLAG( 3)
      WRITE(IU06,*) '4  FRICTION VELOCITY .........', PFLAG( 4),
     & '......', FFLAG( 4),'.....',GFLAG( 4)
      WRITE(IU06,*) '5  WAVE MODEL WIND DIRECTION..', PFLAG( 5),
     & '......', FFLAG( 5),'.....',GFLAG( 5)
      WRITE(IU06,*) '6  WAVE PEAK FREQUENCY .......', PFLAG( 6),
     & '......', FFLAG( 6),'.....',GFLAG( 6)
      WRITE(IU06,*) '7  DRAG COEFFICIENT ..........', PFLAG( 7),
     & '......', FFLAG( 7),'.....',GFLAG( 7)
      WRITE(IU06,*) '8  NORMALISED WAVE STRESS ....', PFLAG( 8),
     & '......', FFLAG( 8),'.....',GFLAG( 8)
      WRITE(IU06,*) '9  MEAN SQUARE SLOPE..........', PFLAG( 9),
     & '......', FFLAG( 9),'.....',GFLAG( 9)
      WRITE(IU06,*) '10 WAVE MODEL WIND SPEED......', PFLAG(10),
     & '......', FFLAG(10),'.....',GFLAG(10)
      WRITE(IU06,*) '11 WIND SEA WAVE HEIGHT.......', PFLAG(11),
     & '......', FFLAG(11),'.....',GFLAG(11)
      WRITE(IU06,*) '12 SWELL WAVE HEIGHT..........', PFLAG(12),
     & '......', FFLAG(12),'.....',GFLAG(12)
      WRITE(IU06,*) '13 WIND SEA DIRECTION.........', PFLAG(13),
     & '......', FFLAG(13),'.....',GFLAG(13)
      WRITE(IU06,*) '14 SWELL WAVE DIRECTION.......', PFLAG(14),
     & '......', FFLAG(14),'.....',GFLAG(14)
      WRITE(IU06,*) '15 WIND SEA MEAN PERIOD (-1) .', PFLAG(15),
     & '......', FFLAG(15),'.....',GFLAG(15)
      WRITE(IU06,*) '16 SWELL MEAN PERIOD (-1) ....', PFLAG(16),
     & '......', FFLAG(16),'.....',GFLAG(16)
      WRITE(IU06,*) ' '
      WRITE(IU06,*) '22 ALTIMETER WAVE HEIGHT......', PFLAG(22),
     & '......', FFLAG(22),'.....',GFLAG(22)
      WRITE(IU06,*) '23 CORRECTED ALT WAVE HEIGHT..', PFLAG(23),
     & '......', FFLAG(23),'.....',GFLAG(23)
      WRITE(IU06,*) '24 ALTIMETER RANGE CORRECTION.', PFLAG(24),
     & '......', FFLAG(24),'.....',GFLAG(24)
      WRITE(IU06,*) ' '
      WRITE(IU06,*) '25 WAVE MEAN PERIOD (1)  .....', PFLAG(25),
     & '......', FFLAG(25),'.....',GFLAG(25)
      WRITE(IU06,*) '26 WAVE MEAN PERIOD (2)  .....', PFLAG(26),
     & '......', FFLAG(26),'.....',GFLAG(26)
      WRITE(IU06,*) '27 WAVE DIRECTIONAL SPREAD ...', PFLAG(27),
     & '......', FFLAG(27),'.....',GFLAG(27)
      WRITE(IU06,*) '28 WIND SEA MEAN PERIOD (1) ..', PFLAG(28),
     & '......', FFLAG(28),'.....',GFLAG(28)
      WRITE(IU06,*) '29 SWELL MEAN PERIOD (1) .....', PFLAG(29),
     & '......', FFLAG(29),'.....',GFLAG(29)
      WRITE(IU06,*) '30 WIND SEA MEAN PERIOD (2) ..', PFLAG(30),
     & '......', FFLAG(30),'.....',GFLAG(30)
      WRITE(IU06,*) '31 SWELL MEAN PERIOD (2) .....', PFLAG(31),
     & '......', FFLAG(31),'.....',GFLAG(31)
      WRITE(IU06,*) '32 WIND SEA DIRECTIONAL SPREAD', PFLAG(32),
     & '......', FFLAG(32),'.....',GFLAG(32)
      WRITE(IU06,*) '33 SWELL DIRECTIONAL SPREAD ..', PFLAG(33),
     & '......', FFLAG(33),'.....',GFLAG(33)
      WRITE(IU06,*) '34 WAVE SPECTRAL KURTOSIS.....', PFLAG(34),
     & '......', FFLAG(34),'.....',GFLAG(34)
      WRITE(IU06,*) '35 BENJAMIN-FEIR INDEX   .....', PFLAG(35),
     & '......', FFLAG(35),'.....',GFLAG(35)
      WRITE(IU06,*) '36 WAVE SPECTRAL PEAKEDNESS...', PFLAG(36),
     & '......', FFLAG(36),'.....',GFLAG(36)
      WRITE(IU06,*) '37 BATHYMETRY (an & fc step 0)', PFLAG(37),
     & '......', FFLAG(37),'.....',GFLAG(37)

      IF(FFLAG21) CALL OUTXT(IU06,IU21,IU31,IU23,IU33,IDELINS,IASSI,
     & NOUTT,'IU21','IU31','IU23','IU33')
      IF(FFLAG25) CALL OUTXT(IU06,IU25,IU35,IU27,IU37,IDELSPT,IASSI,
     & NOUTT,'IU25','IU35','IU27','IU37')
      IF(FFLAG26) CALL OUTXT(IU06,IU26,IU36,IU28,IU38,IDELSPS,IASSI,
     & NOUTT,'IU26','IU36','IU28','IU38')

      IF(ITEST.GT.0) THEN
      WRITE(IU06,*) '                        PRINTER   UNIT    GRIB'
      WRITE(IU06,*) '17 SPECTRA OF TOTAL SEA .....', PFLAG(17),
     & '......', FFLAG(17),'.....',GFLAG(17)
      WRITE(IU06,*) ' SEE INSTEAD OUTPUT OF RESTART SPECTRA IN GRIB !!!'
      WRITE(IU06,*) '  '

      WRITE(IU06,*) '                        PRINTER   UNIT    GRIB'
      WRITE(IU06,*) '18 SWELL SPECTRA ............', PFLAG(18),
     & '......', FFLAG(18),'.....',GFLAG(18)
      WRITE(IU06,*) ' GRIBBING NOT YET IMPLEMENTED !!!'
      WRITE(IU06,*) '  '
      ENDIF

      WRITE(IU06,*) '  '

      WRITE(IU06,*) ' ACCESS TO THE FIELD DATA BASE: '
     & ,'    F = DISABLED   T = ENABLED    ', LFDB
      IF (LFDB) THEN
        ILCFDBSF = LEN_TRIM(CFDBSF)
        IF(ILCFDBSF.GT.0) THEN
          CLFORM='("  ITS ROOT DIRECTORY IS ", A256)'
          WRITE (CLFORM(31:33), '(I3.3)') ILCFDBSF
          WRITE(IU06,CLFORM) CFDBSF(1:ILCFDBSF)
        ELSE
          WRITE(IU06,*)' **********************************************'
          WRITE(IU06,*)'   WARNING : FDB ROOT DIRECTORY NOT SPECIFIED ' 
          WRITE(IU06,*)' **********************************************'
        ENDIF
      ENDIF
      WRITE(IU06,*) '  '


      IF ( LFDB .AND. .NOT. (GFLAG20.OR.GFLAG21))
     &  THEN
        WRITE(IU06,*)' ************************************************'
        WRITE(IU06,*)' *                                              *'
        WRITE(IU06,*)' * ACCESS TO THE FIELD DATA BASE REQUIRES GRIB  *'
        WRITE(IU06,*)' * CODED DATA BUT NO GFLAG WAS SET TO TRUE      *'
        WRITE(IU06,*)' * THIS IS CHANGED AUTOMATICALLY AND THE        *'
        WRITE(IU06,*)' * STANDARD SET OF PARAMETERS IS ENABLED        *'
        WRITE(IU06,*)' * FOR PACKING                                  *'
        WRITE(IU06,*)' *                                              *'
        WRITE(IU06,*)' ************************************************'
        GFLAG( 1) = .TRUE.
        GFLAG( 2) = .TRUE.
        GFLAG( 3) = .TRUE.
        GFLAG( 6) = .TRUE.
        GFLAG( 7) = .TRUE.
        GFLAG(15) = .TRUE.
        GFLAG(16) = .TRUE.
        GFLAG(17) = .TRUE.
        GFLAG(18) = .TRUE.
        GFLAG(19) = .TRUE.
        GFLAG(20) = .TRUE.
        GFLAG(21) = .TRUE.
        GFLAG20 = .TRUE.
        GFLAG21 = .TRUE.
      ENDIF
      IF (LFDB) THEN
        WRITE(IU06,*) ' OUTPUT OF GRIB INTEGRATED PARAMETERS REDIRECTED'
     &   ,' TO THE FIELD DATA BASE'
        WRITE(IU06,*) '                    '
      ENDIF

      WRITE(IU06,*) '  '
      IF (NOUTS.GT.0.AND.IREST.EQ.1) THEN
        WRITE(IU06,*) ' NUMBER OF SPECTRA OUTPUT TIMES IS NOUTS = ',
     &   NOUTS
        WRITE(IU06,*) ' SPECTRA OUTPUT WILL BE PROCESSED AT:'
        WRITE(IU06,'(6(2X,A12))') (COUTS(I),I=1,NOUTS)
        WRITE(IU06,*) '  '

      ELSEIF (IREST.EQ.1) THEN
        IF(.NOT.LGRIBOUT) THEN
          WRITE(IU06,*) ' SPECTRA FILES WILLBE WRITTEN OUT TO DISK '
     &    ,'EVERY ...', IDELRES, ' SECONDS AND AT THE END OF THE RUN'
        ELSE
          WRITE(IU06,*) ' RESTART FILES WILL BE DISPOSED '
     &    ,'EVERY ...', IDELRES, ' SECONDS AND AT THE END OF THE RUN'
        ENDIF
        IF(CDATER.LT.CDATEE.AND..NOT.LGRIBOUT) WRITE(IU06,*)
     &   ' !! HOWEVER BOTH RESTART FILES WILL ONLY BE SAVED',
     &   ' AT ...', CDATER 
        IF(CDATES.LT.CDATEE) WRITE(IU06,*)
     &   ' BUT SPECTRA FILES ALONE  WILL BE SAVED UNTIL '
     &   ,'...', CDATES 
      ELSE
        WRITE(IU06,*) ' SPECTRA FILES WILL NOT BE WRITTEN OUT TO DISK'
      ENDIF

      WRITE(IU06,*) '  '
      IF(LPBIOIN.AND..NOT.LGRIBIN) THEN
        WRITE (IU06,*)
     &   ' PBIO SOFTWARE IS USED TO READ SELECTED INPUT FILES' 
      ENDIF
      IF(LNOCDIN.AND.LGRIBIN) THEN
        WRITE(IU06,*) '  '
        WRITE (IU06,*) 
     &   ' NO DRAG COEFFICIENT FIELD IS PROVIDED AS INPUT'
        WRITE (IU06,*) ' THE FIELD WILL BE INITIALISED BY TAKING' 
        WRITE (IU06,*) ' ZERO WAVE STRESS (TAUW)' 
      ENDIF
      IF(LGRIBIN.AND.ITEST.GT.0) THEN
        WRITE(IU06,*) '  '
        WRITE (IU06,*)
     &   ' GRIB SPECTRA FIELD ARE USED AS INPUT' 
      ENDIF
      IF(LPBIOOUT.AND..NOT.LGRIBOUT) THEN
        WRITE(IU06,*) ' SP '
        WRITE (IU06,*)
     &   ' PBIO SOFTWARE IS USED TO WRITE TO SELECTED OUTPUT FILES' 
      ENDIF
      WRITE(IU06,*) '  '
      IF(LFDBIOIN) THEN
        WRITE (IU06,*)
     &   ' FDB SOFTWARE IS USED TO READ SELECTED INPUT FILES' 
      ELSE
        IF(ITEST.GT.0) WRITE (IU06,*)
     &   ' INPUT FILES ARE READ FROM DISK' 
      ENDIF
      WRITE(IU06,*) '  '
      IF(LFDBIOOUT) THEN
        WRITE (IU06,*)
     &   ' FDB SOFTWARE IS USED TO WRITE OUTPUT SPECTRA FILES' 
      ELSE
        WRITE (IU06,*)
     &   ' OUTPUT SPECTRA FILES ARE WRITTEN OUT TO DISK' 
      ENDIF
!     CFDB2DSP is reset to null since the FDB will be initialised
!     with CFDBSF, provided it is used (ie grib output of integrated
!     parameters (see fld2fdb).
      IF(LFDB.AND.(GFLAG20.OR.GFLAG21)) CFDB2DSP  = "" 
      IF (LFDBIOIN.OR.LFDBIOOUT) THEN
        CLFORM='("  ITS ROOT DIRECTORY IS ", A256)'
        ILCFDB2DSP = LEN_TRIM(CFDB2DSP)
        IF ( ILCFDB2DSP .EQ. 0 ) THEN
          ILCFDBSF = LEN_TRIM(CFDBSF)
          IF(ILCFDBSF.GT.0) THEN
            WRITE (CLFORM(31:33), '(I3.3)') ILCFDBSF
            WRITE(IU06,CLFORM) CFDBSF(1:ILCFDBSF)
          ENDIF
        ELSE
          WRITE (CLFORM(31:33), '(I3.3)') ILCFDB2DSP
          WRITE(IU06,CLFORM) CFDB2DSP(1:ILCFDB2DSP)
        ENDIF
      ENDIF
      WRITE(IU06,*) '  '
      CALL FLUSH(IU06)

      CALL WSTREAM_STRG(ISTREAM,CSTREAM,NENSFNB,NTOTENS,MARSFCTYPE,
     &                  KSTREAM, LASTREAM)

      WRITE(IU06,'("  HARD DRIVE PATH NAME : ",/,5X,A70)') CPATH
      WRITE(IU06,*) '  '
      WRITE(IU06,'("  CURRENT RUN:")')
      WRITE(IU06,'("  GRIB TABLE : ", I4)') NLOCGRB
      WRITE(IU06,'("  STREAM ....: ", A4)') CSTREAM 
      WRITE(IU06,'("  CLASS......: ", A4)') YCLASS
      WRITE(IU06,'("  EXPERIMENT.: ", A4)') YEXPVER

      IF(LGRHDIFS) THEN
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' !!  KSEC1 HEADERS ARE SET BY IFS !! ' 
        WRITE(IU06,*) ' '
      ENDIF

      IF ( ISTREAM .EQ. 1082 .OR. ISTREAM .EQ. 1095 .OR.
     &     ISTREAM .EQ. 1203 .OR. ISTREAM .EQ. 1204 ) THEN
        WRITE(IU06,*) '  '
        IF(ISTREAM.EQ.1095 .OR. ISTREAM.EQ.1203 ) THEN
          WRITE(IU06,'("  MONTHLY FORECAST RUN : ")')
        ELSEIF(ISTREAM.EQ.1204 ) THEN
          WRITE(IU06,'("  MONTHLY FORECAST HINDCAST RUN : ")')
          WRITE(IU06,'("  WITH REFERENCE DATE:      ", I8  )') IREFDATE
        ELSE
          WRITE(IU06,'("  SEASONAL FORECAST RUN : ")')
        ENDIF
        WRITE(IU06,'("  ************************************ ")')
        WRITE(IU06,'("  ENSEMBLE NUMBER:         ", I4  )') NENSFNB
        WRITE(IU06,'("  TOTAL NUMBER OF ENSEMBLE:", I4,/)') NTOTENS
        WRITE(IU06,'("  SYSTEM NUMBER:           ", I4  )') NSYSNB
        WRITE(IU06,'("  METHOD NUMBER:           ", I4  )') NMETNB 
        WRITE(IU06,*) '  '
      ELSE IF ( ISTREAM .EQ. 1083 ) THEN
        WRITE(IU06,'("  MULTI ANALYSIS FORECAST HINDCAST RUN : ")')
        WRITE(IU06,'("  ************************************ ")')
        WRITE(IU06,'("  ENSEMBLE NUMBER:         ", I4  )') NENSFNB
        WRITE(IU06,'("  TOTAL NUMBER OF ENSEMBLE:", I4,/)') NTOTENS
      ELSE IF ( ISTREAM .EQ. 1084 ) THEN
        WRITE(IU06,'("  ENSEMBLE FORECAST HINDCAST RUN : ")')
        WRITE(IU06,'("  WITH REFERENCE DATE:      ", I8  )') IREFDATE
        WRITE(IU06,'("  ENSEMBLE FORECAST RUN : ")')
        WRITE(IU06,'("  ************************************ ")')
        WRITE(IU06,'("  ENSEMBLE NUMBER:         ", I4  )') NENSFNB
        WRITE(IU06,'("  TOTAL NUMBER OF ENSEMBLE:", I4,/)') NTOTENS
      ELSE IF ( ISTREAM .EQ. 1085 ) THEN
        WRITE(IU06,'("  FORECAST HINDCAST RUN : ")')
        WRITE(IU06,'("  WITH REFERENCE DATE:      ", I8  )') IREFDATE
      ELSE IF ( ISTREAM .EQ. 1078 ) THEN
        WRITE(IU06,'("  NEW ENSEMBLE FORECAST HINDCAST RUN OVERLAP: ")')
        WRITE(IU06,'("  WITH REFERENCE DATE:      ", I8  )') IREFDATE
        WRITE(IU06,'("  ENSEMBLE FORECAST RUN : ")')
        WRITE(IU06,'("  ************************************ ")')
        WRITE(IU06,'("  ENSEMBLE NUMBER:         ", I4  )') NENSFNB
        WRITE(IU06,'("  TOTAL NUMBER OF ENSEMBLE:", I4,/)') NTOTENS
      ELSE IF ( ISTREAM .EQ. 1079 ) THEN
        WRITE(IU06,'("  NEW ENSEMBLE FORECAST HINDCAST RUN : ")')
        WRITE(IU06,'("  WITH REFERENCE DATE:      ", I8  )') IREFDATE
        WRITE(IU06,'("  ENSEMBLE FORECAST RUN : ")')
        WRITE(IU06,'("  ************************************ ")')
        WRITE(IU06,'("  ENSEMBLE NUMBER:         ", I4  )') NENSFNB
        WRITE(IU06,'("  TOTAL NUMBER OF ENSEMBLE:", I4,/)') NTOTENS
      ELSE
        IF (NENSFNB.NE.0.OR.NTOTENS.NE.0) THEN
          WRITE(IU06,'("  ENSEMBLE FORECAST RUN : ")')
          WRITE(IU06,'("  ************************************ ")')
          WRITE(IU06,'("  ENSEMBLE NUMBER:         ", I4  )') NENSFNB
          WRITE(IU06,'("  TOTAL NUMBER OF ENSEMBLE:", I4,/)') NTOTENS
        ENDIF
      ENDIF

      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' TEST OUTPUT LEVEL IS .............. ITEST = ',
     & ITEST
      WRITE(IU06,*) ' TEST OUTPUT OF IN BLOCK LOOPS UPTO ITESTB = ',
     & ITESTB
      WRITE(IU06,*) '  '

      IF (.NOT. LWCOU .AND.
     & (  (LFDB .OR. LFDBIOIN .OR. LFDBIOOUT ))) THEN
        WRITE(IU06,*) ' ACCESS TO FIELD DATA BASE REQUESTED '
        WRITE(IU06,*) ' DATA BASE WILL BE INITIALIZED '
        WRITE(IU06,*) ' '
        INITFDB_STAT = IINITFDBSUBS()
        IF (INITFDB_STAT .NE. 0) THEN
          WRITE(IU06,'("****************************************")')
          WRITE(IU06,'("*                                      *")')
          WRITE(IU06,'("*    FATAL ERROR IN SUB. USERIN        *")')
          WRITE(IU06,'("*    ==========================        *")')
          WRITE(IU06,'("*    INITIALIZATIO OF THE              *")')
          WRITE(IU06,'("*    FIELD DATA BASE                   *")')
          WRITE(IU06,'("*     F A I L E D .                    *")')
          WRITE(IU06,'("*                                      *")')
          WRITE(IU06,'("*  ERROR CODE (ISTAT) =",i7, "         *")')
     &     INITFDB_STAT
          WRITE(IU06,'("*                                      *")')
          WRITE(IU06,'("****************************************")')
          CALL ABORT1
        ENDIF
      ENDIF

! ----------------------------------------------------------------------

!*    2. CHECK INTEGER RATIOS BETWEEN TIMESTEPS.
!        ---------------------------------------

      LERROR = .FALSE.

!*    2.1 WIND OUTPUT AND PROPAGATION TIME STEP.
!         --------------------------------------

      IF ((IDELWO.LT.IDELPRO.AND.MOD(IDELPRO,IDELWO).NE.0) .OR.
     &    (IDELWO.GE.IDELPRO.AND.MOD(IDELWO,IDELPRO).NE.0)) THEN
        WRITE(IU06,*) '*******************************************'
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
        WRITE(IU06,*) '*    ==========================           *'
        WRITE(IU06,*) '* WIND OUTPUT TIMSTEP AND PROPAGATION     *'
        WRITE(IU06,*) '* TIME STEP DO NOT HAVE INTEGER RATIO.    *'
        WRITE(IU06,*) '* WIND OUTPUT TIMSTEP    IDELWO = ', IDELWO
        WRITE(IU06,*) '* PROPAGATION TIME STEP IDELPRO = ', IDELPRO
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*******************************************'
        LERROR = .TRUE.
      ENDIF

!*    2.2 SOURCE FUNCTION AND PROPAGATION TIMESTEP.
!         -----------------------------------------

      IF (IDELT.GT.IDELPRO.OR.MOD(IDELPRO,IDELT).NE.0) THEN
        WRITE(IU06,*) '*******************************************'
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
        WRITE(IU06,*) '*    ==========================           *'
        WRITE(IU06,*) '* SOURCE FUNCTION  AND PROPAGATION        *'
        WRITE(IU06,*) '* TIME STEP DO NOT HAVE INTEGER RATIO OR  *'
        WRITE(IU06,*) '* SOURCE FUNCTION TIMESTEP IS GREATER     *'
        WRITE(IU06,*) '* THEN PROPAGATION TIMESTEP.              *'
        WRITE(IU06,*) '* SOURCE FUNCTION TIMESTEP IDELT = ', IDELT
        WRITE(IU06,*) '* PROPAGATION TIMESTEP   IDELPRO = ', IDELPRO
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*******************************************'
        LERROR = .TRUE.
      ENDIF

!*    2.3 SOURCE FUNCTION AND WIND OUTPUT TIMESTEP.
!         -----------------------------------------

      IF (MOD(IDELWO,IDELT).NE.0) THEN
        WRITE(IU06,*) '*******************************************'
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
        WRITE(IU06,*) '*    ==========================           *'
        WRITE(IU06,*) '* SOURCE FUNCTION  AND WIND OUTPUT        *'
        WRITE(IU06,*) '* TIME STEP DO NOT HAVE INTEGER RATIO OR  *'
        WRITE(IU06,*) '* SOURCE FUNCTION TIMESTEP IDELT = ', IDELT
        WRITE(IU06,*) '* WIND OUTPUT TIMESTEP    IDELWO = ', IDELWO
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*******************************************'
        LERROR = .TRUE.
      ENDIF

!*    2.4 WIND INPUT AND WIND OUTPUT TIMESTEP.
!         ------------------------------------

      IF (IDELWO.GT.IDELWI) THEN
        WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*) '+                                         +'
        WRITE(IU06,*) '+   WARNING ERROR IN SUB. USERIN          +'
        WRITE(IU06,*) '+   ============================          +'
        WRITE(IU06,*) '+ WIND INPUT TIME STEP IS LESS THAN       +'
        WRITE(IU06,*) '+ WIND OUTPUT STEP                        +'
        WRITE(IU06,*) '+ WIND INPUT TIMESTEP   IDELWI = ', IDELWI
        WRITE(IU06,*) '+ WIND OUTPUT TIMESTEP  IDELWO = ', IDELWO
        WRITE(IU06,*) '+                                         +'
        WRITE(IU06,*) '+ WIND INPUT CHANGED TO WIND OUTPUT       +'
        WRITE(IU06,*) '+ MODEL WILL USE A NEW WIND FIELD EVERY   +'
        WRITE(IU06,*) '+ WIND OUTPUT TIME STEP AND IGNORE FIELDS +'
        WRITE(IU06,*) '+ IN BETWEEN.                             +'
        WRITE(IU06,*) '+                                         +'
        WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++'
        IDELWI = IDELWO
      ENDIF
      IF ((IDELWO.LE.IDELWI.AND.MOD(IDELWI,IDELWO).NE.0)) THEN
        WRITE(IU06,*) '*******************************************'
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
        WRITE(IU06,*) '*    ==========================           *'
        WRITE(IU06,*) '* WIND INPUT  AND WIND OUTPUT             *'
        WRITE(IU06,*) '* TIME STEP DO NOT HAVE INTEGER RATIO OR  *'
        WRITE(IU06,*) '* WIND INPUT TIMESTEP   IDELWI = ', IDELWI
        WRITE(IU06,*) '* WIND OUTPUT TIMESTEP  IDELWO = ', IDELWO
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*******************************************'
        LERROR = .TRUE.
      ENDIF

!*    2.5 FILE DISPOSE TIMESTEP.
!         ----------------------

      IF (MOD(IDELRES,IDELPRO).NE.0 .AND.
     &    MOD(IDELRES,IDELWI).NE.0) THEN
        WRITE(IU06,*) '*******************************************'
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
        WRITE(IU06,*) '*    ==========================           *'
        WRITE(IU06,*) '* NEW OUTPUT FILES ARE REQUESTED EVERY    *'
        WRITE(IU06,*) '*    IDELRES = ',IDELRES,' SECONDS'
        WRITE(IU06,*) '* IDELRES MUST BE MULTIPLES OF            *'
        WRITE(IU06,*) '* THE WIND INPUT TIMESTEP   IDELWI = ', IDELWI
        WRITE(IU06,*) '* THE PROPAGATION TIMESTEP IDELPRO = ',IDELPRO
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*******************************************'
        LERROR = .TRUE.
      ENDIF

!*    2.5 OUTPUT OPTION.
!         --------------

      IF (NOUTT.GT.0) THEN
        DO J=1,NOUTT
          CALL DIFDATE (CDATEA, COUTT(J), ISHIFT)
          IF (ISHIFT.LE.0 .OR. MOD(ISHIFT,IDELPRO).NE.0) THEN
            WRITE(IU06,*) '++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*) '+                                      +'
            WRITE(IU06,*) '+    WARNING ERROR IN SUB. USERIN      +'
            WRITE(IU06,*) '+    ============================      +'
            WRITE(IU06,*) '+ OUTPUT DATE IS NOT AT THE END OF A   +'
            WRITE(IU06,*) '+ PROPAGATION TIMESTEP.                +'
            WRITE(IU06,*) '+ DATE IS : ', COUTT(J)
            WRITE(IU06,*) '+ PROGRAM WILL IGNORE THIS OUTPUT TIME +'
            WRITE(IU06,*) '+                                      +'
            WRITE(IU06,*) '++++++++++++++++++++++++++++++++++++++++'
          ENDIF
        ENDDO
      ELSE
        IF ((FFLAG20.OR.PFLAG20.OR.GFLAG20) .AND. IDELINT.EQ.0) THEN
          WRITE(IU06,*) '*******************************************'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
          WRITE(IU06,*) '*    ==========================           *'
          WRITE(IU06,*) '* OUTPUT OF INTEGRATED DATA (TOTAL SEA)   *'
          WRITE(IU06,*) '* IS REQUESTED.                           *'
          WRITE(IU06,*) '* OUTPUT TIME STEP IDELINT                *'
          WRITE(IU06,*) '* HAS TO BE SPECIFIED IN NAMELIST         *'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*******************************************'
          LERROR = .TRUE.
        ENDIF
        IF ((FFLAG20.OR.PFLAG20.OR.GFLAG20) .AND.
     &   MOD(IDELINT,IDELPRO).NE.0) THEN
          WRITE(IU06,*) '*******************************************'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
          WRITE(IU06,*) '*    ==========================           *'
          WRITE(IU06,*) '* OUTPUT OF INTEGRATED DATA (TOTAL SEA)   *'
          WRITE(IU06,*) '* IS REQUESTED.                           *'
          WRITE(IU06,*) '* OUTPUT TIME STEP HAS TO BE A MULTIPLE   *'
          WRITE(IU06,*) '* OF THE PROPAGATION TIME STEP.           *'
          WRITE(IU06,*) '* OUTPUT TIME STEP IS      IDELINT = ',
     &     IDELINT
          WRITE(IU06,*) '* PROPAGATION TIME STEP IS IDELPRO = ',
     &     IDELPRO
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*******************************************'
          LERROR = .TRUE.
        ENDIF
        IF ((FFLAG21.OR.PFLAG21.OR.GFLAG21) .AND. IDELINS.EQ.0 ) THEN
          WRITE(IU06,*) '*******************************************'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
          WRITE(IU06,*) '*    ==========================           *'
          WRITE(IU06,*) '* OUTPUT OF INTEGRATED DATA (SEA + SWELL) *'
          WRITE(IU06,*) '* IS REQUESTED.                           *'
          WRITE(IU06,*) '* OUTPUT TIME STEP IDELINS                *'
          WRITE(IU06,*) '* HAS TO BE SPECIFIED IN NAMELIST         *'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*******************************************'
          LERROR = .TRUE.
        ENDIF
        IF ((FFLAG21.OR.PFLAG21.OR.GFLAG21) .AND.
     &   MOD(IDELINS,IDELPRO).NE.0) THEN
          WRITE(IU06,*) '*******************************************'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
          WRITE(IU06,*) '*    ==========================           *'
          WRITE(IU06,*) '* OUTPUT OF INTEGRATED DATA (SEA + SWELL) *'
          WRITE(IU06,*) '* IS REQUESTED.                           *'
          WRITE(IU06,*) '* OUTPUT TIME STEP HAS TO BE A MULTIPLE   *'
          WRITE(IU06,*) '* OF THE PROPAGATION TIME STEP.           *'
          WRITE(IU06,*) '* OUTPUT TIME STEP IS      IDELINS = ',
     &     IDELINS
          WRITE(IU06,*) '* PROPAGATION TIME STEP IS IDELPRO = ',
     &     IDELPRO
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*******************************************'
          LERROR = .TRUE.
        ENDIF
        IF ((FFLAG25.OR.PFLAG25.OR.GFLAG25) .AND. IDELSPT.EQ.0) THEN
          WRITE(IU06,*) '*******************************************'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
          WRITE(IU06,*) '*    ==========================           *'
          WRITE(IU06,*) '* OUTPUT OF SPECTRA (TOTAL SEA)           *'
          WRITE(IU06,*) '* IS REQUESTED.                           *'
          WRITE(IU06,*) '* OUTPUT TIME STEP IDELSPT                *'
          WRITE(IU06,*) '* HAS TO BE SPECIFIED IN NAMELIST         *'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*******************************************'
          LERROR = .TRUE.
        ENDIF
        IF ((FFLAG25.OR.PFLAG25.OR.GFLAG25) .AND.
     &   MOD(IDELSPT,IDELPRO).NE.0) THEN
          WRITE(IU06,*) '*******************************************'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
          WRITE(IU06,*) '*    ==========================           *'
          WRITE(IU06,*) '* OUTPUT OF SPECTRA (TOTAL SEA)           *'
          WRITE(IU06,*) '* IS REQUESTED.                           *'
          WRITE(IU06,*) '* OUTPUT TIME STEP HAS TO BE A MULTIPLE   *'
          WRITE(IU06,*) '* OF THE PROPAGATION TIME STEP.           *'
          WRITE(IU06,*) '* OUTPUT TIME STEP IS      IDELSPT = ',
     &     IDELSPT
          WRITE(IU06,*) '* PROPAGATION TIME STEP IS IDELPRO = ',
     &     IDELPRO
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*******************************************'
          LERROR = .TRUE.
        ENDIF
        IF ((FFLAG26.OR.PFLAG26.OR.GFLAG26) .AND. IDELSPS.EQ.0) THEN
          WRITE(IU06,*) '*******************************************'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
          WRITE(IU06,*) '*    ==========================           *'
          WRITE(IU06,*) '* OUTPUT OF SPECTRA (SEA + SWELL)         *'
          WRITE(IU06,*) '* IS REQUESTED.                           *'
          WRITE(IU06,*) '* OUTPUT TIME STEP IDELSPS                *'
          WRITE(IU06,*) '* HAS TO BE SPECIFIED IN NAMELIST         *'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*******************************************'
          LERROR = .TRUE.
        ENDIF
        IF ((FFLAG26.OR.PFLAG26.OR.GFLAG26) .AND.
     &   MOD(IDELSPS,IDELPRO).NE.0) THEN
          WRITE(IU06,*) '*******************************************'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
          WRITE(IU06,*) '*    ==========================           *'
          WRITE(IU06,*) '* OUTPUT OF SPECTRA (SEA + SWELL)         *'
          WRITE(IU06,*) '* IS REQUESTED.                           *'
          WRITE(IU06,*) '* OUTPUT TIME STEP HAS TO BE A MULTIPLE   *'
          WRITE(IU06,*) '* OF THE PROPAGATION TIME STEP.           *'
          WRITE(IU06,*) '* OUTPUT TIME STEP IS      IDELSPS = ',
     &     IDELSPS
          WRITE(IU06,*) '* PROPAGATION TIME STEP IS IDELPRO = ',
     &     IDELPRO
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*******************************************'
          LERROR = .TRUE.
        ENDIF
      ENDIF

      IF (NOUTS.GT.0) THEN
        DO J=1,NOUTS
          CALL DIFDATE (CDATEA, COUTS(J), ISHIFT)
          IF (ISHIFT.LE.0 .OR. MOD(ISHIFT,IDELPRO).NE.0) THEN
            WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*) '+                                       +'
            WRITE(IU06,*) '+    WARNING ERROR IN SUB. USERIN       +'
            WRITE(IU06,*) '+    ============================       +'
            WRITE(IU06,*) '+ SPECTRA OUTPUT DATE IS NOT AT THE END +'
            WRITE(IU06,*) '+ OF A PROPAGATION TIMESTEP.            +'
            WRITE(IU06,*) '+ DATE IS : ', COUTS(J)
            WRITE(IU06,*) '+ PROGRAM WILL IGNORE THIS OUTPUT TIME  +'
            WRITE(IU06,*) '+                                       +'
            WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++'
          ENDIF
        ENDDO
      ENDIF

      IF (NASS.GT.0.AND.IASSI.EQ.1) THEN
        DO J=1,NASS
          CALL DIFDATE (CDATEA, CASS(J), ISHIFT)
          IF (ISHIFT.LE.0 .OR. MOD(ISHIFT,IDELPRO).NE.0) THEN
            WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*) '+                                       +'
            WRITE(IU06,*) '+    WARNING ERROR IN SUB. USERIN       +'
            WRITE(IU06,*) '+    ============================       +'
            WRITE(IU06,*) '+ ASSIMILATION DATE IS NOT AT THE END   +'
            WRITE(IU06,*) '+ OF A PROPAGATION TIMESTEP.            +'
            WRITE(IU06,*) '+ DATE IS : ', CASS(J)
            WRITE(IU06,*) '+ PROGRAM WILL ABORT ' 
            WRITE(IU06,*) '+                                       +'
            WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++'
            CALL ABORT1
          ENDIF
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    3. ERROR CHECK.
!        ------------

      IF (LERROR) THEN
        WRITE(IU06,*) '*******************************************'
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*    FATAL ERROR(S) IN SUB. USERIN        *'
        WRITE(IU06,*) '*    =============================        *'
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '* CORRECT USER INPUT AS INDICATED ABOVE   *'
        WRITE(IU06,*) '* AND TRY AGAIN!!!!!!!!!!!!!!!!!!!!!!!!   *'
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
        WRITE(IU06,*) '* ---------------   --------------        *'
        WRITE(IU06,*) '*******************************************'
        CALL ABORT1
      ELSE
        RETURN
      ENDIF

! ----------------------------------------------------------------------

      END SUBROUTINE USERIN
