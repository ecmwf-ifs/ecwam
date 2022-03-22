
      SUBROUTINE MPUSERIN

! ----------------------------------------------------------------------

!**** *MPUSERIN* - ROUTINE TO READ NAMELIST INPUTS.

!     J. BIDLOT      ECMWF          JULY    1996
!     B. HANSEN      ECMWF          APRIL   1997
!     B. HANSEN      ECMWF          JANUARY 1998 NAMELIST INPUT.
!     S. ABDALLA     ECMWF          OCTOBER 2001 LGUST & LADEN ADDED

!*    PURPOSE.
!     --------
!     READ LMESSPASS AND LFDB FROM USER INPUT SKIPPING ALL THE OTHER
!     INPUT PARAMETER WHICH WILL BE READ IN WITH USERIN.

!**   INTERFACE.
!     ----------

!      *CALL* *MPUSERIN

!     METHOD.
!     -------

!     SEE USERIN

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWALTAS , ONLY : NUMALT   ,IBUFRSAT  ,ALTSDTHRSH,ALTBGTHRSH, &
     &            HSALTCUT, LALTGRDOUT, LALTPAS, LALTPASSIV,            &
     &            XKAPPA2  ,HSCOEFCOR,HSCONSCOR ,LALTCOR   ,LALTLRGR,   &
     &            LODBRALT ,CSATNAME
      USE YOWCOUP  , ONLY : LWCOU    ,KCOUSTEP  ,LWFLUX ,LWVFLX_SNL,    &
     &            LWCOUNORMS, LLNORMIFS2WAM,LLNORMWAM2IFS,LLNORMWAMOUT, &
     &            LLNORMWAMOUT_GLOBAL,                                  &
     &            LWNEMOCOU, LWNEMOCOUSEND, LWNEMOCOURECV,              &
     &            LWNEMOCOUDEBUG, LWNEMOCOUCIC, LWNEMOCOUCIT,           &
     &            LWNEMOCOUCUR,                                         &
     &            LWNEMOCOUSTK,  LWNEMOCOUSTRN, LWNEMOTAUOC, NEMOFRCO,  &
     &            LLCAPCHNK
      USE YOWCOUT  , ONLY : COUTT    ,COUTS    ,CASS     ,FFLAG    ,    &
     &            FFLAG20  ,GFLAG    ,                                  &
     &            GFLAG20  ,NFLAG    ,                                  &
     &            NFLAGALL ,UFLAG    ,LFDB     ,NOUTT    ,NOUTS    ,    &
     &            NASS     ,JPPFLAG  ,                                  &
     &            IRWDIR   , IRCD    ,IRU10    ,                        &
     &            LRSTPARALW,LRSTPARALR,LRSTINFDAT,                     &
     &            NTRAIN   ,                                            &
     &            IPFGTBL  ,                                            &
     &            LLOUTERS ,                                            &
     &            LWAMANOUT,                                            &
     &            NWRTOUTWAM,                                           &
     &            LSECONDORDER,                                         &
     &            LWAM_USE_IO_SERV
      USE YOWCPBO  , ONLY : GBOUNC_MAX, IBOUNC ,CBCPREF
      USE YOWCURR  , ONLY : IDELCUR  ,CDATECURA, LLCFLCUROFF
      USE YOWFPBO  , ONLY : IBOUNF
      USE YOWGRIBHD, ONLY : LGRHDIFS ,LNEWLVTP ,IMDLGRBID_G, IMDLGRBID_M
      USE YOWGRIB_HANDLES , ONLY : NGRIB_HANDLE_IFS
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE ,LCIWABR  ,              &
     &            CITHRSH  ,CIBLOCK  ,LICETH   ,                        &
     &            CITHRSH_SAT, CITHRSH_TAIL    ,CDICWA
      USE YOWMESPAS, ONLY : LMESSPASS,                                  &
     &            LFDBIOOUT,LGRIBIN  ,LGRIBOUT ,LNOCDIN
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : SWAMPWIND,SWAMPWIND2,DTNEWWIND,LTURN90 ,    &
     &            SWAMPCIFR,SWAMPCITH,LWDINTS  ,LL1D     ,CLDOMAIN
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDATER   ,CDATES   ,    &
     &            IDELPRO  ,IDELT    ,IDELWI   ,                        &
     &            IDELWO   ,IDELALT  ,IREST    ,IDELRES  ,IDELINT  ,    &
     &            IDELBC   ,                                            &
     &            IDELINS  ,IDELSPT  ,IDELSPS  ,ICASE    ,ISHALLO  ,    &
     &            ISNONLIN ,                                            &
     &            IDAMPING ,                                            &
     &            LBIWBK   ,                                            &
     &            IREFRA   ,IPROPAGS ,IASSI    ,                        &
     &            NENSFNB  ,NTOTENS  ,NSYSNB   ,NMETNB   ,CDATEA   ,    &
     &            YCLASS   ,YEXPVER  ,L4VTYPE  ,LFRSTFLD ,LALTAS   ,    &
     &            LSARAS   ,LSARINV  ,ISTREAM  ,NLOCGRB  ,NCONSENSUS,   &
     &            NDWD     ,NMFR     ,NNCEP    ,NUKM     ,IREFDATE ,    &
     &            LGUST    ,LADEN    ,NPROMA_WAM,LSUBGRID ,LLSOURCE ,   &
     &            LNSESTART,                                            &
     &            LSMSSIG_WAM,CMETER ,CEVENT   ,                        &
     &            LRELWIND ,                                            &
     &            IDELWI_LST, IDELWO_LST, CDTW_LST, NDELW_LST
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWTEXT  , ONLY : LRESTARTED,ICPLEN   ,USERID   ,RUNID    ,   &
     &            PATH     ,CPATH    ,CWI
      USE YOWUNPOOL, ONLY : LLUNSTR  ,LPREPROC, LVECTOR, IVECTOR
      USE UNSTRUCT_BOUND , ONLY : LBCWA
      USE UNWAM    , ONLY : USE_DIRECT_WIND_FILE
      USE UNWAM    , ONLY : LIMPLICIT, JGS_DIFF_SOLVERTHR,              &
     &            SOURCE_IMPL, WAE_SOLVERTHR,                           &
     &            LNONL, BLOCK_GAUSS_SEIDEL,                            &
     &            LLIMT, L_SOLVER_NORM, LCHKCONV
      USE YOWWAMI  , ONLY : CBEGDT   ,CENDDT   ,CBPLTDT  ,CEPLTDT  ,    &
     &            CLSPDT   ,CRSTDT   ,IANALPD  ,IFOREPD  ,IDELWIN  ,    &
     &            IASSIM   ,NFCST    ,ISTAT
      USE YOWWIND  , ONLY : CWDFILE  ,LLWSWAVE ,LLWDWAVE ,RWFAC
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "mpcrtbl.intfb.h"
#include "wposnam.intfb.h"

      INTEGER(KIND=JWIM) :: IU05
      INTEGER(KIND=JWIM) :: ISAT, IC, II
      INTEGER(KIND=JWIM) :: IWAM_GET_UNIT

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      CHARACTER(LEN=14), PARAMETER :: ZERO = ' '
      CHARACTER(LEN=1) :: CLMTSU(4), CLOTSU(8)
      CHARACTER(LEN=2) :: C2
      CHARACTER(LEN=70) :: CLHEADER
      CHARACTER(LEN=72) :: NAMELIST_FILENAME 

      LOGICAL :: LLEOF, LLOPENED

! ----------------------------------------------------------------------

      NAMELIST /NALINE/ CLHEADER,                                       &
     &   CBPLTDT, CEPLTDT, CDATEF,                                      &
     &   IDELPRO, IDELT, IDELWO, IDELWI, CLMTSU, IDELALT,               &
     &   IDELINT, IDELINS, IDELSPT, IDELSPS, IDELRES,                   &
     &   IDELCUR, CDATECURA,                                            &
     &   LLCFLCUROFF,                                                   &
     &   CLOTSU, CDATER, CDATES,                                        &
     &   FFLAG,  GFLAG, NFLAG,                                          &
     &   LLOUTERS,                                                      &
     &   LFDB, LGRIBIN, LGRIBOUT, LFDBIOOUT,                            &
     &   LRSTPARALW, LRSTPARALR, LRSTINFDAT,                            &
     &   LWAMANOUT,                                                     &
     &   NWRTOUTWAM,                                                    &
     &   LSECONDORDER,                                                  &
     &   ICASE, ISHALLO, ITEST, ITESTB, IREST, IASSI,                   &
     &   IPROPAGS,                                                      &
     &   IREFRA,                                                        &
     &   ISNONLIN,                                                      &
     &   IDAMPING,                                                      &
     &   LBIWBK  ,                                                      &
     &   LMASKICE,                                                      &
     &   IBOUNC, IBOUNF,                                                &
     &   IDELBC, CBCPREF,                                               &
     &   USERID, RUNID,  PATH, YCLASS, YEXPVER, CPATH,                  &
     &   IMDLGRBID_G, IMDLGRBID_M,                                      &
     &   NENSFNB, NTOTENS, NSYSNB, NMETNB,                              &
     &   LMESSPASS, LWCOU, LNOCDIN, LODBRALT,                           &
     &   LALTCOR, L4VTYPE, LFRSTFLD, LALTAS, LSARAS, LSARINV, XKAPPA2,  &
     &   IBUFRSAT, CSATNAME,                                            &
     &   SWAMPWIND, SWAMPWIND2, SWAMPCIFR, SWAMPCITH,                   &
     &   DTNEWWIND, LTURN90,                                            &
     &   LALTLRGR, HSCOEFCOR, HSCONSCOR,ALTSDTHRSH,ALTBGTHRSH,HSALTCUT, &
     &   ISTREAM, NLOCGRB, IREFDATE,                                    &
     &   NCONSENSUS, NDWD, NMFR, NNCEP, NUKM,                           &
     &   LGUST, LADEN, LRELWIND, LALTGRDOUT, LSUBGRID, LALTPAS,         &
     &   LLSOURCE,                                                      &
     &   LNSESTART,                                                     &
     &   LLUNSTR, LPREPROC, LVECTOR, IVECTOR,                           &
     &   WAE_SOLVERTHR, JGS_DIFF_SOLVERTHR,                             &
     &   USE_DIRECT_WIND_FILE,                                          &
     &   LIMPLICIT,                                                     &
     &   SOURCE_IMPL,                                                   &
     &   LNONL, BLOCK_GAUSS_SEIDEL,                                     &
     &   LLIMT, L_SOLVER_NORM, LCHKCONV,                                &
     &   LBCWA,                                                         &
     &   LSMSSIG_WAM,CMETER,CEVENT,                                     &
     &   LLWSWAVE, LLWDWAVE,                                            &
     &   NPROMA_WAM, LL1D, LGRHDIFS ,LNEWLVTP,                          &
     &   LWCOUNORMS, LLNORMIFS2WAM, LLNORMWAM2IFS, LLNORMWAMOUT,        &
     &   LLNORMWAMOUT_GLOBAL,                                           &
     &   LICERUN, LCIWABR, LICETH,                                      &
     &   LWVFLX_SNL,                                                    &
     &   LWNEMOCOU, NEMOFRCO,                                           &
     &   LWNEMOCOUSEND, LWNEMOCOUSTK, LWNEMOCOUSTRN, LWNEMOTAUOC,       &
     &   LWNEMOCOURECV,                                                 &
     &   LWNEMOCOUCIC, LWNEMOCOUCIT, LWNEMOCOUCUR,                      &
     &   LWNEMOCOUDEBUG,                                                &
     &   LLCAPCHNK,                                                     &
     &   LWAM_USE_IO_SERV


      CHARACTER(LEN=14) :: CLOUT
      NAMELIST /NAOT/ CLOUT

      CHARACTER(LEN=14) :: CLSOUT
      NAMELIST /NAOS/ CLSOUT

      CHARACTER(LEN=14) :: CLAOUT
      NAMELIST /NAAT/ CLAOUT

      INTEGER :: IDWI, IDWO
      CHARACTER(LEN=14) :: CLWOUT
      NAMELIST /NAWI/ IDWI, IDWO, CLWOUT

!     NAMELIST NALINE : 
!     ===============
!     CBPLTDT: USER INPUT START DATE OF RUN.
!     CEPLTDT: USER INPUT END DATE OF RUN.
!     CDATEF: BEGIN DATE OF FORECAST.
!     IDELPRO: PROPAGATION TIME STEP.
!              IF COUPLED TO IFS, THE INPUT VALUE OF IDELPRO WILL BE
!              THE UPPER BOUND ON THE TIME STEP, WHICH WOULD OTHERWISE
!              BE SET TO THE COUPLING TIME STEP.
!     IDELT: SOURCE TERM INTEGRATION TIME STEP.
!              IF COUPLED TO IFS, THE INPUT VALUE OF IDELT WILL BE
!              THE UPPER BOUND ON THE TIME STEP, WHICH WOULD OTHERWISE
!              BE SET TO THE COUPLING TIME STEP.
!     IDELWO: WIND OUTPUT TIME STEP.
!             (if specific values are not specified by namelist NAWI).
!     IDELWI: WIND INPUT TIME STEP.
!             (if specific values are not specified by namelist NAWI).
!     CLMTSU: STEP UNIT (S : seconds or H : hours).
!     IDELALT: ALTIMETER DATA TIME WINDOW (in seconds).
!     IDELINT: OUTPUT TIME STEP FOR INTEGRATED PARAMETER OF TOTAL SEA
!              (if output times are not specified by namelist NAOT).
!     IDELRES: OUTPUT TIME STEP for RESTART SPECTRA (if output times
!              are not specified by namelist NAOS).
!     IDELCUR: CURRENTS INPUT TIME STEP IN SECONDS. 
!     CDATECURA: FIRST DATE FOR THE CURRENT INPUT.
!     LLCFLCUROFF: IF TRUE PREVENT CFL CRITERIA FAILURE DUE TO SURFACE CURRENTS
!                  IT RESETS THE CURRENT REFRACTION TERMS TO 0 FOR THOSE POINTS
!                  WHERE CFL WAS NOT SATISFIED (it should only be used in operational applications)
!     CLOTSU: STEP UNIT (S : seconds or H : hours) FOR IDELINT,
!             IDELRES, IDELBC.
!     CDATER: ONE SPECIFIC OUTPUT TIME FOR RESTART FILES if output times
!             are not specified by namelist NAOS).
!     CDATES: LAST OUTPUT TIME ALLOWED BESIDE CDATER FOR RESTART SPECTRA
!             (if output times are not specified by namelist NAOS).
!     FFLAG: OUTPUT FLAG FOR OUTPUT TO USER OUTPUT OF EACH OUTPUT TYPE.
!     GFLAG: OUTPUT FLAG FOR OUTPUT TO GRIB OF EACH OUTPUT TYPE.
!     NFLAG: OUTPUT FLAG FOR USER OUTPUT DISPLAY OF FIELD NORM FOR
!            EACH OUTPUT TYPE.
!     LLOUTERS : IF TRUE CALL OUTERS: OUTPUT OF SATELLITE COLOCATION SPECTRA
!     TYPE OF INTEGRATED PARAMETERS IN FFLAG GFLAG (see OUTINT) :
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
!     38 : MAXIMUM WAVE HEIGHT (m)
!     39 : MAXIMUM WAVE PERIOD (S)
!     40 : U-COMPONENT OF STOKES DRIFT (m/s) 
!     41 : V-COMPONENT OF STOKES DRIFT (m/s) 
!     42 : U-COMPONENT OF SURFACE CURRENT (m/s) 
!     43 : V-COMPONENT OF SURFACE CURRENT (m/s) 
!     44 : NORMALIZED ENERGY FLUX INTO OCEAN 
!     45 : NORMALIZED ENERGY FLUX INTO WAVES 
!     46 : NORMALIZED STRESS INTO OCEAN 
!     47 : TRAIN 1 WAVE HEIGHT (M)
!     48 : TRAIN 1 MEAN WAVE DIRECTION (DEG.)
!     49 : TRAIN 1 WAVE MEAN PERIOD BASED ON THE 1/f INTEGRATION OF F (s.)
!     50 : TRAIN 2 WAVE HEIGHT (M)
!     51 : TRAIN 2 MEAN WAVE DIRECTION (DEG.)
!     52 : TRAIN 2 WAVE MEAN PERIOD BASED ON THE 1/f INTEGRATION OF F (s.)
!     53 : TRAIN 3 WAVE HEIGHT (M)
!     54 : TRAIN 3 MEAN WAVE DIRECTION (DEG.)
!     55 : TRAIN 3 WAVE MEAN PERIOD BASED ON THE 1/f INTEGRATION OF F (s.)
!     56 : MEAN SQUARE WAVE STRAIN IN THE SEA ICE.  
!     57 : SIGNIFICANT WAVE HEIGHT WITH PERIOD > 10s (m)
!     58 : SURFACE AIR DENSITY (kg m**-3)
!     59 : CONVECTIVE VELOCITY SCALE (m/s)
!     60 : SEA ICE COVER (-)
!     61 : SEA ICE THICKNESS (m)
!     62 : SPECTRAL SKWENESS (-)
!     63 : NEMO SST (K)
!     64 : NEMO SEA ICE COVER (-)
!     65 : NEMO SEA ICE THICKNESS (m)
!     66 : NEMO ZONAL CURRENT (m/s)
!     67 : NEMO MERIDIONAL CURRENT (m/s)
!     68 : WAVE ENERGY FLUX MAGNITUDE  (W/m)
!     69 : WAVE ENERGY FLUX MEAN DIRECTION (Degree true)
!     70 : SIGNIFICANT WAVE HEIGHT OF ALL WAVES WITH PERIOD BETWEEN 10 AND 12 SECONDS (m)
!     71 : SIGNIFICANT WAVE HEIGHT OF ALL WAVES WITH PERIOD BETWEEN 12 AND 14 SECONDS (m)
!     72 : SIGNIFICANT WAVE HEIGHT OF ALL WAVES WITH PERIOD BETWEEN 14 AND 17 SECONDS (m)
!     73 : SIGNIFICANT WAVE HEIGHT OF ALL WAVES WITH PERIOD BETWEEN 17 AND 21 SECONDS (m)
!     74 : SIGNIFICANT WAVE HEIGHT OF ALL WAVES WITH PERIOD BETWEEN 21 AND 25 SECONDS (m)
!     75 : SIGNIFICANT WAVE HEIGHT OF ALL WAVES WITH PERIOD BETWEEN 25 AND 30 SECONDS (m)

!     JPPFLAG-4 : 1st EXTRA EXPERIMENTAL FIELD
!     JPPFLAG-3 : 2nd EXTRA EXPERIMENTAL FIELD
!     JPPFLAG-2 : 3rd EXTRA EXPERIMENTAL FIELD
!     JPPFLAG-1 : 4th EXTRA EXPERIMENTAL FIELD
!     JPPFLAG   : 5th EXTRA EXPERIMENTAL FIELD

!     LFDB: TRUE IF INTEGRATED PARAMETERS ARE WRITTEN OUT TO FDB.
!     LGRIBIN : IF TRUE THE WAVE SPECTRA IS INPUT IN GRIB FORMAT, ELSE
!               THE BINARY RESTART FILES ARE USED.
!     LGRIBOUT : IF TRUE THE WAVE SPECTRA IS OUTPUT IN GRIB FORMAT, ELSE
!                THE BINARY RESTART FILES ARE PRODUCED.
!     LFDBIOOUT : IF TRUE THE GRIB SPECTRA OUTPUT IS SENT TO THE FDB.
!     NWRTOUTWAM: STRIDE WITH WHICH FDB OUTPUT PE's ARE SELECTED.
!     LRSTPARALW: IF TRUE BINARY RESTART FILES WILL BE WRITTEN in PARALLEL
!                 IF THAT OPTION IS USED.
!     LRSTPARALR: IF TRUE BINARY RESTART FILES WILL BE READ in PARALLEL
!                 IF THAT OPTION IS USED.
!     LRSTINFDAT:  F TRUE WILL WRITE AN ADDITIONAL wamfile WITH DATE/TIME INFO
!     LSECONDORDER : IF TRUE THEN SECOND ORDER CORRECTIOn IS APPLIED TO
!                    OUTPUT INTEGRATED PARAMETERS.
!     ICASE: 1 FOR SPHERICAL COORDINATES ELSE CARTESIAN COORDINATES.
!     ISHALLO: 1 FOR DEEP WATER MODE ELSE SHALLOW WATER MODEL.
!     IREFRA: 0 MODEL RUNS WITHOUT REFRACTION.
!             1 MODEL RUNS WITH DEPTH REFRACTION ONLY.
!             2 MODEL RUNS WITH CURRENT REFRACTION ONLY.
!             3 MODEL RUNS WITH DEPTH AND CURRENT REFRACTION.
!     ITEST: TEST OUTPUT LEVEL.
!     ITESTB: MAX BLOCK NUMBER FOR OUTPUT IN BLOCK LOOPS.
!     IREST: 1 FOR THE PRODUCTION OF RESTART FILE(S).
!     IASSI: 1 ASSIMILATION IS DONE IF ANALYSIS RUN.
!     ISNONLIN : 0 FOR OLD SNONLIN, 1 FOR NEW SNONLIN.
!     IDAMPING : 0 NO WAVE DAMPING, 1 WAVE DAMPING ON.
!                ONLY MEANINGFUl FOR IPHYS=0
!     IBOUNC: 1 FOR RUN WHICH INCLUDES BOUNDARY POINTS (COARSE GRID).
!     IBOUNF: 1 FOR RUN WHICH INCLUDES BOUNDARY POINTS (FINE GRID).
!     IDELBC: TIMESTEP FOR THE DISPOSAL OF OUTPUT BC FILE(S) INTO PERMANENT
!             FILE(S). IF SET TO ZERO THEN IT WILL BE RESET TO THE LENGTH
!             OF THE RUN (I.E. THE FILE(S) WILL DISPOSED AT THE END OF THE
!             RUN.
!     CBCPREF: PREFIXES OF BC FILES
!     USERID, RUNID, PATH : OUT OF DATE, SHOULD ONLY BE USED WHENEVER
!                           USER WANTS TO WRITE OUTPUT DIRECTLY TO ECFS
!                           WHICH IS HIGHLY NOT RECOMMENDATED ON THE VPP
!     YCLASS: DATA CLASS USED TO CODE DATA IN GRIB.
!     YEXPVER: EXPERIMENT VERSION USED TO CODE DATA IN GRIB.
!     CPATH: PATH FOR OUTPUT TO DISK (USED WHENEVER OUTPUT IS NOT
!            TO THE FIELD DATA BASE (SEE LFDBIOOUT)).
!     IMDLGRBID_G: GLOBAL MODEL IDENTIFICATION NUMBER FOR GRIB CODING.
!                  DEFAULT VALUE SET IN YOWGRIBHD.
!     IMDLGRBID_M: LAW MODEL IDENTIFICATION NUMBER FOR GRIB CODING.
!                  DEFAULT VALUE SET IN YOWGRIBHD.
!     NTOTENS: TOTAL ENSEMBLE FORECAST MEMBERS (DEFAULT=0).
!     NENSFNB: ENSEMBLE FORECAST NUMBER (DEFAULT=0).
!     NSYSNB : SYSTEM NUMBER TO BE USED FOR GRIBBING OF SEASONAL DATA.
!              or MONTHLY FORECAST RUNS
!     NMETNB : METHOD NUMBER TO BE USED FOR GRIBBING OF SEASONAL DATA.
!              or MONTHLY FORECAST RUNS
!     LMESSPASS: TRUE FOR MESSAGE PASSING ARCHITECHTURE.
!     LWCOU: FALSE FOR UNCOUPLED RUN.
!     NEMO COUPLING FLAGS:
!     LWNEMOCOU: FALSE FOR NO COUPLING TO NEMO RUN.
!     NEMOFRCO: NEMO COUPLING FREQ IN WAM TIME STEPS
!     LWNEMOCOUSEND: IF FALSE THEN NO DATA WILL BE SENT TO NEMO
!     LWNEMOCOUSTK:  IF TRUE SEND SURFACE STOKES DRIFT TO NEMO
!     LWNEMOCOUSTRN: IF TRUE SEND SEA ICE WAVE STRAIN TO NEMO
!     LWNEMOTAUOC IF TRUE THEN THE ATMOSPHERIC STRESS THAT IS PASSED TO
!                 NEMO VIA THE WAVE MODEL WILL ALSO BE MODULATED BY IT.
!
!     LWNEMOCOURECV: IF FALSE THEN NO DATA WILL BE RECEIVED FROM NEMO
!     LWNEMOCOUCIC: IF TRUE THEN SEA ICE CONCENTRATION WILL BE OBTAINED FROM FROM NEMO 
!     LWNEMOCOUCIT: IF TRUE THEN SEA ICE THICKNESS WILL BE OBTAINED FROM FROM NEMO 
!     LWNEMOCOUCUR: IF TRUE THEN SURFACE CURRENTS WILL BE OBTAINED FROM FROM NEMO 
!
!     LWNEMOCOUDEBUG: FALSE IF NO DEBUGGING OUTPUT IN WAM<->NEMO COUPLING
!
!     LLCAPCHNK : CAP CHARNOCK FOR HIGH WINDS.
!     LWAM_USE_IO_SERV: TRUE IF SPECTRAL AND INTEGRATED PARAMETER OUTPUT SHOULD BE
!              DONE USING IFS IO SERVER
!
!     LNOCDIN: IF TRUE THEN GRIB INPUT OF A DRAG COEFFICIENT FIELD IS
!              NOT REQUIRED.  
!     LODBRALT: IF TRUE THEN THE ALTIMETER DATA WILL BE READ AND PASSED 
!               THROUGH OBSERVATION DATABASE (ODB)
!     LALTCOR: IF TRUE THEN THE ALTIMETER DATA WILL BE CORRECTED
!              SEE GRFIELD. 
!     L4VTYPE: IF TRUE THEN MARS TYPE 4V IS USED INSTEAD OF FG.
!     LFRSTFLD: IF TRUE THEN INITIAL INTEGRATED PARAMETER FIELDS WILL
!               OUTPUT PROVIDED THEY WERE NOT INPUT.
!     XKAPPA2: KAPPA2 PARAMETER USED IN THE FORMULA FOR THE
!              DETERNMINATION OF ALTIMETER WAVE HEIGHT.
!     LWAMANOUT: CONTROLS WHETHER FIELDS WILL BE WRITTEN OUT AT ANALYSIS TIME OR NOT.
!                (USEFUL WHEN WAVE DATA ASSIMILATION IS DONE IN MORE THAN ONE TRAJECTORY.)
!     LALTAS: CONTROLS WHETHER ALTIMETER DATA ARE ASSIMILATED
!     LSARAS: CONTROLS WHETHER SAR DATA ASSIMILATED
!     LSARINV: CONTROLS WHETHER SAR INVERSION IS DONE
!     SWAMPWIND : CONSTANT WIND SPEED USED BY SWAMP CASE.
!     SWAMPWIND2 : SECOND WIND SPEED VALUE TO BE USED (if different than 0)
!                  IN SWAMP CASE. IT WILL BE APPLIED DTNEWWIND HOURS
!                  AFTER STATING TIME.
!     SWAMPCIFR : CONSTANT SEA ICE FRACTION FOR THE 1/2 NORTHERN PART OF SWAMP DOMAIN
!     SWAMPCITH : CONSTANT SEA ICE THICKNESS FOR THE 1/2 NORTHERN PART OF SWAMP DOMAIN
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
!     LALTPAS: FLAG USED TO CONTROL FEED THROUGH ALT DATA PASSIVELY. 
!                (SEE GRFIELD). 
!     LGUST:   FLAG USED TO ACTIVATE COMPUTATIONS RELATED TO GUSTINESS 
!     LADEN:   FLAG USED TO ACTIVATE COMPUTATIONS RELATED TO AIR DENSITY
!     LRELWIND: IF TRUE THEN RELATIVE WINDS ARE USED WITH RESPECT TO
!               SURFACE CURRENTS.
!     LLWSWAVE: FLAG USE TO ACTIVATE USE OF WAVE PARAMETER WIND SPEED AS
!               INPUT TO CONSTRUCT THE WIND FORCING - UNCOUPLED RUNS ONLY
!     LLWDWAVE: FLAG USE TO ACTIVATE USE OF WAVE PARAMETER WIND DIRECTION
!               AS INPUT TO CONSTRUCT THE WIND FORCING - UNCOUPLED RUNS ONLY.
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
!     NPROMA_WAM: NUMBER OF THE GRID POINTS PER LOOP WHEN CUT INTO CHUNKS
!     LL1D: IF TRUE THEN THE DOMAIN IS ONLY DIVIDED INTO LATITUDINAL BANDS
!     LWCOUNORMS : FLAG THAT CONTROLS WHETHER OR NOT GLOBAL NORMS
!                  ARE PRODUCED FOR THE DIFFERENT FIELDS EXCHANGED 
!                  BETWEEN WAM AND THE COUPLED ATMOSPHERIC MODEL.
!                  IT IS FALSE BY DEFAULT.
!     LLNORMIFS2WAM : IF TRUE NORMS FOR FIELDS PASSED FROM IFS TO WAM WILL BE PRODUCED
!     LLNORMWAM2IFS : IF TRUE NORMS FOR FIELDS PASSED FROM WAM TO IFS WILL BE PRODUCED
!     LLNORMWAMOUT : IF TRUE NORMS OF SELECTED OUPTUT FIELDS WILL BE PRODUCED
!     LLNORMWAMOUT_GLOBAL : IF TRUE NORMS FOR SELECTED OUPTUT FIELDS WILL BE GLOBAL (see above)
!     LGRHDIFS : FLAGS CONTROLLING WHETHER OR NOT GRIB HEADER INFORMATION
!                IS COPIED FROM THE ATMOSPHERIC MODEL (ONLY USEFULL IF
!                COUPLED TO THE IFS).
!     LNEWLVTP : FLAG CONTROLLING WHETHER OR NOT THE NEW GRIB LEVTYPE
!                DEFINITIONS ARE USED.
!     LICERUN : FLAG CONTROLLING WHETHER OR NOT SEA ICE FRACTION (OR SST)
!               FIEDS ARE PROVIDED WITH THE WIND FIELDS TO GENERATE THE
!               SEA ICE MASK (TRUE BY DEFAULT). 
!     LCIWABR : FLAG CONTROLLING THE USE OF SEA ICE BOTTOM FRICTION ATTENUATION  
!     LICETH  : FLAG CONTROLLING WHETHER OR NOT SEA ICE THICKNESS FILEDS
!               ARE PROVIDED WITH THE WIND FIELDS (FALSE BY DEFAULT). 
!     LLSOURCE : FLAG CONTROLLING WHETHER OR NOT THE SOURCE TERM CONTRIBUTION 
!                IS COMPUTED.
!     LNSESTART : FLAG CONTROLLING WHETHER OR NOT THE INITIAL SPECTRA ARE
!                 RESET TO NOISE LEVEL.
!     LSMSSIG_WAM : .T. = send signals to SMS or ECFLOW (ECMWF supervisor)
!     CMETER :  SMS or ECFLOW meter command (ECMWF supervisor)
!     CEVENT :  SMS or ECFLOW event command (ECMWF supervisor)


!     NAMELIST NAOT : 
!     ===============
!     CLOUT : LIST OF OUTPUT TIME FOR INTEGRATED PARAMETERS.
      
!     NAMELIST NAOS : 
!     ===============
!     CLSOUT : LIST OF OUTPUT TIME FOR SPECTRA. 

!     NAMELIST NAAT : 
!     ===============
!     CLAOUT : LIST OF TIME FOR ALTIMETER DATA ASSIMILATION. 

!     NAMELIST NAWI (only meaningful for uncoupled runs): 
!     ===============
!     IDWI : LIST OF IDELWI
!     IDWO : LIST OF IDELWO
!     CLWOUT : LIST OF TIME UNTIL WHICH IDWI AND IDWO ARE VALID. 

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('MPUSERIN',0,ZHOOK_HANDLE)

!*    0. SET DEFAULT VALUES FOR THE NAMELIST ELEMENTS.
!        ---------------------------------------------

      CLMTSU    = 'S'
      CLOTSU    = 'H'
      CLHEADER  =  ZERO
      CBPLTDT   =  ZERO
      CEPLTDT   =  ZERO
      CDATEF    =  ZERO
      IDELPRO   =  0
      IDELT     =  0
      IDELWO    =  0
      IDELWI    =  0
      IDELALT   =  21600 
      IDELINT   =  0

!!!! obsolete should be removed as soon as I can update scripts
      IDELINS   =  0
      IDELSPT   =  0
      IDELSPS   =  0
!!!!

      IDELRES   =  0
      IDELCUR   =  0
      CDATECURA = ZERO
      LLCFLCUROFF = .TRUE.
      CDATER    = ZERO
      CDATES    = ZERO
      FFLAG(:)  = .FALSE. 
!     note: if IREST>0 then GFLAG(IR) for IR=IRWDIR, IRCD, IRU10 will always be reset
!     true if grib restart files are requested (LGRIBOUT=.true.).
!     see below.
      GFLAG(:)  = .FALSE.
      NFLAG(:)  = .FALSE. 
!     patial initialisation of the output table
      CALL MPCRTBL
      NFLAG(1)  = .TRUE. 
      NFLAG(6)  = .TRUE. 
      NFLAG(IRWDIR)= .TRUE. 
      NFLAG(IRCD)  = .TRUE. 
      NFLAG(IRU10) = .TRUE. 

      LLOUTERS = .FALSE.

      LSECONDORDER = .TRUE.

      LFDB      = .TRUE. 
      LGRIBIN   = .TRUE. 
      LGRIBOUT  = .TRUE. 
      LFDBIOOUT = .TRUE. 
      NWRTOUTWAM= 1
      LRSTPARALW= .TRUE.
      LRSTPARALR= .TRUE.
      LRSTINFDAT= .FALSE.
      LODBRALT  = .FALSE.
      ICASE     = 1 
      ISHALLO   = 0 
      ISNONLIN  = 1 
      IDAMPING  = 1 
      IPROPAGS  = 0 
      IREFRA    = 0 
      ITEST     = 0 
      ITESTB    = 0 
      IREST     = 0 
      IASSI     = 0 
      IBOUNC    = 0 
      IBOUNF    = 0 
      IDELBC    = 0
      DO II=1,GBOUNC_MAX
        WRITE(C2,'(I2.2)') II
        CBCPREF(II)='B'//C2
      ENDDO
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
      LMESSPASS = .TRUE.

      NOUTT     = 0

      LNOCDIN   = .FALSE. 

      DO ISAT=1,NUMALT
        IBUFRSAT(ISAT) = 0
        CSATNAME(ISAT) = 'SATELLITE NAME MISSING'
        LALTPAS(ISAT)   = .FALSE. 
        LALTCOR(ISAT)   = .FALSE. 
        XKAPPA2(ISAT)   = 0.0700_JWRB
        LALTLRGR(ISAT)  = .FALSE.
        HSCOEFCOR(ISAT) = 1.0_JWRB
        HSCONSCOR(ISAT) = 0.0_JWRB

        ALTBGTHRSH(ISAT) = 1.0_JWRB
!       if no value is provided in the namelist ALTSDTHRSH will
!       be set in grfield.
        ALTSDTHRSH(ISAT) = -1.0_JWRB

!       HSALTCUT is used in combination with the error estimate of
!       the altimeter data to determine the minimum Hs allowed for
!       altimeter data.
        HSALTCUT(ISAT) = 999999.
        LALTGRDOUT(ISAT)  = .FALSE.
      ENDDO

      L4VTYPE   = .FALSE. 

      LFRSTFLD  = .FALSE.

      LALTAS    = .FALSE. 
      LSARAS    = .FALSE. 
      LSARINV   = .FALSE. 

      LWAMANOUT = .TRUE.

      SWAMPWIND = 18.45_JWRB
      SWAMPWIND2 = 0.0_JWRB
      DTNEWWIND = 0.0_JWRB
      SWAMPCIFR = 0.0_JWRB
      SWAMPCITH = 0.0_JWRB
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

      LLWSWAVE = .FALSE.
      LLWDWAVE = .FALSE.

      LSUBGRID = .TRUE.

      NPROMA_WAM = 1
      NPROMA_WAM = HUGE(NPROMA_WAM)/2

      LL1D = .TRUE.

      LWCOU = .FALSE.

      LWVFLX_SNL = .TRUE.

      LWNEMOCOU = .FALSE.

      NEMOFRCO = 0

      LWNEMOCOUSEND = .TRUE.

      LWNEMOCOUSTK=.TRUE.

      LWNEMOCOUSTRN=.FALSE.

      LWNEMOTAUOC = .TRUE.

      LWNEMOCOURECV = .FALSE.

      LWNEMOCOUCIC=.FALSE.

      LWNEMOCOUCIT=.FALSE.

      LWNEMOCOUCUR=.FALSE.

      LWNEMOCOUDEBUG = .FALSE.

      LLCAPCHNK = .FALSE.

      LWAM_USE_IO_SERV = .FALSE.

      LWCOUNORMS = .FALSE.
      LLNORMIFS2WAM = .FALSE.
      LLNORMWAM2IFS = .FALSE.
      LLNORMWAMOUT = .FALSE.
      LLNORMWAMOUT_GLOBAL = .FALSE.

      LGRHDIFS = .FALSE.

      LNEWLVTP = .FALSE.

      LBIWBK = .TRUE.

      LICERUN = .TRUE.

      LCIWABR = .TRUE.

      LMASKICE = .FALSE.

      LICETH = .FALSE.

      LLSOURCE = .TRUE.

      LNSESTART = .FALSE.

      LRELWIND = .TRUE.

      NDELW_LST = 0

      LSMSSIG_WAM =.FALSE.

      CMETER='smsmeter'

      CEVENT='smsevent'

      LLUNSTR=.FALSE.
      LVECTOR=.FALSE.
      IVECTOR=1
      LPREPROC=.FALSE.
      USE_DIRECT_WIND_FILE=.FALSE.
      JGS_DIFF_SOLVERTHR = 1.E-5_JWRU
      WAE_SOLVERTHR = 1.E-10_JWRU
      LIMPLICIT = .FALSE.
      SOURCE_IMPL = .FALSE.
      LNONL = .FALSE.
      BLOCK_GAUSS_SEIDEL = .TRUE.
      LLIMT = .FALSE.
      L_SOLVER_NORM = .FALSE.
      LCHKCONV = .TRUE.
      LBCWA = .FALSE.

      LLOPENED = .FALSE.


! ----------------------------------------------------------------------

!*    1. READ NAMELIST NALINE.
!        ---------------------

      NAMELIST_FILENAME='wam_namelist' 
      WRITE(IU06,*) ' '
      IU05 =  IWAM_GET_UNIT (IU06, NAMELIST_FILENAME, 's', 'f', 0, 'READ')
  
      CALL WPOSNAM (IU05, 'NALINE', LLEOF)
      IF (.NOT. LLEOF) THEN
        READ (IU05, NALINE)
      ELSE
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'+ SUBROUTINE MPUSERIN :                    +'
        WRITE(IU06,*)'+ READ NAMELIST FAILED                     +' 
        WRITE(IU06,*)'+ NAMELIST FILENAME: ',NAMELIST_FILENAME
        WRITE(IU06,*)'+ PROGRAM WILL ABORT                       +'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        CALL ABORT1
      ENDIF

      IF (LWCOU) LSMSSIG_WAM=.FALSE.


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
      IF (.NOT.ALLOCATED(COUTT)) THEN
        ALLOCATE(COUTT(NOUTT))
      ENDIF
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
      IF (.NOT.ALLOCATED(COUTS)) THEN
        ALLOCATE(COUTS(NOUTS))
      ENDIF
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
      IF (.NOT.ALLOCATED(CASS)) THEN
        ALLOCATE(CASS(NASS))
      ENDIF
      SPA: DO IC=1,NASS
        CALL WPOSNAM (IU05, 'NAAT', LLEOF)
        IF (LLEOF) EXIT SPA
        READ (IU05, NAAT, END=1912)
        CASS(IC) = CLAOUT
      ENDDO SPA
1912  CONTINUE

      IF (.NOT.LWCOU) THEN
!           **** FORCING TIME STEPPING AT SPECIFIED TIMES ****
        REWIND(IU05)
        SPWI: DO
          CALL WPOSNAM (IU05, 'NAWI', LLEOF)
          IF (LLEOF) EXIT SPWI
          READ (IU05, NAWI, END=1903)
          NDELW_LST = NDELW_LST + 1 
        ENDDO SPWI
1903    CONTINUE
        REWIND(IU05)
        IF (.NOT.ALLOCATED(IDELWI_LST)) THEN
          ALLOCATE(IDELWI_LST(NDELW_LST))
        ENDIF
        IF (.NOT.ALLOCATED(IDELWO_LST)) THEN
          ALLOCATE(IDELWO_LST(NDELW_LST))
        ENDIF
        IF (.NOT.ALLOCATED(CDTW_LST)) THEN
          ALLOCATE(CDTW_LST(NDELW_LST))
        ENDIF
        SPW: DO IC=1,NDELW_LST
          CALL WPOSNAM (IU05, 'NAWI', LLEOF)
          IF (LLEOF) EXIT SPW
          READ (IU05, NAWI, END=1913)
          IDELWI_LST(IC) = IDWI 
          IDELWO_LST(IC) = IDWO 
          CDTW_LST(IC) = CLWOUT
        ENDDO SPW
1913    CONTINUE
      ENDIF

      CLOSE(IU05)

!           **** MODEL TIME STEPS ****
      IF (CLMTSU(1) .EQ. 'H') IDELPRO = IDELPRO*3600
      IF (CLMTSU(2) .EQ. 'H') IDELT   = IDELT*3600
      IF (CLMTSU(3) .EQ. 'H') IDELWO  = IDELWO*3600
      IF (CLMTSU(4) .EQ. 'H') IDELWI  = IDELWI*3600
!           **** OUTPUT TIME IN FIXED INTERVALS ****
      IF (CLOTSU(1) .EQ. 'H') IDELINT = IDELINT*3600
      IF (CLOTSU(7) .EQ. 'H') IDELRES = IDELRES*3600
      IF (CLOTSU(8) .EQ. 'H') IDELBC  = IDELBC*3600


!     RESET CERTAIN FLAGS:

!     WE SHOULD RECEIVE DATA FROM NEMO
      IF (LWNEMOCOUCIC.OR.LWNEMOCOUCIT.OR.LWNEMOCOUCUR)                 &
     &   LWNEMOCOURECV = .TRUE.


! Here we set LL1D = .TRUE. for the case of LLUNSTR in order to omit the mapping for the parallel strucutured grid
      IF (LLUNSTR) LL1D = .TRUE.
      IF (LWCOU) LSMSSIG_WAM =.FALSE. ! by definition

!     Most of the namelist selection will be written to the logfiles in userin.

!     Some are printed below

      IF (IRANK.EQ.1) THEN
        WRITE(6,*) '==============================================='
        WRITE(6,*) '*** MPUSERIN has read the following settings'
        WRITE(6,*) '*** LMESSPASS = ',LMESSPASS
        WRITE(6,*) '*** LFDB = ',LFDB
        WRITE(6,*) '*** LRSTPARALW = ',LRSTPARALW
        WRITE(6,*) '*** LRSTPARALR = ',LRSTPARALR
        WRITE(6,*) '*** LRSTINFDAT = ',LRSTINFDAT
        WRITE(6,*) '*** LWCOU= ',LWCOU
        WRITE(6,*) '*** LWNEMOCOU= ',LWNEMOCOU
        WRITE(6,*) '*** LWNEMOCOUSEND  = ',LWNEMOCOUSEND
        WRITE(6,*) '*** LWNEMOCOURECV  = ',LWNEMOCOURECV
        WRITE(6,*) '*** LWNEMOCOUDEBUG = ',LWNEMOCOUDEBUG
        WRITE(6,*) '*** LWNEMOCOUCIC   = ',LWNEMOCOUCIC
        WRITE(6,*) '*** LWNEMOCOUCIT   = ',LWNEMOCOUCIT
        WRITE(6,*) '*** LWNEMOCOUCUR   = ',LWNEMOCOUCUR
        WRITE(6,*) '*** LWNEMOCOUSTK   = ',LWNEMOCOUSTK
        WRITE(6,*) '*** LWNEMOCOUSTRN  = ',LWNEMOCOUSTRN
        WRITE(6,*) '*** LWNEMOTAUOC    = ',LWNEMOTAUOC
        WRITE(6,*) '*** LSUBGRID= ',LSUBGRID
        WRITE(6,*) '*** IPROPAGS= ',IPROPAGS
        WRITE(6,*) '*** IREFRA= ',IREFRA
        WRITE(6,*) '*** LLUNSTR= ',LLUNSTR
        WRITE(6,*) '*** LVECTOR= ',LVECTOR
        WRITE(6,*) '*** IVECTOR= ',IVECTOR
        WRITE(6,*) '*** LPREPROC= ',LPREPROC
        WRITE(6,*) '*** LWCOUNORMS= ',LWCOUNORMS
        WRITE(6,*) '*** LLNORMIFS2WAM= ',LLNORMIFS2WAM
        WRITE(6,*) '*** LLNORMWAM2IFS= ',LLNORMWAM2IFS
        WRITE(6,*) '*** LLNORMWAMOUT= ',LLNORMWAMOUT
        WRITE(6,*) '*** LLNORMWAMOUT_GLOBAL= ',LLNORMWAMOUT_GLOBAL
        WRITE(6,*) '*** LSMSSIG_WAM= ',LSMSSIG_WAM
        WRITE(6,*) '*** LWAM_USE_IO_SERV = ',LWAM_USE_IO_SERV
        WRITE(6,*) '==============================================='
      ENDIF

      IF (LHOOK) CALL DR_HOOK('MPUSERIN',1,ZHOOK_HANDLE)

      END SUBROUTINE MPUSERIN
