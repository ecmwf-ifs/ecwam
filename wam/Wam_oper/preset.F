      PROGRAM preset 

! ----------------------------------------------------------------------

!**** *PRESET* - GENERATES ALL BINARY FILES REQUIRED FOR A WAMODEL START
!             OR GENERATES SPECTRAL GRIB FILE REQUIRED FOR A WAMODEL
!             START. IN THAT CASE NO DRAG COEFFICIENT FIELD IS REQUIRED
!             TO RUN WAMODEL (SEE LNOCDIN IN WAMODEL).
!             OR CONVERTS AND INTERPOLATES GRIB SPECTRA PARAMETER 250 
!             INTO PARAMETER 251. THE CORRESPONDING DRAG COEFFICIENT 
!             FIELD SHOULD BE PROVIDED WITH AN APPROPRIATE MARS REQUEST.

!     SUSANNE HASSELMANN  MPI     JULY 1986.
!     ANNEGRET SPEIDEL    MPI      MAY 1988 PARAMETER STATEMENTS.
!     ANNEGRET SPEIDEL    MPI NOVEMBER 1988 CRAY-2 VERSION.
!     CYCLE_3 MODICIFATIONS:
!     ----------------------
!     RENATE PORTZ       MPI      JUNE 1990 COMPUTATION OF INITIAL
!                                           JONSWAP SPECTRA FROM
!                                           INITIAL WIND FIELD.
!     CYCLE_4 MODICIFATIONS:
!     ----------------------
!     H. GUNTHER  GKSS/ECMWF  DECEMBER 1990
!     J. BIDLOT    ECMWF   FEBRUARY 1996  MESSAGE PASSING
!     J. BIDLOT    ECMWF   MARCH 1997  MODIFY ROUTINES FOR OUTPUT OF
!                          RESTART FILE 
!     B. HANSEN    ECMWF   APRIL    1997  RESTART FACILITY.
!     B. HANSEN    ECMWF   JANUARY  1998  NAMELIST INPUT.
!     B. HANSEN    ECMWF   FEBRUARY 1998  WRITE 2DSP TO FDB.
!     J. BIDLOT    ECMWF   MARCH 1998 PRODUCE GRIB SPECTRA (215) AND
!                                     INTERPOLATE TO OUTPUT GRID.
!                                     NO LAW FILE IS PRODUCED THEN.
!     J. BIDLOT    ECMWF   OCTOBER 1998   MODULES.
!*    PURPOSE.
!     --------
!       TO INITIALISE ALL FILES REQUESTED BY THE WAMODEL.
!**   INTERFACE.
!     ----------
!       *IU01*   INTEGER    INPUT  UNIT UNBLOCKED WIND FILE.
!                           (SEE SUB READWIND).
!       *IU05*   INTEGER    USER INPUT UNIT.
!       *IU06*   INTEGER    PRINTER OUTPUT.
!       *IU07*   INTEGER    INPUT  UNIT PREPROC GRID OUTPUT.
!       *IU12*   INTEGER    OUTPUT UNIT BLOCKS OF SPECTRA.
!       *IU14*   INTEGER    OUTPUT UNIT SECOND LAT OF BLOCKS.
!       *IU15*   INTEGER    OUTPUT UNIT LAST WINDFIELDS.
!                 !OR!      OUTPUT UNIT LAST DRAG COEFFICIENT.
!     METHOD.
!     -------
!       IF USED TO PRODUCED A COLD START (LGRIBIN=FALSE) :

!       A JONSWAP SPECTRAL SHAPE IS ASSUMED.
!       JONSWAP PARAMETERS ARE DEFINED EITHER BY USER INPUT (IOPTI=0) OR
!       BY FETCH LAWS (IOPTI=1, or 2). THE 2-D SPECTRA ARE COMPUTED FOR
!       EACH POINT IN A BLOCK, THE WAMODEL BLOCKS ARE INITIALISED BY
!       THESE SPECTRA AND ALL BLOCKS AND OVERLAPPING LATITUDES ARE SAVED
!       IF FETCH LAWS ARE USED TO DEFINE PARAMETERS THE FIRST WIND
!       FIELD IS GENERATED OTHERWISE THE WIND FIELD IS INITIALISED
!       WITH ZEROS. THE MAIN MODEL WILL RECONSTRUCT THE WIND ANY HOW.
!       THE FILE HANDLING OF THE RESTART FILES IS COMPUTER DEPENDENT.
!       SUB GSFILE HAS TO BE MODIFIED, TO COPY THE UNIT ALIAS FILES
!       (UNITS IU12,IU14, AND IU15) TO PERMANENT FILES.


!       IF USED TO CONVERT 250 INTO 251 (LGRIBIN=TRUE) :

!       CONVERT AND INTERPOLATE SPECTRAL COMPONENTS SEE SPLITGRS.
!       GLOBAL FIELDS PRIOR 19961205 AND MEDITERANEAN FIELDS WERE CODED
!       IN GRIB AS 0NE SINGLE FIELD. GLOBAL FIELDS AFTER 19961205 WERE
!       SPLIT INTO TWO FIELDS. USE THE INPUT NAMELIST IOPTI TO SPECIFY
!       THE NUMBER OF FIELDS 250 IS PROVIDED WITH.
!       
!     EXTERNALS.
!     ----------
!       *ABORT1*    - TERMINATES PROCESSING.
!       *BUILDSTRESS- RECONSTRUCT THE LAW FILE FROM WINDS AND CD. 
!       *CREWFN*    - CREATES A WIND FILE NAME.
!       *GRSTNAME*  - BUILD NAME FOR INPUT/OUTPUT FILES 
!       *GSFILE*    - GETS OR SAVES FILES (COMPUTER DEPENDENT).
!       *INCDATE*   - UPDATES A DATE TINE GROUP.
!       *INTGRS*    - PERFORMS INTERPOLATIONS OF INPUT SPECTRA 
!       *LOCINT*    - INTERPOLATES TO MODEL GRID.
!       *MSTART*    - GENERATES THE RESTART FILES.
!       *NOTIM*     - CONTROLS WIND GENERATION (NO TIME INTERPOLATION).
!       *PEAK*      - COMPUTES PARAMETERS BY FETCH LAWS.
!       *PREWIND*   - PREPARES WINDS.
!       *SAVSPEC*   - SAVES THE RESTART SPECTRA.
!       *SAVSTRESS* - SAVES THE RESTART WIND AND STRESS FIELDS.
!       *SPECTRA*   - COMPUTES SPECTRA FROM PARAMETERS.
!       *SPLITBL*   - SAVES SECOND TO LAST LATITUDES.
!       *SPLITGRS*  - PERFORMS SPLITTING OF SPECTRA ON SUB DOMAINS 
!       *SPR*       - DIRECTIONAL DISTRIBUTION.
!       *READGRS*   - READS GRIBBED SPECTRUM DATA
!       *READPRE*   - READS PREPROC OUTPUT.
!       *READWIND*   - READS A WIND FIELD.
!       *TIMIN*     - CONTROLS WIND GENERATION (TIME INTERPOLATION).

!     REFERENCES
!     ----------
!       NONE.
! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3   ,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK0    ,KSEK1    ,
     &            KSEK2    ,KSEK3    ,KSEK4    ,PSEK2    ,PSEK3    ,
     &            PPMISS   ,PPEPS    ,PPREC    ,NTENCODE 
      USE YOWGRID  , ONLY : DELPHI   ,NLONRGG  ,IJL
      USE YOWJONS  , ONLY : FM       ,GAMMA    ,SA       ,SB       ,
     &            THETAQ
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IRGG     ,AMOWEP   ,
     &            AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO
      USE YOWMESPAS, ONLY : LMESSPASS,LPBIOIN  ,LPBIOOUT ,LFDBIOIN ,
     &            LFDBIOOUT,LGRIBIN  ,LGRIBOUT
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,
     &            KTAG     ,NPRECR   ,NPRECI
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,
     &            NBLO     ,NIBLO    ,CLDOMAIN
      USE YOWPCONS , ONLY : G        ,RAD      ,DEG      ,ZMISS
      USE YOWSTAT  , ONLY : MARSTYPE ,YCLASS   ,YEXPVER  ,CDATEA   ,
     &            CDATEE   ,CDATEF   ,CDTPRO   ,CDATER   ,CDATES   ,
     &            CFDB2DSP ,IDELPRO  ,IDELWI   ,IDELWO   ,IREFRA   ,
     &            NENSFNB  ,NTOTENS
      USE YOWSPEC  , ONLY : NSTART   ,NEND     ,U10OLD   ,THWOLD   ,
     &            USOLD    ,Z0OLD    ,TAUW     ,FL1
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWTEXT  , ONLY : ICPLEN   ,USERID   ,RUNID    ,PATH     ,
     &            CPATH
      USE YOWUNIT  , ONLY : IU12     ,IU14     ,IU15     ,IU04     ,
     &            IUSCR
      USE YOWWIND  , ONLY : CDA      ,CDAWIFL  ,CDATEWO  ,CDATEFL

! -------------------------------------------------------------------
!* 0.    COMMON OUTPUT CONTROL.
!  ----------------------------
#include <outctrl.h>
      DATA csubna /"PRESET "/

! ----------------------------------------------------------------------
      REAL X4(2)
      INTEGER I4(2)

      CHARACTER*  1 CLTUNIT
      CHARACTER* 12 ZERO, CDUM
      CHARACTER* 70 HEADER
      CHARACTER* 80 SFILENAME, ISFILENAME

      PARAMETER (ZERO=' ', CDUM='000000000000')

      LOGICAL LINTGRS

! DIMENSION DUMMY COUPLED VARIABLES

      INTEGER IUGB(1),IVGB(1)

! ----------------------------------------------------------------------

      NAMELIST /NALINE/ HEADER,
     &          IOPTI, ITEST, ITESTB,
     &          ALFA, FM, GAMMA, SA, SB, THETA, FETCH,
     &          USERID, RUNID, PATH, CPATH,
     &          CDATEA, IDELWI, CLTUNIT,
     &          LPBIOOUT, LFDBIOOUT, LGRIBIN, LGRIBOUT,
     &          MARSTYPE, YCLASS, YEXPVER, CFDB2DSP

!     IOPTI : WHEN LGRIBIN = FALSE, IT SELECTS COLD START SPECTRAL FORM
!             (0, 1, or 2 see MSTART)
!             OR INDICATES WHETHER INPUT SPECTRA IS CONTAINED IN 1 OR 2
!             FIELD(S) (WHEN LGRIBIN = TRUE). 
!     ITEST, ITESTB : TEST OUTPUT LEVEL SEE MODULE YOWTEST
!     ALFA, FM, GAMMA, SA, SB, THETA, FETCH : JONSWAP PARAMETERS
!                                             SEE MODULE YOWINTS
!     USERID, RUNID, PATH : OUT OF DATE, SHOULD ONLY BE USED WHENEVER 
!                           USER WANTS TO WRITE OUTPUT DIRECTLY TO ECFS
!                           WHICH IS HIGHLY NOT RECOMMENDATED ON THE VPP
!     CPATH : PATH FOR OUTPUT TO DISK (USED WHENEVER OUTPUT IS NOT 
!             TO THE FIELD DATA BASE (SEE LFDBIOOUT)).
!     CDATEA : INPUT DATE.
!     IDELWI : INPUT WIND TIME STEP IN UNIT DEFINED BY CLTUNIT
!     CLTUNIT : INPUT WIND TIME STEP UNIT (S : seconds or H : hours)
!     LPBIOOUT : IF TRUE THE PBIO SOFTWARE WILL BE USED FOR BINARY
!                OUTPUT (see LGRIBOUT).
!     LFDBIOOUT : IF TRUE THE OUTPUT WILL BE SEND TO THE FDB.
!     LGRIBIN : IF TRUE GRIB SPECTRA 250 ARE INPUT FOR INTERPOLATION AND
!               CONVERSION TO PARAMETER 251. OTHERWISE A COLD START WILL
!               BE CREATED BASED ON THE CHOICE OF IOPTI.
!     LGRIBOUT : IF TRUE THE WAVE SPECTRA IS OUTPUT IN GRIB FORMAT, ELSE
!                THE BINARY RESTART FILES ARE PRODUCED.
!     MARSTYPE : DATA TYPE USED TO CODE DATA IN GRIB
!     YCLASS   : DATA CLASS USED TO CODE DATA IN GRIB.
!     YEXPVER : EXPERIMENT VERSION USED TO CODE DATA IN GRIB.
!     CFDB2DSP : FDB ROOT USED TO WRITE GRIB OUTPUT TO FDB.

! ----------------------------------------------------------------------


!*    0. SET DEFAULT VALUES FOR THE NAMELIST ELEMENTS.
!        ---------------------------------------------

      msub = 0
      csuvi = ' '
      ctest = ' '
      mtest = 0
      msuvi = 0

      HEADER = ZERO
      IOPTI  =    1
      ITEST  =   -9
      ITESTB =   -9
      ALFA   =    0.0
      FM     =    0.0
      GAMMA  =    0.0
      SA     =    0.0
      SB     =    0.0
      THETA  =    0.0
      FETCH  =    0.0
      USERID = ZERO
      RUNID  = ZERO
      PATH   = ZERO
      CPATH  = ZERO
      CDATEA = ZERO
      CLTUNIT= ZERO
      IDELWI =    0

      LMESSPASS = .FALSE.

      LPBIOIN   = .FALSE.
      LPBIOOUT  = .FALSE.
!     LFDBIOIN : obsolete, was used to read directly from fdb
!                it should not be used
      LFDBIOIN  = .FALSE. 
      LFDBIOOUT = .FALSE. 
      LGRIBIN  = .FALSE.
      LGRIBOUT = .FALSE.

      MARSTYPE = 'an'
      YCLASS   = 'rd'
      YEXPVER  = USERID//'a'
      CFDB2DSP = "/vpp700/fdb18"

!*    1. DEFINE UNIT NAMES.
!        ------------------

      IU04 = 4
      IU05 = 5
      IU06 = 6
      IU07 = 7
      IU12 = 12
      IU14 = 14
      IU15 = 15

!     1.1 INITIALISATION OF MPP VARIABLES WITH DEFAULT VALUES
!         ---------------------------------------------------

      IRANK=1
      NPROC=1
      X4=1.
      NPRECR = N_PRECISION(X4)
      I4=1
      NPRECI = N_PRECISION(I4)
      KTAG=100

! ----------------------------------------------------------------------

! ALLOCATE NECESSARY ARRAYS

      ALLOCATE (NSTART(NPROC),NEND(NPROC))

! ----------------------------------------------------------------------

!*    2. READ NAMELIST NALINE.
!        ---------------------

      CDATEWO = ' '
      CDAWIFL = ' '
      CDATEFL = ' '

      READ (IU05, NALINE)

      IF (CLTUNIT .EQ. 'H') IDELWI = IDELWI*3600
      CDATEF = CDATEA
      CDATER = '000000000000'
      CDATES = '000000000000'
      ICPLEN=LEN_TRIM(CPATH)

      NTOTENS = 0
! ----------------------------------------------------------------------

!*    3. READ PREPROC OUTPUT.
!        --------------------

      CALL READPRE (IU07, IREFRA)

      NINF=1
      NSUP=NIBLO

      NSTART=1
      NEND=IJL(1)

      ALLOCATE(IUSCR(NBLO))
      DO IG=1,NBLO
        IUSCR(IG) = 39+IG
      ENDDO

!*    3.* SET GRIB HEADERS FOR INPUTS/OUTPUTS
!         -----------------------------------
      IF(LGRIBIN.OR.LGRIBOUT) THEN
        CALL SETWGRIBHD
      ENDIF
! ----------------------------------------------------------------------

!*    4. PRINTER PROTOCOL OF INPUT.
!        --------------------------

      WRITE (IU06,'(" PROGRAM PRESET:", / , "
     & WAM Cycle4 Release8 CI 4 (common with IFS CY19R2)
     &    CC=wam_wab_CY19R1_modules", /,
     & "  USER INPUT:" )')

      WRITE (IU06,'(A70)') HEADER
      WRITE (IU06,'('' MODEL OPTIONS  :'',/)')

      IF(LGRIBIN) THEN
        WRITE (IU06,*)'  '
        WRITE (IU06,*)' ******************************************'
        WRITE (IU06,*)' THE INPUT SPECTRA ARE GRIB FILES '
        WRITE (IU06,*)' ******************************************'
        WRITE (IU06,*)'  '
      ENDIF
      WRITE (IU06,*)'  '
      WRITE (IU06,*)' ******************************************'
      IF(LGRIBOUT) THEN
        WRITE (IU06,*)' THE OUTPUT SPECTRA WILL BE GRIBBED '
        WRITE (IU06,*)' INTO ',NFRE*NANG,' FIELDS'
      ELSE
        WRITE (IU06,*)' THE OUTPUT OF SPECTRA WILL BE BINARY  '
      ENDIF
      WRITE (IU06,*)' ******************************************'
      WRITE (IU06,*)'  '

      IF (.NOT.LGRIBIN.AND.IOPTI.EQ.0) THEN
        WRITE (IU06,'('' INITIAL VALUES ARE COMPUTED FROM'',
     &   '' INPUT PARAMETERS.'')')
      ELSEIF (.NOT.LGRIBIN.AND.IOPTI.EQ.1) THEN
        WRITE (IU06,'('' INITIAL VALUES ARE COMPUTED FROM'',
     &   '' LOCAL WIND.'')')
        WRITE (IU06,'('' WAVE ENERGY IS ZERO IN CALM WIND AREAS.'')')
      ELSEIF (.NOT.LGRIBIN.AND.IOPTI.EQ.2) THEN
        WRITE (IU06,'('' INITIAL VALUES ARE COMPUTED FROM'',
     &   '' LOCAL WIND.'')')
        WRITE (IU06,'('' PARAMETERS USED IN CALM WIND AREAS.'')')
      ELSEIF (LGRIBIN.AND.IOPTI.EQ.1) THEN
        WRITE (IU06,'('' INITIAL SPECTRAL VALUES ARE COMPUTED FROM'',
     &   '' ARCHIVED SPECTRUM AND DRAG COEF. FIELDS.'')')
        WRITE (IU06,'('' AND ARCHIVED WAVE MODEL ANALYSED WIND SPEED'',
     &   '' IF AVAILABLE.'')')
      ELSEIF (LGRIBIN.AND.IOPTI.EQ.2) THEN
        WRITE (IU06,'('' INITIAL SPECTRAL VALUES ARE COMPUTED FROM'',
     &   '' ARCHIVED SLIT SPECTRUM AND DRAG COEF. FIELDS.'')')
        WRITE (IU06,'('' AND ARCHIVED WAVE MODEL ANALYSED WIND SPEED'',
     &   '' IF AVAILABLE.'')')
      ELSE
        WRITE (IU06,'('' INVALID INPUT OPTION. PROGRAM WILL ABORT '')')
        CALL ABORT1
      ENDIF

      WRITE(IU06,*) '  '
      IF(.NOT.LGRIBIN.AND.LPBIOIN) THEN
        WRITE (IU06,*)
     &   ' PBIO SOFTWARE IS USED TO READ SELECTED INPUT FILES' 
      ENDIF
      IF(.NOT.LGRIBOUT.AND.LPBIOOUT) THEN
        WRITE (IU06,*)
     &   ' PBIO SOFTWARE IS USED TO WRITE TO SELECTED OUTPUT FILES' 
      ENDIF
      WRITE(IU06,*) '  '

      IF(LFDBIOIN) THEN
        WRITE (IU06,*)
     &   ' FDB SOFTWARE IS USED TO READ SELECTED INPUT FILES' 
      ENDIF
      IF(LFDBIOOUT) THEN
        WRITE (IU06,*)
     &   ' FDB SOFTWARE IS USED TO WRITE TO SELECTED OUTPUT FILES' 
          WRITE(IU06,'("   ITS ROOT DIRECTORY IS ", A)')
     &     CFDB2DSP(1:LEN_TRIM(CFDB2DSP))
        IF(.NOT.LGRIBOUT) THEN
          WRITE (IU06,*) ' OUTPUT OF BINARY DATA ' 
          WRITE (IU06,*) ' HOWEVER THIS OPTION IS OBSOLETE AND '
          WRITE (IU06,*) ' SHOULD NOT BE RUN. PROGRAM WILL ABORT.'
          CALL ABORT1
        ENDIF
      ENDIF

      WRITE (IU06,*) ' TEST OUTPUT LEVEL IS .......... ITEST = ', ITEST
      WRITE (IU06,*) ' TEST OUTPUT IN BLOCK LOOP UPTO ITESTB = ', ITESTB

      IF (.NOT.LGRIBIN) THEN
        WRITE (IU06,'('' JONSWAP PARAMETERS  :'',/)')
        WRITE (IU06,'('' ALFA : '',F10.5,'' FM : '',F10.5,'' GAMMA : '',
     &                F10.5,'' SA : '',F10.5,'' SB : '',F10.5)')
     &                ALFA, FM, GAMMA, SA, SB
        WRITE (IU06,'('' MEAN WAVE DIRECTION :  THETA = '',F10.5,
     &                '' DEGREE'')')  THETA
      ENDIF
      WRITE (IU06,*) '  '
      WRITE (IU06,*) ' WIND INPUT TIMESTEP (SECONDS)      : ',IDELWI
      WRITE (IU06,*) '  '
      WRITE (IU06,*) ' END OF USER INPUT PROTOCOLL'
      WRITE (IU06,'(''  NUMBER OF DIRECTION BINS  NANG = '',I4)') NANG
      WRITE (IU06,'(''  NUMBER OF FREQUENCY BINS  NFRE = '',I4)') NFRE

      IF (LFDBIOOUT .OR. LFDBIOIN) THEN
        WRITE(IU06,*) ' ACCESS TO FIELDS DATA BASE REQUESTED '
        WRITE(IU06,*) ' DATA BASE WILL BE INITIALIZED '
        INITFDB_STAT = IINITFDB()
        IF (INITFDB_STAT .NE. 0) THEN
          WRITE(IU06,'("****************************************")')
          WRITE(IU06,'("*                                      *")')
          WRITE(IU06,'("*    FATAL ERROR IN SUB. USERIN        *")')
          WRITE(IU06,'("*    ==========================        *")')
          WRITE(IU06,'("*    INITIALIZATIO OF THE              *")')
          WRITE(IU06,'("*    FIELDS DATA BASE                  *")')
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

!*    5. PREPARE WINDFIELD.
!        ------------------

      CDA = ZERO
      CDTPRO  = CDATEA
      CDATEE  = CDATEA
      CDAWIFL = CDATEA
      IDELPRO = IDELWI
      IDELWO  = IDELWI

      LWCOU=.FALSE.
      ILEN=1

      IF (.NOT.LGRIBIN) THEN
        IF(.NOT.ALLOCATED(FL1)) ALLOCATE (FL1(0:NIBLO,NANG,NFRE))
      ENDIF
      IF (.NOT.LGRIBIN.OR..NOT.LGRIBOUT) THEN
        IF(.NOT.ALLOCATED(U10OLD)) ALLOCATE(U10OLD(NIBLO,NBLO))
        IF(.NOT.ALLOCATED(THWOLD)) ALLOCATE(THWOLD(NIBLO,NBLO))
        IF(.NOT.ALLOCATED(USOLD)) ALLOCATE(USOLD(NIBLO,NBLO))
        IF(.NOT.ALLOCATED(Z0OLD)) ALLOCATE(Z0OLD(NIBLO,NBLO))
        IF(.NOT.ALLOCATED(TAUW)) ALLOCATE(TAUW(NIBLO,NBLO))
        U10OLD = 0.
        THWOLD = 0.
        USOLD = 0.
        TAUW = 0.
        Z0OLD = 0.
      ENDIF

      IF (.NOT.LGRIBIN.AND.IOPTI.GT.0) THEN
        CALL PREWIND (U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     &   IUGB, IVGB, ILEN)
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. PREWIND DONE'

      ELSEIF (LGRIBIN) THEN

        IF(.NOT.LGRIBOUT) THEN
!         RECONSTRUCT THE BINARY STRESS AND WIND FILE FROM THE WINDS
!         AND THE DRAG COEFFICIENT FIELDS.
          CALL BUILDSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,NIBLO)
        ENDIF

!       GET SPECTRA FROM GRIBBED DATA.

        ISPLIT=IOPTI

        CALL GRSTNAME(CDTPRO,CDATEF,'GRS',ICPLEN,CPATH,SFILENAME)
        CALL GRSTNAME(CDTPRO,CDATEF,'IGS',ICPLEN,CPATH,ISFILENAME)

!       INTERPOLATE THE GLOBAL FIELD IF IT IS NECESSARY

        LINTGRS=.FALSE. 
        CALL INTGRS (IU06,SFILENAME,ISPLIT,CDTPRO,LINTGRS,ISFILENAME)

        IF(LINTGRS) SFILENAME=ISFILENAME 

        IF(.NOT.LGRIBOUT) THEN
!         THE SPECTRA WILL BE SAVED AS A BINARY FILE
!CCC      note !! FL1 is allocated in readgrs !!!!
          CALL READGRS(IU06,SFILENAME,ISPLIT,CDTPRO)
          IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. READGRS DONE'
        ELSE
!         SPECTRA WILL BE SAVED AS A GRIB DATA (parameter 251) 
          CALL SPLITGRS (IU06,SFILENAME,ISPLIT,CDTPRO,ITEST,LFDBIOOUT,
     &                   MARSTYPE, YCLASS, YEXPVER, CFDB2DSP,CLDOMAIN)
          IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. SPLITGRS DONE'
        ENDIF
      ELSE
        IF (ITEST.GT.0) WRITE (IU06,*) ' WIND SET TO ZERO'
      ENDIF

! ----------------------------------------------------------------------

!*    6. DEFINE FETCH AND MAXIMUM PEAK FREQUENCY.
!        ----------------------------------------


      IF (.NOT.LGRIBIN) THEN
        IF (FETCH.LT.0.1E-5) FETCH = 0.5*DELPHI
        FRMAX = FM
        IF (IOPTI.NE.0) THEN
          WRITE (IU06,*) ' FETCH USED (METRES)       : ', FETCH
          WRITE (IU06,*) ' MAXIMUM PEAK FREQUENCY IS : ', FRMAX
        ENDIF
      ENDIF

! ----------------------------------------------------------------------

!*    7. GENERATE AND WRITE START FILES.
!        -------------------------------

      IF (.NOT.LGRIBIN) THEN
        THETAQ = THETA * RAD
        CALL MSTART (IU12, IU14, IU15, IOPTI, FETCH, FRMAX,
     &              FL1,U10OLD,THWOLD,USOLD,Z0OLD,TAUW)
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. MSTART DONE'
      ENDIF

! ----------------------------------------------------------------------

!*    8. DISPOSE START FILES.
!        --------------------

!     STRESS RELATED FIELDS :
      IF (.NOT.LGRIBOUT) THEN
        CALL SAVSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,NSTART,NEND,
     &                 CDATEA,CDATEA)
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. SAVSTRESS DONE'
        IF (.NOT.LPBIOOUT)
     &     CALL GSFILE (IU06, IU15, 0, CDATEA, CDATEA, 'LAW', 'S')
      ENDIF

!     SPECTRA :
      IF (.NOT.LGRIBIN.AND.LGRIBOUT) THEN
!       THE COLD START SPECTRA WILL BE SAVED AS GRIB FILES.
        CALL SPEC2FDB (IU06, NPROC, IRANK, LWCOU, TH, FR,
     &                 FL1, NINF, NSUP, NANG, NFRE, NGX, NGY,
     &                 NLONRGG, KXLT, IXLG, NIBLO, NBLO,
     &                 AMONOP, XDELLA, NSTART, NEND,
     &                 NKSEK0, NKSEK1, NKSEK2, NKSEK3, NKSEK4,
     &                 NPSEK2, NPSEK3, KSEK0, KSEK1, KSEK2, KSEK3,
     &                 KSEK4, PSEK2, PSEK3, NTENCODE, 
     &                 DEG, PPMISS, PPEPS, PPREC,
     &                 CDATEA, CDATEF, CDATEA,
     &                 CLDOMAIN, MARSTYPE, YCLASS, YEXPVER,
     &                 CFDB2DSP, NENSFNB, NTOTENS, LFDBIOOUT,
     &                 ICPLEN, CPATH, ITEST, IFDB)
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. SPEC2FDB DONE'
      ELSEIF (.NOT.LGRIBOUT) THEN
        CALL SAVSPEC(FL1,NSTART,NEND,CDATEA,CDATEA,CDUM,
     &               MARSTYPE,YCLASS, YEXPVER, CFDB2DSP)
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. SAVSPEC DONE'
      ELSE
        IF (ITEST.GT.0) WRITE (IU06,*) ' OUTPUT TOOK PLACE IN SPLITGRS' 
      ENDIF

! ----------------------------------------------------------------------

!*    9. END OF JOB: DELETE WORK FILES.
!        ------------------------------

! DEALLOCATE ARRAYS

      IF(ALLOCATED(U10OLD)) DEALLOCATE(U10OLD)
      IF(ALLOCATED(THWOLD)) DEALLOCATE(THWOLD)
      IF(ALLOCATED(USOLD)) DEALLOCATE(USOLD)
      IF(ALLOCATED(Z0OLD)) DEALLOCATE(Z0OLD)
      IF(ALLOCATED(TAUW)) DEALLOCATE(TAUW)
      IF(ALLOCATED(FL1)) DEALLOCATE(FL1)
      IF(ALLOCATED(NSTART)) DEALLOCATE(NSTART)
      IF(ALLOCATED(NEND)) DEALLOCATE(NEND)

      WRITE (IU06,*) ' '
      WRITE (IU06,*) ' PROGRAM PRESET: ALL DONE'

      END
