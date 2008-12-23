      SUBROUTINE INITMDL (NADV)

! ----------------------------------------------------------------------

!**** *INITMDL* - INITIALIZES THE WAM MODEL.

!     L. ZAMBRESKY   GKSS/ECMWF    JULY 1988

!     MODIFIED BY:   H. GUNTHER    NOVEMBER 1989
!                    J. BIDLOT     FEBRUARY 1996-1997
!                    J. DOYLE      SEPTEMBER 1996
!                    J. BIDLOT     APRIL 1997
!                    J. BIDLOT     SEPTEMBER : MODIFY PARALLEL INPUT 
!                    S. ABDALLA    OCTOBER 2001: INCLUSION OF AIR 
!                                                DENSITY AND Zi/L

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

!       INITIALIZE THE WAM MODEL.

!**   INTERFACE.
!     ----------

!          ---- FORMAL PARAMETERS ----

!    *CALL* *INITMDL (NADV)*
!       *NADV*     INTEGER   NUMBER OF ADVECTION ITERATIONS
!                            PER CALL OF WAMODEL, OUTPUT PARAMETER.

!          ---- INPUT/OUTPUT UNITS ---

!          THE NAMES ARE DEFINED IN SECTION 1. OF THIS PROGRAM,
!          IF IT IS NOT MENTIONED OTHERWISE.

!           *IU01*   - INPUT  UNIT UNBLOCKED WIND FILE.
!                      (SEE SUB READWIND).
!           *IU02*   - INPUT  UNIT OF BOUDARY VALUES FROM A PREVIOUS
!                      COARSE GRID IF THIS A FINE GRID RUN.
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'FBI'
!                      (OUTPUT OF BOUINT).
!           *IU05*   - USER INPUT UNIT.
!           *IU06*   - PRINTER OUTPUT.
!           *IU07*   - INPUT  UNIT OF PRECOMPUTED GRID PARAMETERS.
!                      (OUTPUT OF PREPROC).
!           *IU08*   - INPUT  UNIT OF MODULE YOWUBUF.
!                      (OUTPUT OF PREPROC).
!NEST
!           *IU09*   - INPUT  UNIT MODULE YOWCBOU (OUTPUT OF PREPROC).
!           *IU10*   - INPUT  UNIT MODULE YOWFBOU (OUTPUT OF PREPROC).
!NEST
!           *IU11*   - INPUT UNIT OF SPECTRA AT ALL GRID POINTS.
!                      EACH PROPAGATION STEP THE FILES CONNECTED
!                      TO IU11 AND IU12 ARE INTERCHANGED.
!           *IU12*   - OUTPUT UNIT BLOCKS OF SPECTRA (SEE IU11).
!C                     THIS UNIT IS NOT USED ANYMORE IF PBIO SOFTWARE IS
!C                     USED TO READ THE SPECTRA. (FLAG LPBIOIN)
!           *IU13*   - INPUT UNIT OF SPECTRA ON LAST LATITUDE
!                      OF A BLOCK. SPECTRA ARE SAVED FROM THE
!                      SECOND LATITUDE OF NEXT BLOCK.
!                      EACH PROPAGATION STEP THE FILES CONNECTED
!                      TO IU13 AND IU14 ARE INTERCHANGED.
!           *IU14*   - OUTPUT UNIT SECOND LATITUDES (SEE IU13).
!           *IU15*   - OUTPUT UNIT LAST WINDFIELDS.
!C                     THIS UNIT IS NOT USED ANYMORE IF PBIO SOFTWARE IS
!C                     USED TO READ THE LAST WINDFIELDS. (FLAG LPBIOIN)
!REFRA
!           *IU16*   - INPUT/OUTPUT UNIT OF MODULE YOWREFD.
!REFRA
!           *IU17*   - INPUT  UNIT OF BLOCKED WINDS.
!                      THIS FILE IS DYNAMICALLY ASSIGNED IN SUB
!                      IMPLSCH. SEE SUB CREWFN TOO.
!           *IU18*   - INPUT  UNIT OF BLOCKED WINDS.
!                      THIS FILE IS DYNAMICALLY ASSIGNED IN SUB
!                      IMPLSCH. SEE SUB CREWFN TOO.
!NEST
!           *IU19*   - OUTPUT UNIT OF BOUNDARY VALUES IF
!                      THIS IS A FINE GRID RUN.
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'CBO'
!NEST
!           *IU20*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
!                      OF THE TOTAL SPECTRUM (FIRST GUESS).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'MAP'
!           *IU21*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
!                      OF SWELL AND WIND WAVES (FIRST GUESS).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'SWE'
!ASSI
!           *IU22*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
!                      OF THE TOTAL SPECTRUM (ANALYSED).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'AMP'
!           *IU23*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
!                      OF SWELL AND WIND WAVES (ANALYSED)
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'ASS'
!ASSI
!           *IU25*   - OUTPUT UNIT OF SPECTRA AT CERTAIN GRID POINTS
!                      (FIRST GUESS).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'OUT'
!           *IU26*   - OUTPUT UNIT OF SWELL SPECTRA AT CERTAIN POINTS
!                      (FIRST GUESS).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'SWS'
!ASSI
!           *IU27*   - OUTPUT UNIT OF SPECTRA AT CERTAIN GRID POINTS
!                      (ANALYSED).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'AUT'
!           *IU28*   - OUTPUT UNIT OF SWELL SPECTRA AT CERTAIN POINTS
!                      (ANALYSED).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'ASW'
!ASSI
!           *IUSCR*  - SCRATCH UNITS FOR WIND PROCESSING.
!           *IUVELO* - OUTPUT UNIT OF BLOCKED WINDFILEDS.
!                      FILES ARE DYNAMICALLY ASSIGNED IN SUB
!                      NOTIM OR TIMIN. SEE SUB CREWFN TOO.

!___PACK  THESE UNITS ARE USED FOR FILES CONTAINING GRIB DATA.

!         (IN CASE OF SPECTRA ALL SPECTRA OF THE MODEL ARE WRITTEN)
!           *IU30*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
!                      OF THE TOTAL SPECTRUM (FIRST GUESS).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'MPP'
!           *IU31*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
!                      OF SWELL AND WIND WAVES (FIRST GUESS).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'SWP'
!ASSI
!           *IU32*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
!                      OF THE TOTAL SPECTRUM (ANALYSED).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'APP'
!           *IU33*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
!                      OF SWELL AND WIND WAVES (ANALYSED)
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'AWP'
!ASSI
!           *IU35*   - OUTPUT UNIT OF SPECTRA AT ALL GRID POINTS
!                      (FIRST GUESS).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'OUP'
!           *IU36*   - OUTPUT UNIT OF SWELL SPECTRA AT ALL POINTS
!                      (FIRST GUESS).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'WSP'
!ASSI
!           *IU37*   - OUTPUT UNIT OF SPECTRA AT ALL GRID POINTS
!                      (ANALYSED).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'AUP'
!           *IU38*   - OUTPUT UNIT OF SWELL SPECTRA AT ALL POINTS
!                      (ANALYSED).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'ASP'
!ASSI
!___PACK

!          FOR A START THE RESTART FILES WILL DYNAMICALLY BE ASSIGNED
!          AND COPIED TO THE *OUTPUT* UNITS (IU12, IU14, IU15).
!          IF IT IS REQUESTED TO SAVE RESTART FILES THESE WILL BE COPIED
!          IN REGULAR INTERVALS FROM THE UNIT ALIAS FILES OF IU11 OR
!          IU12, IU13 OR IU14, AND IU15 TO THE PREMANENT RESTART FILES.
!          FOR DETAILS OF THE FILE NAMES SEE SUB GSFILE.

!          SUB NOTIM OR TIMIN OPENS FILES AND ASSIGNS THEM TO UNIT
!          IUVELO FOR THE BLOCKED WINDS. THESE FILES ARE READ AND
!          DELETED IN SUB IMPLSCH (IU17 AND IU18). THE FILE NAMES
!          ARE CREATED IN SUB CREWFN.

!          THE PROGRAM CLOSES AND DELETES ALL WORK FILES.

!          ALL PARAMETERS HAVE TO BE THE VALUES GIVEN AT THE END
!          OF THE PREPROC OUTPUT IN COLUMN 'REQUIRED'.

!     METHOD.
!     -------

!          THIS ROUTINE INITIALISES THE WAVEMODEL:
!            -  DEFINES THE UNITS FOR INPUT/OUTPUT,
!            -  READS THE USER INPUT FILE,
!            -  INITIALIZES SOME MODEL CONSTANTS,
!            -  GETS THE RECOVERY FILES,
!            -  READS THE YOWMON BLOCKS PRECOMPUTED BY PROG PREPROC,
!            -  DOES SOME GENERAL BOOKEEPING REGARDING
!               DATES, INTEGRATION TIME STEPS AND OUTPUT TIME STEPS.
!            -  READ MODULE YOWUBUF AND SPECTRA IF ONE BLOCK VERSION.
!REFRA
!            -  PRECOMPUTES AND WRITES TO IU16 REFRACTION TERMS.
!REFRA
!            -  OPENS THE FIRST RESULT FILES.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.
!       *DIFDATE*   - COMPUTES A TIME DIFFERENCE.
!       *GETREST*   - GETS RESTART FILES.
!REFRA
!       *GRADI*     - COMPUTES DEPTH AND CURRENT GRADIENTS.
!REFRA
!       *GSFILE*    - ROUTINE TO DYNAMICALLY FETCH OR DISPOSE FILES.
!       *GETSPEC*   - READS SPECTRA
!       *GETSTRESS  - READS RESTART STRESS/WIND FILE
!       *HEADBC*    - WRITE BOUNDARY OUTPUT FILE HEADER.
!       *INCDATE*   - UPDATE DATE TIME GROUP.
!       *INIWCST*   - SETS GLOBAL CONSTANTS.
!REFRA
!       *PROPDOT*   - PRECOMPUTE REFRACTION.
!REFRA
!NEST
!       *READBOU*   - READS PREPROC BOUNDARY FILES.
!NEST
!       *SETMARSTYPE* SETS VARIABLE MARSTYPE
!       *SETWGRIBHD*- SETS DEFAULT GRIB HEADERS  
!       *USERIN*    - READS USER INPUT.

!     REFERENCE
!     ---------

!          A MORE DETAILED  DISCUSSION MAY BE FOUND IN SUB WAMODEL.

! ----------------------------------------------------------------------

      USE YOWCBOU  , ONLY : IBOUNC   ,NBOUNC   ,IJARC    ,IGARC
      USE YOWCOUP  , ONLY : LWCOU    ,KCOUSTEP
      USE YOWCOUT  , ONLY : COUTT    ,COUTLST  ,FFLAG20  ,FFLAG21  ,
     &            FFLAG25  ,FFLAG26  ,PFLAG20  ,PFLAG21  ,PFLAG25  ,
     &            PFLAG26  ,GFLAG20  ,GFLAG21  ,GFLAG25  ,GFLAG26  ,
     &            NGOUT    ,IGAR     ,IJAR     ,NOUTT    ,LOUTINT  ,
     &            LOUTSPP
      USE YOWCURR  , ONLY : CDTCUR   ,IDELCUR  ,CDATECURA
      USE YOWFBOU  , ONLY : IBOUNF
      USE YOWFRED  , ONLY : FR       ,TH       ,FR5     ,
     &             FRM5    ,COFRM4   ,COEF4    ,FRATIO 
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,IGL      ,IJLT     ,
     &             COSPH
      USE YOWICE   , ONLY : ICEMASK
      USE YOWINDN  , ONLY : ENH
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,AMOWEP   ,AMOSOP   ,
     &            AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,
     &            KMNOP    ,KMSOP
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,
     &            KTAG
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGY      ,NBLO     ,
     &            NIBLO    ,NIBLD    ,NBLD     ,NIBLC    ,NBLC
      USE YOWPCONS , ONLY : G        ,CIRC     ,ZPI      ,RAD
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTH    ,DEPTHA   ,DEPTHD   ,
     &            INDEP    ,TCGOND
      USE YOWSPEC  , ONLY : NSTART   ,NEND     ,KLENTOP  ,KLENBOT  ,
     &            U10OLD   ,THWOLD   ,USOLD    ,Z0OLD    ,TAUW     ,
     &            ROAIRO   ,ZIDLOLD  ,
     &            FL1      ,SL
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDTPRO   ,CDTRES   ,
     &            CDTINTT  ,CDTINTS  ,CDTSPT   ,CDTSPS   ,IDELPRO  ,
     &            IDELT    ,IDELWI   ,IDELWO   ,IDELRES  ,IDELINT  ,
     &            IDELINS  ,IDELSPT  ,IDELSPS  ,ISHALLO  ,IREFRA   ,
     &            CDATEA   ,MARSTYPE ,LANAONLY ,ISNONLIN ,IPROPAGS
      USE YOWTABL  , ONLY : FAC0     ,FAC1     ,FAC2     ,FAC3     ,
     &            FAK      ,FRHF     ,DFIMHF
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWTEXT  , ONLY : LRESTARTED
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KLAT2    ,KCOR     ,
     &            KRLAT    ,KRLON    ,KRLAT2   ,KRLON2   ,LSAMEDEPTH
      USE YOWUNIT  , ONLY : IU02     ,IU11     ,IU12     ,
     &            IU13     ,IU14     ,IU15     ,IU16     ,IU17     ,
     &            IU18     ,IU19     ,IU20     ,IU21     ,IU22     ,
     &            IU23     ,IU25     ,IU26     ,IU27     ,IU28     ,
     &            IUVELO   ,IUSCR    ,IU30     ,IU31     ,
     &            IU32     ,IU33     ,IU35     ,IU36     ,IU37     ,
     &            IU38
      USE YOWWAMI  , ONLY : CBPLTDT
      USE YOWWIND  , ONLY : CDA      ,CDAWIFL  ,CDATEWO  ,CDATEFL
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! -------------------------------------------------------------------

      CHARACTER*14 ZERO,CDUM
      CHARACTER*24 FILNM

!----------------------------------------------------------------------
      REAL ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('INITMDL',0,ZHOOK_HANDLE)

      CALL INIWCST


!     0.1 SET IREAD TO THE PROCESSOR WHICH WILL READ THE INPUT FILES
!         ----------------------------------------------------------
      IREAD=1

      IF(NPROC.EQ.1) IREAD=1

!*    1. DEFINITION OF MODEL PARAMETERS.
!        -------------------------------

      ZERO = ' '

!*    1.1  DEFINE UNIT NAMES.
!          ------------------

!NEST
      IU02 = 2
!NEST


! STANDARD INPUT UNIT (now a set of namelists)
! --------------------

      IU05 =  I_GET_UNIT (IU06, 'wam_namelist', 's', 'f', 0)

! STANDARD OUTPUT UNIT
! --------------------
!    DEFINITION HAS  BEEN MOVED TO MPDECOMP

! GRID AND UBUF FILES UNITS
! -------------------------

!    DEFINITION HAS  BEEN MOVED TO MPDECOMP

!NEST
      IU09 = 9
      IU10 = 10
!NEST
      IU11 = 11


! RESTART FILES UNITS
! -------------------

! NOTE IF MESSAGE PASSING
! NOTE IU12 WILL ONLY BE CONNECTED TO UNIT 12 FOR PURE BINARY INPUT
! NOTE OTHERWISE IT WILL CONNECTED TO THE FILENAME BLS WITH
! NOTE AN EXTENSION FUNCTION OF THE PE.

      IF(LWCOU) THEN
        IU12 = 36
        IU15 = 37
      ELSE
        IU12 = 12
        IU15 = 15
      ENDIF


      IU13 = 13
      IU14 = 14
!REFRA
      IU16 = 16
!REFRA
      IU17 = 17
      IU18 = 18
!NEST
      IU19 = 19
!NEST
      IU20 = 20
      IU21 = 21
!ASSI
      IU22 = 22
      IU23 = 23
!ASSI
      IU25 = 25
      IU26 = 26
!ASSI
      IU27 = 27
      IU28 = 28
!ASSI

!__PACK
      IU30 = 30
      IU31 = 31
!ASSI
      IU32 = 32
      IU33 = 33
!ASSI
      IU35 = 35
      IU36 = 36
!ASSI
      IU37 = 37
      IU38 = 38
!ASSI
!__PACK

      IUVELO = 39

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. INITMDL: UNITS DEFINED'
        CALL FLUSH (IU06)
      ENDIF

      LOUTINT=.FALSE.
      LOUTSPP=.FALSE.

!*    1.2 INPUT OF USER PARAMETER.
!         ------------------------

      CALL USERIN (IU05, IBOUNC, IBOUNF, IFORCA)
      CLOSE(IU05)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '    SUB. INITMDL: USERIN DONE'
        CALL FLUSH(IU06)
      ENDIF

!     DEFINE A FEW CONSTANTS FOR USE IN IMPLSCH

      IF(.NOT.ALLOCATED(FR5)) ALLOCATE(FR5(NFRE))
      IF(.NOT.ALLOCATED(FRM5)) ALLOCATE(FRM5(NFRE))
      IF(.NOT.ALLOCATED(COFRM4)) ALLOCATE(COFRM4(NFRE))

      DO M=1,NFRE
        FR5(M) = FR(M)**5
        FRM5(M) = 1./FR5(M)
        IF(ISHALLO.EQ.1) THEN
          COFRM4(M) = COEF4*G/FR(M)**4
        ELSE
          COFRM4(M) = 0.5*COEF4*ZPI**4/SQRT(G)
        ENDIF
      ENDDO

      IF(LMESSPASS) THEN

        KTAG=100

!      1.3 OUTPUT MODEL DECOMPOSITION DETAILS 
!          ----------------------------------

        IF (ITEST.GE.1) THEN
          WRITE(IU06,*)  ' MODEL DOMAIN DECOMPOSITION  : '
          WRITE(IU06,*)  ' =========================='
          DO IP=1,NPROC
            WRITE(IU06,*)
            WRITE(IU06,*) ' PROCESS NUMBER : ',IP
            WRITE(IU06,*) ' NSTART : ',NSTART(IP)
            WRITE(IU06,*) ' NEND   : ',NEND(IP)
            WRITE(IU06,*) ' N      : ',NEND(IP)-NSTART(IP)+1
            WRITE(IU06,*) ' KLENBOT: ',KLENBOT(IP)
            WRITE(IU06,*) ' KLENTOP: ',KLENTOP(IP)
            WRITE(IU06,*) ' ----------------------- '
          ENDDO
          CALL FLUSH(IU06)
        ENDIF
      ENDIF


!     1.4 SET INTEGRATED OUTPUT PARAMETER TABLE FOR MESSAGE PASSING
!         IN CASE OF NO MESSAGE PASSING THE TABLE IS STILL USED TO
!         POINT TO THE PARAMETRE WHICH ARE OUTPUT.
!         ----------------------------------------------------------

      CALL MPCRTBL
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '    SUB. INITMDL: MPCRTBL DONE'
        CALL FLUSH (IU06)
      ENDIF

!     1.5 DETERMINE LAST OUTPUT DATE
!         --------------------------

      IF (NOUTT.GT.0) THEN
        COUTLST=COUTT(NOUTT)
      ELSE
        COUTLST=CDATEA
        IF (GFLAG20) THEN
          DO WHILE (COUTLST.LE.CDATEE.AND.IDELINT.GT.0)
            CALL INCDATE (COUTLST,IDELINT)
          ENDDO
          CALL INCDATE (COUTLST,-IDELINT)
        ENDIF

        CDUM=CDATEA
        IF (GFLAG21) THEN
          DO WHILE (CDUM.LE.CDATEE.AND.IDELINS.GT.0)
            CALL INCDATE (CDUM,IDELINS)
          ENDDO
          CALL INCDATE (CDUM,-IDELINS)
        ENDIF
        IF(CDUM.GT.COUTLST) COUTLST=CDUM

        IF(COUTLST.GT.CDATEE) COUTLST=CDATEE
      ENDIF

! ----------------------------------------------------------------------

!*    2. INPUT FROM PREPROCESSING PROGRAMS.
!        ----------------------------------

!     THE ACTUAL READING HAS BEEN MOVED TO MPDECOMP.

      ITESTB = MIN(ITESTB,IGL)

      IF (NGOUT.EQ.0.AND.(PFLAG25.OR.FFLAG25.OR.PFLAG26.OR.FFLAG26))THEN
        WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*) '+                                           +'
        WRITE(IU06,*) '+    WARNING ERROR IN SUB. INITMDL          +'
        WRITE(IU06,*) '+    =============================          +'
        WRITE(IU06,*) '+ NUMBER OF OUTPUT POINT FOR SPECTRA        +'
        WRITE(IU06,*) '+ READ FROM PREPROC OUTPUT IS ZERO.         +'
        WRITE(IU06,*) '+ OUTPUT FLAGS CHANGED TO NO OUTPUT.        +'
        WRITE(IU06,*) '+                                           +'
        WRITE(IU06,*) '+ MODEL CONTINUES.                          +'
        WRITE(IU06,*) '+                                           +'
        WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++'
        PFLAG25 = .FALSE.
        FFLAG25  = .FALSE.
        PFLAG26 = .FALSE.
        FFLAG26  = .FALSE.
      ENDIF
!NEST

      IF (.NOT.ALLOCATED(IUSCR)) THEN
        ALLOCATE(IUSCR(NBLO))
      ENDIF
      DO IG=1,NBLO
        IUSCR(IG) = 39+IG
      ENDDO

!*    2.1 READ MODULE YOWCBOU AND YOWFBOU.
!     ------------------------------------

      IF (IBOUNC.EQ.1 .OR. IBOUNF.EQ.1) THEN
        CALL READBOU (IU09, IU10, IU06)
        CLOSE (UNIT=IU09, STATUS='KEEP')
        CLOSE (UNIT=IU10, STATUS='KEEP')
        IF (ITEST.GE.2) THEN
          WRITE(IU06,*)'    SUB. INITMDL: BOUNDARY POINTS READ ',
     &     ' AND FILES CLOSED '
          CALL FLUSH (IU06)
        ENDIF
      ENDIF
!NEST

!*    2.2.* SET GRIB HEADERS FOR INPUTS/OUTPUTS
!          ------------------------------------
      LANAONLY=.FALSE.
      IF((CDATEA.EQ.CDATEE).AND.(CDATEA.EQ.CDATEF)) LANAONLY=.TRUE.

      CALL SETMARSTYPE

      CALL SETWGRIBHD
      IF (ITEST.GE.2)
     1   WRITE(IU06,*) '   SUB. INITMDL: DEFAULT  GRIB HEADERS SET'

      IF(MARSTYPE.EQ.'cf'.OR.MARSTYPE.EQ.'pf'.OR.MARSTYPE.EQ.'fc')THEN
        IF(ALLOCATED(FAC0)) DEALLOCATE(FAC0)
        IF(ALLOCATED(FAC1)) DEALLOCATE(FAC1)
        IF(ALLOCATED(FAC2)) DEALLOCATE(FAC2)
        IF(ALLOCATED(FAC3)) DEALLOCATE(FAC3)
        IF(ALLOCATED(FAK)) DEALLOCATE(FAK)
        IF(ALLOCATED(FRHF)) DEALLOCATE(FRHF)
        IF(ALLOCATED(DFIMHF)) DEALLOCATE(DFIMHF)
      ENDIF

! ----------------------------------------------------------------------

!*    3. PRINT INITIAL CONDITIONS AS READ FROM PERPROCESSING.
!        ----------------------------------------------------

      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' WAVE MODEL GRID ORGANISATION:'
      WRITE(IU06,3002) ' SOUTHERNMOST LATITUDE IN GRID IS .......: ',
     & AMOSOP, ' DEGREE'
      WRITE(IU06,3002) ' NORTHERNMOST LATITUDE IN GRID IS .......: ',
     & AMONOP, ' DEGREE'
      WRITE(IU06,3002) ' WESTERNMOST LONGITUDE IN GRID IS .......: ',
     & AMOWEP, ' DEGREE'
      WRITE(IU06,3002) ' EASTERNMOST LONGITUDE IN GRID IS .......: ',
     & AMOEAP, ' DEGREE'
      WRITE(IU06,3002) ' LATITUDE INCREMENT IS ..................: ',
     & XDELLA, ' DEGREE'
      WRITE(IU06,3002) ' LONGITUDE INCREMENT IS .................: ',
     & XDELLO, ' DEGREE'
      WRITE(IU06,*) '  '
      WRITE(IU06,3003) ' TOTAL NUMBER OF BLOCKS IS...............: ',
     & IGL
      WRITE(IU06,3003) ' TOTAL LENGTH OF EACH BLOCK .............: ',
     & NIBLO
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' SPECTRAL RESOLUTION:'
      WRITE(IU06,3003) ' TOTAL NUMBER OF DIRECTIONS .............: ',
     & NANG 
      WRITE(IU06,3003) ' TOTAL NUMBER OF FREQUENCIES ............: ',
     & NFRE 
      IF(FFLAG25.OR.PFLAG25.OR.GFLAG25.OR.
     & FFLAG26.OR.PFLAG26.OR.GFLAG26) THEN
        WRITE(IU06,*) '  '
        WRITE(IU06,*) ' OUTPUT POINTS:'
        WRITE(IU06,*) ' TOTAL NUMBER OF OUTPUT POINTS IS: ',NGOUT
        WRITE(IU06,'(/,4X,''BLOCK NO'',6X,''INDEX NO'',
     &   8X,''LONGITUDE'',6X,''LATITUDE'')')
        DO I=1,NGOUT
          IX  = IXLG(IJAR(I),IGAR(I))
          KX  = KXLT(IJAR(I),IGAR(I))
          XLO = AMOWEP+REAL(IX-1)*ZDELLO(KX)
          XLA = AMOSOP+REAL(KX-1)*XDELLA

          WRITE(IU06,'((6X,I3,10X,I6,7X,F10.3,4X,F10.3))')
     &     IGAR(I), IJAR(I), XLO, XLA
        ENDDO
      ENDIF
 3002 FORMAT(3x,a,f9.3,a)
 3003 FORMAT(3x,a,i6,  a)
!NEST
      CALL FLUSH (IU06)
      IF (IBOUNC .EQ. 1) THEN
        WRITE(IU06,*)
        WRITE(IU06,*) ' COARSE GRID: BOUNDARY OUTPUT POINTS :'
        WRITE(IU06,*) ' TOTAL NUMBER OF BOUNDARY POINTS IS: ',NBOUNC
        IF (ITEST.GE.2) THEN
          WRITE(IU06,'(/,4X,''BLOCK NO'',6X,''INDEX NO'',
     &     8X,''LONGITUDE'',6X,''LATITUDE'')')
          DO I=1,NBOUNC
            IX    = IXLG(IJARC(I),IGARC(I))
            KX    = KXLT(IJARC(I),IGARC(I))
            XLO = AMOWEP+REAL(IX-1)*ZDELLO(KX)
            XLA = AMOSOP+REAL(KX-1)*XDELLA

            WRITE(IU06,'((6X,I3,10X,I5,7X,F10.3,4X,F10.3))')
     &       IGARC(I), IJARC(I), XLO, XLA
          ENDDO
        ENDIF
      ENDIF
!NEST

! ----------------------------------------------------------------------

!*    4. CONNECT RESTART FIELDS TO OUTPUT UNITS (IF PBIO SOFTWARE NOT
!        USED). AND READ IN LAST WINDFIELDS FROM RESTARTFILE.
!        ---------------------------------------------------------------

      IF (ITEST.GE.2)CALL FLUSH(IU06)

      CALL GETSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     &               ROAIRO,ZIDLOLD,ICEMASK,
     &               NSTART,NEND,IREAD)

      IF (LWCOU) THEN
        IDELWO = KCOUSTEP 
        IDELWI = KCOUSTEP 
      ENDIF
 
      CDA = CDTPRO

      IF (CDTPRO.NE.ZERO .OR. LRESTARTED) THEN

!*    4.1 MODEL STARTS FROM FILES OUT OF A PREVIOUS MODEL RUN.
!         ----------------------------------------------------

         CDATEWO = CDTPRO
         IF (IDELT.LT.IDELWO) CALL INCDATE(CDATEWO,IDELWO/2)
         CDAWIFL = CDTPRO
         IDELWH = MAX(IDELWI,IDELPRO)
         CALL INCDATE(CDAWIFL,IDELWH)
         CDATEFL = CDATEWO

        IF (CDTPRO.LT.CDATEA.OR.CDTPRO.GT.CDATEE.OR.
     &   (IFORCA.EQ.1.AND.CDTPRO.GT.CDATEF) .OR.
     &   (IFORCA.NE.1.AND.CDTPRO.LE.CDATEF)) THEN
          WRITE(IU06,*) ' *******************************************'
          WRITE(IU06,*) ' *    FATAL ERROR IN SUB. INITMDL          *'
          WRITE(IU06,*) ' *    ===========================          *'
          WRITE(IU06,*) ' * START DATE FROM RESTART FIELD IS NOT    *'
          WRITE(IU06,*) ' * MODEL PERIOD.                           *'
          IF (IFORCA.EQ.1) THEN
            WRITE(IU06,*) ' *  IN ANALYSIS PERIOD AS REQUESTED.       *'
          ELSE
            WRITE(IU06,*) ' *  IN FORECAST PERIOD AS REQUESTED.       *'
          ENDIF
          WRITE(IU06,*) ' * START DATE OF RUN       IS CDATEA = ',
     &     CDATEA
          WRITE(IU06,*) ' * START DATE OF FORECAST  IS CDATEF = ',
     &     CDATEF
          WRITE(IU06,*) ' * END   DATE OF RUN       IS CDATEE = ',
     &     CDATEE
          WRITE(IU06,*) ' * START DATE FROM RESTART IS CDTPRO = ',
     &     CDTPRO
          WRITE(IU06,*) ' *                                         *'
          WRITE(IU06,*) ' * PROGRAM ABORTS     PROGRAM ABORTS       *'
          WRITE(IU06,*) ' *                                         *'
          WRITE(IU06,*) ' *******************************************'
          CALL ABORT1
        ENDIF
      ELSE

!*    4.2 MODEL STARTS FROM FIELDS CREATED BY PRESET.
!         -------------------------------------------

        CDTPRO = CDATEA
        CDATEWO = CDATEA
        IF (IDELT.LT.IDELWO) CALL INCDATE(CDATEWO,IDELWO/2)
        CDAWIFL = CDATEA
        IDELWH = MAX(IDELWI,IDELPRO)
        CALL INCDATE(CDAWIFL,IDELWH)
        CDATEFL = CDATEWO
      ENDIF

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '    SUB. INITMDL: WIND FIELD AND ',
     &   ' COUNTER INITIALIZED'
        WRITE(IU06,*) '      NEXT WINDFIELD WILL BE READ AT     ',
     &   'CDATEWO = ',CDATEWO
        WRITE(IU06,*) '      NEXT WIND FILE WILL BE ACCESSED AT ',
     &   'CDATEFL = ',CDATEFL
        WRITE(IU06,*) '      NEXT WIND FILE NAME IS FROM        ',
     &   'CDAWIFL = ',CDAWIFL
      ENDIF

! ----------------------------------------------------------------------

!*    5. INITIALIZE MODEL TIME VARIABLES
!        -------------------------------

!*    5.1 OUTPUT TIME VARIABLES.
!         ----------------------

      IF (PFLAG20.OR.FFLAG20.OR.GFLAG20) THEN
        CDTINTT= CBPLTDT
        IF(LANAONLY) THEN
          CDTINTT=CDATEA
        ELSE
 5102     CONTINUE
          CALL INCDATE (CDTINTT, IDELINT)
          IF (CDTINTT.LE.CDTPRO) GOTO 5102
        ENDIF
      ELSE
        CDTINTT = ZERO
      ENDIF
      IF (PFLAG21.OR.FFLAG21.OR.GFLAG21) THEN
        CDTINTS= CBPLTDT
        IF(LANAONLY) THEN
          CDTINTS=CDATEA
        ELSE
 5103     CONTINUE
          CALL INCDATE (CDTINTS, IDELINS)
          IF (CDTINTS.LE.CDTPRO) GOTO 5103
        ENDIF
      ELSE
        CDTINTS = ZERO
      ENDIF
      IF (PFLAG25.OR.FFLAG25.OR.GFLAG25) THEN
        CDTSPT= CBPLTDT
        IF(LANAONLY) THEN
          CDTSPT=CDATEA
        ELSE
 5104     CONTINUE
          CALL INCDATE (CDTSPT, IDELSPT)
          IF (CDTSPT.LE.CDTPRO) GOTO 5104
        ENDIF
      ELSE
        CDTSPT = ZERO
      ENDIF
      IF (PFLAG26.OR.FFLAG26.OR.GFLAG26) THEN
        CDTSPS= CBPLTDT
        IF(LANAONLY) THEN
          CDTSPS=CDATEA
        ELSE
 5105     CONTINUE
          CALL INCDATE (CDTSPS, IDELSPS)
          IF (CDTSPS.LE.CDTPRO) GOTO 5105
        ENDIF
      ELSE
        CDTSPS = ZERO
      ENDIF

!*    5.2 FILE DISPOSE TIME AND RESTART TIME.
!         -----------------------------------

      CDTRES = CBPLTDT
      CALL INCDATE(CDTRES, IDELRES)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '    SUB. INITMDL: TIME COUNTER INITIALIZED'
        CALL FLUSH(IU06)
      ENDIF

! ----------------------------------------------------------------------

!*    6. CONSISTENCY CHECK ACCORDING TO CFL CRITERION.
!        ---------------------------------------------

      GVE = G/(ZPI*FR(1)*2.)
      DPH = DELPHI
      IF(DPH.EQ.0.0) THEN
        CFLP= 0.
      ELSE
        CFLP= IDELPRO*GVE/DPH
      ENDIF
      DLH = CIRC 
!     FIND THE SMALLEST DLH
      DO KX=KMSOP,KMNOP
        DLH_KX = DELLAM(KX)*COSPH(KX)
        DLH = MIN(DLH_KX,DLH)
      ENDDO

      IF(DLH.EQ.0.0) THEN
        CFLL= 0.
      ELSE
        CFLL= IDELPRO*GVE/DLH
      ENDIF
      IF (CFLP.GT.1. .OR. CFLL.GT.1.) THEN
        WRITE(IU06,*) ' **********************************************'
        WRITE(IU06,*) ' *                                            *'
        WRITE(IU06,*) ' *       FATAL ERROR IN SUB. INITMDL          *'
        WRITE(IU06,*) ' *       ===========================          *'
        WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED.               *'
        WRITE(IU06,*) ' * CFLP: ',CFLP,'  GROUP VELOCITY: ',GVE
        WRITE(IU06,*) ' * CFLL: ',CFLL,'  GROUP VELOCITY: ',GVE
        WRITE(IU06,*) ' * PROPAGATION TIME: ',IDELPRO
        WRITE(IU06,*) ' * GRID Y-DISTANCE AT MOST NORTH LATITUDE: ',DPH
        WRITE(IU06,*) ' * GRID X-DISTANCE AT MOST NORTH LATITUDE: ',DLH
        WRITE(IU06,*) ' *     PROGRAM ABORTS   PROGRAM ABORTS        *'
        WRITE(IU06,*) ' *                                            *'
        WRITE(IU06,*) ' **********************************************'
        CALL ABORT1
      ENDIF

! ----------------------------------------------------------------------

!*    7. NUMBER OF PROPAGATION TIME STEPS PER CALL.
!        ------------------------------------------

      NADV = IDELWI/IDELPRO
      NADV = MAX(NADV,1)
      IF(LANAONLY) THEN
        NADV = 0
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' **** NOTE ****'
        WRITE(IU06,*) ' NO FORWARD TIME INTEGRATION WILL BE PERFORMED'
        WRITE(IU06,*) ' ONLY THE DATA ASSIMILATION.'
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' '
      ENDIF

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '    SUB. INITMDL: NUMBER OF PROPAGATION STEPS'
        WRITE(IU06,*) '                  IN ONE CALL OF SUB WAVEMDL'
        WRITE(IU06,*) '                  WILL BE NADV = ', NADV
        CALL FLUSH(IU06)
      ENDIF

! ----------------------------------------------------------------------

!*    8. IF ONE BLOCK VERSION READ SPECTRA
!        ---------------------------------

      IF (IGL.EQ.1) THEN

!*    8.1 READ BLOCK OF SPECTRA.
!         ----------------------

        CALL GETSPEC(FL1,NSTART,NEND,IREAD,U10OLD,THWOLD)

        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '    SUB. INITMDL: READ',
     &     ' BLOCK OF SPECTRA'
          CALL FLUSH (IU06)
        ENDIF

!SHALLOW
!     8.2 COMPUTE THE NONLINEAR TRANSFER FUNCTION COEFFICIENTS FOR SNL
!         VALID FOR THE NEW FORMULATION
!         ------------------------------------------------------------ 

      IF(.NOT.ALLOCATED(ENH))
     &   ALLOCATE(ENH(NSTART(IRANK):NEND(IRANK),NFRE+4,IGL))

      IF (ISNONLIN.NE.0) THEN
        ENH_MAX=10.
        IF (ISHALLO.NE.1) THEN
          DO IG=1,IGL
            DO M=1,NFRE
               DO IJ=NSTART(IRANK),NEND(IRANK)
                 D = DEPTH(IJ,IG)
                 OM = ZPI*FR(M)
                 XK = AKI(OM,D)
                 ENH(IJ,M,IG) = MIN(ENH_MAX,TRANSF(XK,D))
               ENDDO
            ENDDO
            DO M=NFRE+1,NFRE+4
               DO IJ=NSTART(IRANK),NEND(IRANK)
                 D = DEPTH(IJ,IG)
                 OM = ZPI*FR(NFRE)*FRATIO**(M-NFRE)
!                NOTE THAT TFAK IS NOT DEFINED BEYOND M=NFRE
!                HENCE THE USE OF FUNCTIOn AKI.
                 XK = AKI(OM,D)
                 ENH(IJ,M,IG) = MIN(ENH_MAX,TRANSF(XK,D))
               ENDDO
            ENDDO
          ENDDO
        ELSE
          DO IG=1,IGL
            DO M=1,NFRE+4
               DO IJ=NSTART(IRANK),NEND(IRANK)
                 ENH(IJ,M,IG) = 1.
               ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!*    8.3 COMPUTE SHALLOW WATER TABLE INDICES.
!         ------------------------------------

        IF (ISHALLO.NE.1) THEN

          IF(.NOT.ALLOCATED(INDEP)) ALLOCATE(INDEP(NINF-1:IJLT(IGL)))
          INDEP(NINF-1)=NDEPTH
          DO IJ=NINF,IJLT(IGL)
            XD = LOG(DEPTH(IJ,IGL)/DEPTHA)/LOG(DEPTHD)+1.
            ID = NINT(XD)
            ID = MAX(ID,1)
            INDEP(IJ) = MIN(ID,NDEPTH)
          ENDDO

!     8.4 DETERMINE WHETHER A GRID POINT IS SURROUNDED BY
!         POINTS WITH THE SAME DEPTH INDEX (EXCLUDING LAND POINTS).
!         ---------------------------------------------------------

          IF(.NOT.ALLOCATED(LSAMEDEPTH)) 
     &        ALLOCATE(LSAMEDEPTH(NSTART(IRANK):NEND(IRANK)))

          IF(IPROPAGS.EQ.2) THEN
            DO IJ = NSTART(IRANK),NEND(IRANK)
              IF( INDEP(IJ).EQ.INDEP(KLON(IJ,1)) .AND.
     &            INDEP(IJ).EQ.INDEP(KLON(IJ,2)) .AND.
     &            INDEP(IJ).EQ.INDEP(KLAT(IJ,1)) .AND.
     &            INDEP(IJ).EQ.INDEP(KLAT(IJ,2)) .AND.
     &            INDEP(IJ).EQ.INDEP(KLAT2(IJ,1)).AND.
     &            INDEP(IJ).EQ.INDEP(KLAT2(IJ,2)) ) THEN
                LSAMEDEPTH(IJ) = .TRUE.
              ELSE
                LSAMEDEPTH(IJ) = .FALSE.
              ENDIF
              IF(LSAMEDEPTH(IJ)) THEN
                OUTER : DO ICL=-1,1
                  DO IC=1,4
                    IF( INDEP(IJ).NE.INDEP(KCOR(IJ,IC,ICL))) THEN
                      LSAMEDEPTH(IJ) = .FALSE.
                      EXIT OUTER
                    ENDIF 
                  ENDDO
                ENDDO OUTER
              ENDIF

            ENDDO
          ELSE IF(IPROPAGS.EQ.1) THEN
            DO IJ = NSTART(IRANK),NEND(IRANK)
              IF( INDEP(IJ).EQ.INDEP(KLON(IJ,1)) .AND.
     &            INDEP(IJ).EQ.INDEP(KLON(IJ,2)) .AND.
     &            INDEP(IJ).EQ.INDEP(KLAT(IJ,1)) .AND.
     &            INDEP(IJ).EQ.INDEP(KLAT(IJ,2)) .AND.
     &            INDEP(IJ).EQ.INDEP(KLAT2(IJ,1)).AND.
     &            INDEP(IJ).EQ.INDEP(KLAT2(IJ,2)).AND.
     &            INDEP(IJ).EQ.INDEP(KRLON(IJ,1)) .AND.
     &            INDEP(IJ).EQ.INDEP(KRLON(IJ,2)) .AND.
     &            INDEP(IJ).EQ.INDEP(KRLAT(IJ,1)) .AND.
     &            INDEP(IJ).EQ.INDEP(KRLAT(IJ,2)) .AND.
     &            INDEP(IJ).EQ.INDEP(KRLAT2(IJ,1)).AND.
     &            INDEP(IJ).EQ.INDEP(KRLAT2(IJ,2)).AND.
     &            INDEP(IJ).EQ.INDEP(KRLON2(IJ,1)).AND.
     &            INDEP(IJ).EQ.INDEP(KRLON2(IJ,2)) ) THEN
                LSAMEDEPTH(IJ) = .TRUE.
              ELSE
                LSAMEDEPTH(IJ) = .FALSE.
              ENDIF
            ENDDO
          ELSE
            DO IJ = NSTART(IRANK),NEND(IRANK)
              IF( INDEP(IJ).EQ.INDEP(KLON(IJ,1)) .AND.
     &            INDEP(IJ).EQ.INDEP(KLON(IJ,2)) .AND.
     &            INDEP(IJ).EQ.INDEP(KLAT(IJ,1)) .AND.
     &            INDEP(IJ).EQ.INDEP(KLAT(IJ,2)) .AND.
     &            INDEP(IJ).EQ.INDEP(KLAT2(IJ,1)).AND.
     &            INDEP(IJ).EQ.INDEP(KLAT2(IJ,2)) ) THEN
                LSAMEDEPTH(IJ) = .TRUE.
              ELSE
                LSAMEDEPTH(IJ) = .FALSE.
              ENDIF
            ENDDO
          ENDIF

        ENDIF
!SHALLOW
        IF (ITEST.GE.2)
     &   WRITE(IU06,*) '    SUB. INITMDL: MODULE YOWUBUF ',
     &   'READ AND INITIALIZED '
      ENDIF

! ----------------------------------------------------------------------
!REFRA
!*    9. PRECOMPUTE BOTTOM REFRACTION TERMS.
!        -----------------------------------
      IF (IREFRA.NE.0) THEN
        NIBLD=NIBLO
        NBLD=NBLO
        NIBLC=NIBLO
        NBLC=NBLO
      ELSE
        NIBLD=0
        NBLD=0
        NIBLC=0
        NBLC=0
      ENDIF

!     INITIALISE CDTCUR
      CDTCUR=CDATECURA
      IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
        CALL INCDATE(CDTCUR,-3600*IDELCUR)
      ENDIF

!     COMPUTE BOTTOM REFRACTION TERMS
      IF (IREFRA.EQ.1) THEN
        CALL PROPDOT
        IF (ITEST.GE.2)
     &   WRITE(IU06,*) '    SUB. INITMDL: REFRACTION TERMS',
     &   'INITIALIZED '
         CALL FLUSH(IU06)
      ENDIF
!     THE CURRENT REFRACTION TERMS ARE NOW COMPUTED IN WAVEMDL ONCE
!     THE CURRENT FIELD IS INPUT.

!REFRA
! ----------------------------------------------------------------------
!CMFA
!*    10. OPEN THE FIRST RESULT FILES.
!         ----------------------------

!           INTEGRATED PARAMETER DATA FILE.

!     IF (FFLAG20)
!    1   CALL GSFILE (IU06, IU20, 0, CDTPRO, CDATEF, 'MAP', 'S')

!            INTEGRATED SWELL AND SEA DATA FILE.

!     IF (FFLAG21)
!    1   CALL GSFILE (IU06, IU21, 0, CDTPRO, CDATEF, 'SWE', 'S')

!             SPECTRA FILE.

!     IF (FFLAG25)
!    1   CALL GSFILE (IU06, IU25, 0, CDTPRO, CDATEF, 'OUT', 'S')

!             SWELL SPECTRA FILE.

!     IF (FFLAG26)
!    1   CALL GSFILE (IU06, IU26, 0, CDTPRO, CDATEF, 'SWS', 'S')
!NEST
!             BOUNDARY VALUE FILE.

!     IF (IBOUNC.EQ.1)
!    1   CALL GSFILE (IU06, IU19, 0, CDTPRO, CDATEF, 'CBO', 'S')
!NEST

!     IF (ITEST.GE.2)
!    1   WRITE(IU06,*) '   SUB. INITMDL: OUTPUT FILES OPENED'
!CMFE
!NEST

!        WRITE BOUNDARY VALUE FILE HEADER.

      IF (IBOUNC.EQ.1) THEN
        CALL HEADBC (NBOUNC, IDELPRO, TH(1), FR(1), IU19, IU06)
        IF (ITEST.GE.2)
     &   WRITE(IU06,'(''    SUB. INITMDL: HEADER FOR  '',
     &   ''COARSE GRID WAS WRITTEN OF UNIT IU19 = '',A8)') IU19
      ENDIF
!NEST

      IF (ITEST.GE.2)
     1   WRITE(IU06,*) '   SUB. INITMDL: DEFAULT  GRIB HEADERS SET'

      IF (LHOOK) CALL DR_HOOK('INITMDL',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE INITMDL
