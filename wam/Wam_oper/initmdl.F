      SUBROUTINE INITMDL (KDELWI,NADV)
C
C ----------------------------------------------------------------------
C
C**** *INITMDL* - INITIALIZES THE WAM MODEL.
C
C     L. ZAMBRESKY   GKSS/ECMWF    JULY 1988
C
C     MODIFIED BY:   H. GUNTHER    NOVEMBER 1989
C                    J. BIDLOT     FEBRUARY 1996-1997
C                    J. DOYLE      SEPTEMBER 1996
C                    J. BIDLOT     APRIL 1997
C                    J. BIDLOT     SEPTEMBER : MODIFY PARALLEL INPUT 
C
CSHALLOW
C          DIFFERENCES FOR SHALLOW WATER RUNS TO DEEP WATER RUNS
C          ARE ENCLOSED IN COMMENT LINES : 'CSHALLOW'.
CSHALLOW
CNEST
C          DIFFERENCES FOR NESTED GRID RUNS TO NORMAL RUNS
C          ARE ENCLOSED IN COMMENT LINES : 'CNEST'.
CNEST
CREFRA
C          DIFFERENCES FOR REFRACTION RUNS TO NORMAL RUNS
C          ARE ENCLOSED IN COMMENT LINES : 'CREFRA'.
CREFRA
C
C*    PURPOSE.
C     --------
C
C       INITIALIZE THE WAM MODEL.
C
C**   INTERFACE.
C     ----------
C
C          ---- FORMAL PARAMETERS ----
C
C    *CALL* *INITMDL (KDELWI,NADV)
C      *NADV*      INTEGER   NUMBER OF ADVECTION ITERATIONS
C                            PER CALL OF WAMODEL, OUTPUT PARAMETER.
C      *KDELWI*              WIND INPUT TIME STEP IN SECONDS IF
C                            COUPLED WITH THE IFS.
C
C          ---- INPUT/OUTPUT UNITS ---
C
C          THE NAMES ARE DEFINED IN SECTION 1. OF THIS PROGRAM,
C          IF IT IS NOT MENTIONED OTHERWISE.
C
C           *IU01*   - INPUT  UNIT UNBLOCKED WIND FILE.
C                      (SEE SUB READWND).
C           *IU02*   - INPUT  UNIT OF BOUDARY VALUES FROM A PREVIOUS
C                      COARSE GRID IF THIS A FINE GRID RUN.
C                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'FBI'
C                      (OUTPUT OF BOUINT).
C           *IU04*   - INPUT AND OUTPUT UNIT FOR WAMINFO FILE
C           *IU05*   - USER INPUT UNIT.
C           *IU06*   - PRINTER OUTPUT.
C           *IU07*   - INPUT  UNIT OF PRECOMPUTED COMMON BLOCKS.
C                      (OUTPUT OF PREPROC).
C           *IU08*   - INPUT  UNIT OF COMMON UBUF.
C                      (OUTPUT OF PREPROC).
CNEST
C           *IU09*   - INPUT  UNIT COMMON CBOUND (OUTPUT OF PREPROC).
C           *IU10*   - INPUT  UNIT COMMON FBOUND (OUTPUT OF PREPROC).
CNEST
C           *IU11*   - INPUT UNIT OF SPECTRA AT ALL GRID POINTS.
C                      EACH PROPAGATION STEP THE FILES CONNECTED
C                      TO IU11 AND IU12 ARE INTERCHANGED.
C           *IU12*   - OUTPUT UNIT BLOCKS OF SPECTRA (SEE IU11).
CC                     THIS UNIT IS NOT USED ANYMORE IF PBIO SOFTWARE IS
CC                     USED TO READ THE SPECTRA. (FLAG LPBIOIN)
C           *IU13*   - INPUT UNIT OF SPECTRA ON LAST LATITUDE
C                      OF A BLOCK. SPECTRA ARE SAVED FROM THE
C                      SECOND LATITUDE OF NEXT BLOCK.
C                      EACH PROPAGATION STEP THE FILES CONNECTED
C                      TO IU13 AND IU14 ARE INTERCHANGED.
C           *IU14*   - OUTPUT UNIT SECOND LATITUDES (SEE IU13).
C           *IU15*   - OUTPUT UNIT LAST WINDFIELDS.
CC                     THIS UNIT IS NOT USED ANYMORE IF PBIO SOFTWARE IS
CC                     USED TO READ THE LAST WINDFIELDS. (FLAG LPBIOIN)
CREFRA
C           *IU16*   - INPUT/OUTPUT UNIT OF COMMON REFDOT.
CREFRA
C           *IU17*   - INPUT  UNIT OF BLOCKED WINDS.
C                      THIS FILE IS DYNAMICALLY ASSIGNED IN SUB
C                      IMPLSCH. SEE SUB CREWFN TOO.
C           *IU18*   - INPUT  UNIT OF BLOCKED WINDS.
C                      THIS FILE IS DYNAMICALLY ASSIGNED IN SUB
C                      IMPLSCH. SEE SUB CREWFN TOO.
CNEST
C           *IU19*   - OUTPUT UNIT OF BOUNDARY VALUES IF
C                      THIS IS A FINE GRID RUN.
C                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'CBO'
CNEST
C           *IU20*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
C                      OF THE TOTAL SPECTRUM (FIRST GUESS).
C                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'MAP'
C           *IU21*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
C                      OF SWELL AND WIND WAVES (FIRST GUESS).
C                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'SWE'
CASSI
C           *IU22*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
C                      OF THE TOTAL SPECTRUM (ANALYSED).
C                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'AMP'
C           *IU23*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
C                      OF SWELL AND WIND WAVES (ANALYSED)
C                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'ASS'
CASSI
C           *IU25*   - OUTPUT UNIT OF SPECTRA AT CERTAIN GRID POINTS
C                      (FIRST GUESS).
C                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'OUT'
C           *IU26*   - OUTPUT UNIT OF SWELL SPECTRA AT CERTAIN POINTS
C                      (FIRST GUESS).
C                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'SWS'
CASSI
C           *IU27*   - OUTPUT UNIT OF SPECTRA AT CERTAIN GRID POINTS
C                      (ANALYSED).
C                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'AUT'
C           *IU28*   - OUTPUT UNIT OF SWELL SPECTRA AT CERTAIN POINTS
C                      (ANALYSED).
C                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'ASW'
CASSI
C           *IUSCR*  - SCRATCH UNITS FOR WIND PROCESSING.
C           *IUVELO* - OUTPUT UNIT OF BLOCKED WINDFILEDS.
C                      FILES ARE DYNAMICALLY ASSIGNED IN SUB
C                      NOTIM OR TIMIN. SEE SUB CREWFN TOO.
C
C___PACK  THESE UNITS ARE USED FOR FILES CONTAINING GRIB DATA.
C
C         (IN CASE OF SPECTRA ALL SPECTRA OF THE MODEL ARE WRITTEN)
C           *IU30*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
C                      OF THE TOTAL SPECTRUM (FIRST GUESS).
C                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'MPP'
C           *IU31*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
C                      OF SWELL AND WIND WAVES (FIRST GUESS).
C                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'SWP'
CASSI
C           *IU32*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
C                      OF THE TOTAL SPECTRUM (ANALYSED).
C                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'APP'
C           *IU33*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
C                      OF SWELL AND WIND WAVES (ANALYSED)
C                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'AWP'
CASSI
C           *IU35*   - OUTPUT UNIT OF SPECTRA AT ALL GRID POINTS
C                      (FIRST GUESS).
C                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'OUP'
C           *IU36*   - OUTPUT UNIT OF SWELL SPECTRA AT ALL POINTS
C                      (FIRST GUESS).
C                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'WSP'
CASSI
C           *IU37*   - OUTPUT UNIT OF SPECTRA AT ALL GRID POINTS
C                      (ANALYSED).
C                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'AUP'
C           *IU38*   - OUTPUT UNIT OF SWELL SPECTRA AT ALL POINTS
C                      (ANALYSED).
C                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'ASP'
CASSI
C___PACK
C
CLMESSPASS
C           *IU73*   - OUTPUT UNIT OF TABLE CONTAINING MODEL 
C                      DECOMPOSITION PARAMETERS
CLMESSPASS
C
CICE
C           *IU90*   - INPUT UNIT OF SEA SURFACE TEMPERATURE
CICE
C          FOR A START THE RESTART FILES WILL DYNAMICALLY BE ASSIGNED
C          AND COPIED TO THE *OUTPUT* UNITS (IU12, IU14, IU15).
C          IF IT IS REQUESTED TO SAVE RESTART FILES THESE WILL BE COPIED
C          IN REGULAR INTERVALS FROM THE UNIT ALIAS FILES OF IU11 OR
C          IU12, IU13 OR IU14, AND IU15 TO THE PREMANENT RESTART FILES.
C          FOR DETAILS OF THE FILE NAMES SEE SUB GSFILE.
C
C          THE PROGRAM USES OPEN TO ASSIGN FILES WHICH HAVE A FILEID
C          (IU02,IU19,IU20,IU21,IU25,IU26). FOR DETAILS OF THE
C          FILE NAME CONVENTION SEE SUB OPENFIL.
C
C          SUB NOTIM OR TIMIN OPENS FILES AND ASSIGNS THEM TO UNIT
C          IUVELO FOR THE BLOCKED WINDS. THESE FILES ARE READ AND
C          DELETED IN SUB IMPLSCH (IU17 AND IU18). THE FILE NAMES
C          ARE CREATED IN SUB CREWFN.
C
C          THE FILE HANDLING SUBS OPENFIL AND GSFILE ARE COMPUTER
C          DEPENDENT AND MAY BE ADOPTED BY THE USER.
C          IF A FRONT /END SYSTEM IS USED THE CODE BETWEEN 'CCMFA'
C          AND 'CCMFE' HAS TO BE COMMENTED.
C
C          THE PROGRAM CLOSES AND DELETES ALL WORK FILES.
C
C          ALL PARAMETERS HAVE TO BE THE VALUES GIVEN AT THE END
C          OF THE PREPROC OUTPUT IN COLUMN 'REQUIRED'.
C
C     METHOD.
C     -------
C
C          THIS ROUTINE INITIALISES THE WAVEMODEL:
C            -  DEFINES THE UNITS FOR INPUT/OUTPUT,
C            -  READS THE USER INPUT FILE,
C            -  INITIALIZES SOME MODEL CONSTANTS,
C            -  GETS THE RECOVERY FILES,
C            -  READS THE COMMON BLOCKS PRECOMPUTED BY PROG PREPROC,
C            -  DOES SOME GENERAL BOOKEEPING REGARDING
C               DATES, INTEGRATION TIME STEPS AND OUTPUT TIME STEPS.
C            -  READ COMMON UBUF AND SPECTRA IF ONE BLOCK VERSION.
CREFRA
C            -  PRECOMPUTES AND WRITES TO IU16 REFRACTION TERMS.
CREFRA
C            -  OPENS THE FIRST RESULT FILES.
C
C     EXTERNALS.
C     ----------
C
C       *ABORT1*     - TERMINATES PROCESSING.
C       *DIFDATE*   - COMPUTES A TIME DIFFERENCE.
C       *GETREST*   - GETS RESTART FILES.
CREFRA
C       *GRADI*     - COMPUTES DEPTH AND CURRENT GRADIENTS.
CREFRA
C       *GSFILE*    - ROUTINE TO DYNAMICALLY FETCH OR DISPOSE FILES.
C       *GETSPEC*   - READS SPECTRA
C       *GETSTRESS  - READS RESTART STRESS/WIND FILE
C       *HEADBC*    - WRITE BOUNDARY OUTPUT FILE HEADER.
C       *INCDATE*   - UPDATE DATE TIME GROUP.
C       *MPDISTRIBINTFLD* - DISTRIBUTE INTEGER ARRAY ACROSS PE'S
C       *GETSTRESS  - READS RESTART STRESS/WIND FILE
CREFRA
C       *PROPDOT*   - PRECOMPUTE REFRACTION.
CREFRA
CNEST
C       *READBOU*   - READS PREPROC BOUNDARY FILES.
CNEST
CICE
C       *READICE*   - READS SST FROM UNIT 90 AND PREPARES THE ICE EDGE.
CICE
C       *READPRE*   - READS PREPROC OUTPUT.
C       *MPBCASTGRID* - BROACAST GRID FILE READ BY READPRE
C       *USERIN*    - READS USER INPUT.
C
C     REFERENCE
C     ---------
C
C          A MORE DETAILED  DISCUSSION MAY BE FOUND IN SUB WAMODEL.
C
C ----------------------------------------------------------------------
      USE WAVE_MODULE
C
C -------------------------------------------------------------------
C
#include "param.h"
C
#include "txtpara.h"
C
CNEST
#include "comcbou.h"
C
#include "txtcbou.h"
CNEST
C
#include "comcoup.h"
C
#include "txtcoup.h"
C
#include "comcout.h"
C
#include "txtcout.h"
C
CREFRA
#include "comcurr.h"
C
#include "txtcurr.h"
CREFRA
C
CNEST
#include "comfbou.h"
C
#include "txtfbou.h"
CNEST
C
#include "comfred.h"
C
#include "txtfred.h"
C
#include "comgrid.h"
C
#include "txtgrid.h"
C
CICE
#include "comice.h"
C
#include "txtice.h"
CICE
C
#include "comindn.h"
C
#include "txtindn.h"
C
#include "commap.h"
C
#include "txtmap.h"
C
#include "comover.h"
C
#include "txtover.h"
C
CREFRA
#include "comrefd.h"
C
#include "txtrefd.h"
CREFRA
C
#include "comshal.h"
C
#include "txtshal.h"
C
#include "comstat.h"
C
#include "txtstat.h"
C
#include "comtabl.h"
C
#include "txttabl.h"
C
#include "comtemp.h"
C
#include "comtest.h"
C
#include "txttest.h"
C
#include "comtext.h"
C
#include "txttext.h"
C
#include "comunit.h"
C
#include "txtunit.h"
C
#include "comubuf.h"
C
#include "txtubuf.h"
C
#include "comwind.h"
C
#include "txtwind.h"
C
#include "parcons.h"
C
#include "txtcons.h"
C
      NAMELIST /WTR/ PMCOEF, KFLAGWS, KDANG, NWTRMX, KREOSP
C
#include "comwtr.h"
C
#include "txtwtr.h"
C
#include "commesspass.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comwami.h"


      CHARACTER*12 ZERO,CDUM
      CHARACTER*24 FILNM
      CHARACTER*80 LOGFILENAME 
C
C----------------------------------------------------------------------
      REAL X4(2)
      INTEGER I4(2)
      INTEGER,ALLOCATABLE :: KWORK(:,:)
C
C     0. DETERMINE BYTE STORAGE REPRESENTATION OF REAL NUMBERS
C        -----------------------------------------------------
C
      X4=1.
      NPRECR = N_PRECISION(X4)
      I4=1
      NPRECI = N_PRECISION(I4)
C
C     0.1 SET IREAD TO THE PROCESSOR WHICH WILL READ THE INPUT FILES
C         ----------------------------------------------------------
      IREAD=1
C
      IF(NPROC.EQ.1) IREAD=1
C
C*    1. DEFINITION OF MODEL PARAMETERS.
C        -------------------------------
C
 1000 CONTINUE
C 
      ZERO = ' '
C
C*    1.1  DEFINE UNIT NAMES.
C          ------------------
C
CNEST
C 
      IU02 = 2
CNEST
C
C WAMINFO FILE UNIT
C -----------------
      IF(LWCOU) THEN
        IU04 = 27
      ELSE
        IU04 = 4
      ENDIF

C
C STANDARD INPUT UNIT
C --------------------
      IF(LWCOU) THEN
        IU05 = 26
      ELSE
        IU05 = 3
      ENDIF
C
C STANDARD OUTPUT UNIT
C --------------------
      IF(LMESSPASS.AND.LWCOU) THEN
        IU06=20
      ELSE IF (LMESSPASS.AND..NOT.LWCOU) THEN
        IU06=66
        LOGFILENAME='logfile.%p'
        CALL EXPAND_STRING(IRANK,NPROC,0,0,LOGFILENAME,1) 
        OPEN(IU06,FILE=LOGFILENAME,STATUS='UNKNOWN')
        WRITE(*,*) ' STDOUT FILE ', LOGFILENAME, ' OPENED ON PE ', IRANK
        WRITE(IU06,*) ' STDOUT OF PE ', IRANK
        WRITE(IU06,*) ' ============ '
        CALL FLUSH (IU06)
      ELSE
        IU06=6
      ENDIF
C
C GRID AND UBUF FILES UNITS
C -------------------------
      IF(LWCOU) THEN
        IU07 = 29
        IU08 = 31
      ELSE
        IU07 = 7
        IU08 = 8
      ENDIF

CNEST
      IU09 = 9
      IU10 = 10
CNEST
      IU11 = 11

C
C RESTART FILES UNITS
C -------------------
C
C NOTE IF MESSAGE PASSING
C NOTE IU12 WILL ONLY BE CONNECTED TO UNIT 12 FOR THE INPUT
C NOTE OTHERWISE IT WILL CONNECTED TO THE FILENAME BLS WITH
C NOTE AN EXTENSION FUNCTION OF THE PE.
C
      IF(LWCOU) THEN
        IU12 = 36
        IU15 = 37
      ELSE
        IU12 = 12
        IU15 = 15
      ENDIF


      IU13 = 13
      IU14 = 14
CREFRA
      IU16 = 16
CREFRA
      IU17 = 17
      IU18 = 18
CNEST
      IU19 = 19
CNEST
      IU20 = 20
      IU21 = 21
CASSI
      IU22 = 22
      IU23 = 23
CASSI
      IU25 = 25
      IU26 = 26
CASSI
      IU27 = 27
      IU28 = 28
CASSI

C__PACK
      IU30 = 30
      IU31 = 31
CASSI
      IU32 = 32
      IU33 = 33
CASSI
      IU35 = 35
      IU36 = 36
CASSI
      IU37 = 37
      IU38 = 38
CASSI
C__PACK

      IUVELO = 39
      DO 1101 IG=1,NBLO
         IUSCR(IG) = 39+IG
 1101 CONTINUE
C
C TABLE CONTAINING DECOMPOSITION LIMITS UNIT
C -----------------------------------------
CLMESSPASS
C NOTE : THIS UNIT IS ONLY NEEDED IF THE RESTART FILE OPTION IS SET TO
C        IREST = 2 OR 4
      IF(LMESSPASS.AND.LWCOU) THEN
        IU73 = 28
      ELSE IF (LMESSPASS.AND..NOT.LWCOU) THEN
        IU73 = 73
      ENDIF
CLMESSPASS
CICE
C
C SEA SURFACE TEMPERATURE FILE UNIT
C ---------------------------------
      IF(LWCOU) THEN
        IU90 = 38
      ELSE
        IU90 = 90
      ENDIF
CICE

      IF (ITEST.GE.2) THEN
         WRITE(IU06,*) '   SUB. INITMDL: UNITS DEFINED'
         CALL FLUSH (IU06)
      ENDIF
C
C*    1.2 INPUT OF USER PARAMETER.
C         ------------------------
C
      CALL USERIN (IU05, IBOUNC, IBOUNF, IFORCA)
      IF (ITEST.GE.2) THEN
         WRITE(IU06,*) '    SUB. INITMDL: USERIN DONE'
         CALL FLUSH (IU06)
      ENDIF

      IF(LMESSPASS) THEN

       KTAG=100
C
C      1.3 OUTPUT MODEL DECOMPOSITION DETAILS 
C          ----------------------------------
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
          END DO
          CALL FLUSH (IU06)
        ENDIF
C
C
C      1.4 SAVE TABLE CONTAINING DECOMPOSITION LIMITS (IF NEEDED)
C          ------------------------------------------------------
        IF(IREST.EQ.2.OR.IREST.EQ.4) THEN
          IF(IRANK.EQ.1) THEN
            NMAX=NEND(1)-NSTART(1)+1
            DO IP=2,NPROC
              NMAX=MAX(NMAX,NEND(IP)-NSTART(IP)+1)
            ENDDO
            OPEN(IU73,FILE='split_table')
               WRITE(73,*) NPROC,NMAX,NIBLO,NANG,NFRE 
            DO IP=1,NPROC
               WRITE(IU73,*) NSTART(IP)-KLENBOT(IP),
     1                     NEND(IP)+KLENTOP(IP),NSTART(IP),NEND(IP)
            ENDDO
            CLOSE(IU73)
          ENDIF
        ENDIF

C     1.5 SET INTEGRATED OUTPUT PARAMETER TABLE FOR MESSAGE PASSING
C         ----------------------------------------------------------
        CALL MPCRTBL
        IF (ITEST.GE.2) THEN
           WRITE(IU06,*) '    SUB. INITMDL: MPCRTBL DONE'
           CALL FLUSH (IU06)
        ENDIF
      ENDIF
C
C     1.6 DETERMINE LAST OUTPUT DATE
C         --------------------------
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
C

C ----------------------------------------------------------------------
C
C*    2. INPUT FROM PREPROCESSING PROGRAMS.
C        ----------------------------------
C
 2000 CONTINUE
      IF((LMESSPASS.AND.IRANK.EQ.IREAD).OR.(.NOT.LMESSPASS)) THEN
         CALL READPRE (IU07, IREFRA)
         CLOSE (UNIT=IU07, STATUS='KEEP')
      ENDIF
      IF(LMESSPASS) CALL MPBCASTGRID(IU06,IREAD,KTAG,NSTART,NEND,IREFRA)

      IF (ITEST.GE.2) THEN
         IF (IU07 .LT. 111) THEN
           WRITE(IU06,*) '    SUB. INITMDL: PREPROC OUTPUT READ ',
     1                   'AND FILE CLOSE UNIT IU07 = ', IU07
         ELSE
           WRITE(IU06,'(2a,a3)') '    SUB. INITMDL: PREPROC OUTPUT ',
     1                   'READ AND FILE CLOSE UNIT IU07 = ', IU07
         ENDIF
         CALL FLUSH (IU06)
      ENDIF
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
CNEST
!
      IF(LMESSPASS) THEN
         IF(IGL.GT.1) THEN
           WRITE (IU06,*) '******************************************'
           WRITE (IU06,*) '*                                        *'
           WRITE (IU06,*) '*    FATAL ERROR                         *'
           WRITE (IU06,*) '*    ===========                         *'
           WRITE (IU06,*) '* A MULTI BLOCK VERSION RUN WAS REQUESTED*'
           WRITE (IU06,*) '* THIS OPTION IS NOT IMPLEMENTED FOR     *'
           WRITE (IU06,*) '* MESSAGE PASSING PROTOCOL               *'
           WRITE (IU06,*) '*                                        *'
           WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
           WRITE (IU06,*) '*                                        *'
           WRITE (IU06,*) '******************************************'
           CALL ABORT1
         ELSE
! SPECIFY THE LIMITS OF THE SEGMENT OVER WHICH THE PROCESS HAS DIRECT ACCESS
           IJS(1)=NSTART(IRANK)
           IJL(1)=NEND(IRANK)
           IJLT(1)=NSUP
         ENDIF
      ENDIF

C
C*    2.1 READ COMMON CBOUND AND COMMON FBOUND.
C     ----------------------------------------
C
      IF (IBOUNC.EQ.1 .OR. IBOUNF.EQ.1) THEN
         CALL READBOU (IU09, IU10, IU06)
         CLOSE (UNIT=IU09, STATUS='KEEP')
         CLOSE (UNIT=IU10, STATUS='KEEP')
         IF (ITEST.GE.2) THEN
            WRITE(IU06,*)'    SUB. INITMDL: BOUNDARY POINTS READ ',
     1           ' AND FILES CLOSED '
            CALL FLUSH (IU06)
         ENDIF
      ENDIF
CNEST
CICE
C
C*    2.2. READ COMMON ICE
C     --------------------
C
      LICERUN=.FALSE.
C
      FILNM = 'sfctempin'
      LIU   = LEN_TRIM(FILNM)
      FILNM=FILNM(1:LIU)
      INQUIRE(FILE=FILNM,EXIST=LICERUN)

      IF (LICERUN) THEN
        CALL READICE (IU06,FILNM,ITEST,NSTART,NEND,IREAD,KTAG)
      ELSE
        ICEMASK=1
        WRITE(IU06,*)' '
        WRITE(IU06,*)'    ********************************************'
        WRITE(IU06,*)'    SUB. INITMDL: FILE ',FILNM,' NOT FOUND '
        WRITE(IU06,*)'    SUB. INITMDL: NO ICE MASK WILL BE USED '
        WRITE(IU06,*)'    ********************************************'
        WRITE(IU06,*)' '
        CALL FLUSH(IU06)
      ENDIF
C
CICE
C
C ----------------------------------------------------------------------
C
C*    3. PRINT INITIAL CONDITIONS AS READ FROM PERPROCESSING.
C        ----------------------------------------------------
C
 3000 CONTINUE
        WRITE(IU06,*) '  '
        WRITE(IU06,*) ' WAVE MODEL GRID ORGANISATION:'
        WRITE(IU06,3002) ' SOUTHERNMOST LATITUDE IN GRID IS .......: ',
     1                 AMOSOP, ' DEGREE'
        WRITE(IU06,3002) ' NORTHERNMOST LATITUDE IN GRID IS .......: ',
     1                 AMONOP, ' DEGREE'
        WRITE(IU06,3002) ' WESTERNMOST LONGITUDE IN GRID IS .......: ',
     1               AMOWEP, ' DEGREE'
        WRITE(IU06,3002) ' EASTERNMOST LONGITUDE IN GRID IS .......: ',
     1               AMOEAP, ' DEGREE'
        WRITE(IU06,3002) ' LATITUDE INCREMENT IS ..................: ',
     1                 XDELLA, ' DEGREE'
        WRITE(IU06,3002) ' LONGITUDE INCREMENT IS .................: ',
     1               XDELLO, ' DEGREE'
        WRITE(IU06,*) '  '
        WRITE(IU06,3003) ' TOTAL NUMBER OF BLOCKS IS...............: ',
     1                   IGL
        WRITE(IU06,3003) ' TOTAL LENGTH OF EACH BLOCK .............: ',
     1                    NIBLO
        IF(FFLAG25.OR.PFLAG25.OR.GFLAG25.OR.
     &     FFLAG26.OR.PFLAG26.OR.GFLAG26) THEN
          WRITE(IU06,*) '  '
          WRITE(IU06,*) ' OUTPUT POINTS:'
          WRITE(IU06,*) ' TOTAL NUMBER OF OUTPUT POINTS IS: ',NGOUT
          WRITE(IU06,'(/,4X,''BLOCK NO'',6X,''INDEX NO'',
     1                 8X,''LONGITUDE'',6X,''LATITUDE'')')
          DO 3001 I=1,NGOUT
             IX  = IXLG(IJAR(I),IGAR(I))
             KX  = KXLT(IJAR(I),IGAR(I))
             XLO = AMOWEP+REAL(IX-1)*ZDELLO(KX)
             XLA = AMOSOP+REAL(KX-1)*XDELLA

             WRITE(IU06,'((6X,I3,10X,I6,7X,F10.3,4X,F10.3))')
     1            IGAR(I), IJAR(I), XLO, XLA
 3001     CONTINUE
        ENDIF
 3002 FORMAT(3x,a,f8.2,a)
 3003 FORMAT(3x,a,i6,  a)
CNEST
        CALL FLUSH (IU06)
        IF (IBOUNC .EQ. 1) THEN
           WRITE(IU06,*)
           WRITE(IU06,*) ' COARSE GRID: BOUNDARY OUTPUT POINTS :'
           WRITE(IU06,*) ' TOTAL NUMBER OF BOUNDARY POINTS IS: ',NBOUNC
           IF (ITEST.GE.2) THEN
              WRITE(IU06,'(/,4X,''BLOCK NO'',6X,''INDEX NO'',
     1                    8X,''LONGITUDE'',6X,''LATITUDE'')')
              DO 3004 I=1,NBOUNC
                 IX    = IXLG(IJARC(I),IGARC(I))
                 KX    = KXLT(IJARC(I),IGARC(I))
                 XLO = AMOWEP+REAL(IX-1)*ZDELLO(KX)
                 XLA = AMOSOP+REAL(KX-1)*XDELLA

                 WRITE(IU06,'((6X,I3,10X,I5,7X,F10.3,4X,F10.3))')
     1              IGARC(I), IJARC(I), XLO, XLA
 3004         CONTINUE
           ENDIF
        ENDIF
CNEST
C
C ----------------------------------------------------------------------
C
C*    4. CONNECT RESTART FIELDS TO OUTPUT UNITS (IF PBIO SOFTWARE NOT
C        USED). AND READ IN LAST WINDFIELDS FROM RESTARTFILE.
C        ---------------------------------------------------------------
C
 4000 CONTINUE
      IF (ITEST.GE.2)CALL FLUSH(IU06)
C
      CALL GETSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,NSTART,NEND,IREAD)
C
       IF (LWCOU) THEN
         IDELWO = KDELWI
         IDELWI = KDELWI
       ENDIF
C
      CDA = CDTPRO
      IF (CDTPRO.NE.ZERO .OR. LRESTARTED) THEN
C
C*    4.1 MODEL STARTS FROM FILES OUT OF A PREVIOUS MODEL RUN.
C         ----------------------------------------------------
C
         CDATEWO = CDTPRO
         IF (IDELT.LT.IDELWO) CALL INCDATE(CDATEWO,IDELWO/2)
         CDAWIFL = CDTPRO
         IDELWH = MAX(IDELWI,IDELPRO)
         CALL INCDATE(CDAWIFL,IDELWH)
         CDATEFL = CDATEWO

         IF (CDTPRO.LT.CDATEA.OR.CDTPRO.GE.CDATEE.OR.
     1      (IFORCA.EQ.1.AND.CDTPRO.GT.CDATEF) .OR.
     2      (IFORCA.NE.1.AND.CDTPRO.LE.CDATEF)) THEN
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
     1                        CDATEA
            WRITE(IU06,*) ' * START DATE OF FORECAST  IS CDATEF = ',
     1                        CDATEF
            WRITE(IU06,*) ' * END   DATE OF RUN       IS CDATEE = ',
     1                        CDATEE
            WRITE(IU06,*) ' * START DATE FROM RESTART IS CDTPRO = ',
     1                        CDTPRO
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' * PROGRAM ABORTS     PROGRAM ABORTS       *'
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' *******************************************'
            CALL ABORT1
         END IF
      ELSE
C
C*    4.2 MODEL STARTS FROM FIELDS CREATED BY PRESET.
C         -------------------------------------------
C
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
     1                  ' COUNTER INITIALIZED'
         WRITE(IU06,*) '      NEXT WINDFIELD WILL BE READ AT     ',
     1                  'CDATEWO = ',CDATEWO
         WRITE(IU06,*) '      NEXT WIND FILE WILL BE ACCESSED AT ',
     1                  'CDATEFL = ',CDATEFL
         WRITE(IU06,*) '      NEXT WIND FILE NAME IS FROM        ',
     1                  'CDAWIFL = ',CDAWIFL
         CALL FLUSH(IU06)
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    5. INITIALIZE MODEL TIME VARIABLES
C        -------------------------------
C
C*    5.1 OUTPUT TIME VARIABLES.
C         ----------------------
C
      IF (PFLAG20.OR.FFLAG20.OR.GFLAG20) THEN
         CDTINTT= CBPLTDT
 5102    CONTINUE
            CALL INCDATE (CDTINTT, IDELINT)
         IF (CDTINTT.LE.CDTPRO) GOTO 5102
      ELSE
         CDTINTT = ZERO
      ENDIF
      IF (PFLAG21.OR.FFLAG21.OR.GFLAG21) THEN
         CDTINTS= CBPLTDT
 5103    CONTINUE
            CALL INCDATE (CDTINTS, IDELINS)
         IF (CDTINTS.LE.CDTPRO) GOTO 5103
      ELSE
         CDTINTS = ZERO
      ENDIF
      IF (PFLAG25.OR.FFLAG25.OR.GFLAG25) THEN
         CDTSPT= CBPLTDT
 5104    CONTINUE
            CALL INCDATE (CDTSPT, IDELSPT)
         IF (CDTSPT.LE.CDTPRO) GOTO 5104
      ELSE
         CDTSPT = ZERO
      ENDIF
      IF (PFLAG26.OR.FFLAG26.OR.GFLAG26) THEN
         CDTSPS= CBPLTDT
 5105    CONTINUE
            CALL INCDATE (CDTSPS, IDELSPS)
         IF (CDTSPS.LE.CDTPRO) GOTO 5105
      ELSE
         CDTSPS = ZERO
      ENDIF
C
C*    5.2 FILE DISPOSE TIME AND RESTART TIME.
C         -----------------------------------
C
      CDTRES = CBPLTDT
      CALL INCDATE(CDTRES, IDELRES)
      IF (ITEST.GE.2) THEN
         WRITE(IU06,*) '    SUB. INITMDL: TIME COUNTER INITIALIZED'
         CALL FLUSH(IU06)
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    6. CONSISTENCY CHECK ACCORDING TO CFL CRITERION.
C        ---------------------------------------------
C
 6000 CONTINUE
      GVE = G/(ZPI*FR(1)*2.)
      DPH = DELPHI
      CFLP= IDELPRO*GVE/DPH
      DLH = DELLAM(NGY)*COS(AMONOP*RAD)
      CFLL= IDELPRO*GVE/DLH
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
C
C ----------------------------------------------------------------------
C
C*    7. NUMBER OF PROPAGATION TIME STEPS PER CALL.
C        ------------------------------------------
C
 7000 CONTINUE
      NADV = IDELWI/IDELPRO
      NADV = MAX(NADV,1)
      IF (ITEST.GE.2) THEN
         WRITE(IU06,*) '    SUB. INITMDL: NUMBER OF PROPAGATION STEPS'
         WRITE(IU06,*) '                  IN ONE CALL OF SUB WAVEMDL'
         WRITE(IU06,*) '                  WILL BE NADV = ', NADV
         CALL FLUSH(IU06)
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    8. IF ONE BLOCK VERSION READ SPECTRA AND COMMON UBUF.
C        --------------------------------------------------
C
 8000 CONTINUE
      IF (IGL.EQ.1) THEN
C
C*    8.1 READ BLOCK OF SPECTRA.
C         ----------------------
C
         CALL GETSPEC(FL1,NSTART,NEND,IREAD)
C
         IF (ITEST.GE.2) THEN
            WRITE(IU06,*) '    SUB. INITMDL: READ',
     2                    ' BLOCK OF SPECTRA'
            CALL FLUSH (IU06)
          ENDIF
C
C*    8.2.1 INPUT NEIGHBOURING GRID POINT INDICES (COMMON BLOCK UBUF).
C           ----------------------------------------------------------
C
        IF(LMESSPASS) THEN
           ALLOCATE(KWORK(NIBLO,4)) 
           IF(IRANK.EQ.IREAD) READ (IU08) KWORK 
           DO IC=1,4
              CALL MPDISTRIBINTFLD(IREAD,KTAG,NSTART,NEND,KWORK(1,IC))
              KTAG=KTAG+1
           ENDDO
           DO IC=1,2
             DO IJ = NSTART(IRANK),NEND(IRANK)
                KLAT(IJ,IC) = KWORK(IJ,IC) 
                KLON(IJ,IC) = KWORK(IJ,IC+2) 
             ENDDO
           ENDDO
           DEALLOCATE(KWORK)
        ELSE IF(.NOT.LMESSPASS) THEN
           READ (IU08) KLAT, KLON
        ENDIF
C
C     8.2.2 MODIFY KLAT AND KLON SUCH THAT POINT INDICES FOR LAND IS
C           NINF-1.
C           ---------------------------------------------------------
        IF(LMESSPASS) THEN
           DO IC=1,2
             DO IJ = NSTART(IRANK),NEND(IRANK)
                IF(KLAT(IJ,IC).EQ.0) KLAT(IJ,IC) = NINF-1 
                IF(KLON(IJ,IC).EQ.0) KLON(IJ,IC) = NINF-1 
             ENDDO
           ENDDO
        ENDIF
CSHALLOW
C
C*    8.3 COMPUTE SHALLOW WATER TABLE INDICES.
C         ------------------------------------
C
         IF (ISHALLO.NE.1) THEN
            DO 8301 IJ=NINF,IJLT(IGL)
               XD = LOG(DEPTH(IJ,IGL)/DEPTHA)/LOG(DEPTHD)+1.
               ID = NINT(XD)
               ID = MAX(ID,1)
               INDEP(IJ) = MIN(ID,NDEPTH)
 8301       CONTINUE
         ENDIF
CSHALLOW
         IF (ITEST.GE.2)
     1      WRITE(IU06,*) '    SUB. INITMDL: COMMON UBUF ',
     2                 'READ AND INITIALIZED '
      ENDIF
C
C ----------------------------------------------------------------------
CREFRA
C
C*    9. PRECOMPUTE REFRACTION TERMS.
C        ----------------------------
C
 9000 CONTINUE
      IF (IREFRA.NE.0) THEN
         CALL PROPDOT(SL)
         IF (ITEST.GE.2)
     1      WRITE(IU06,*) '    SUB. INITMDL: COMMON REFDOT ',
     2                    'INITIALIZED '
      ENDIF
CREFRA
C ----------------------------------------------------------------------
CCMFA
C
C*    10. OPEN THE FIRST RESULT FILES.
C         ----------------------------
C
C           INTEGRATED PARAMETER DATA FILE.
C
C     IF (FFLAG20)
C    1   CALL GSFILE (IU06, IU20, 0, CDTPRO, CDATEF, 'MAP', 'S')
C
C            INTEGRATED SWELL AND SEA DATA FILE.
C
C     IF (FFLAG21)
C    1   CALL GSFILE (IU06, IU21, 0, CDTPRO, CDATEF, 'SWE', 'S')
C
C             SPECTRA FILE.
C
C     IF (FFLAG25)
C    1   CALL GSFILE (IU06, IU25, 0, CDTPRO, CDATEF, 'OUT', 'S')
C
C             SWELL SPECTRA FILE.
C
C     IF (FFLAG26)
C    1   CALL GSFILE (IU06, IU26, 0, CDTPRO, CDATEF, 'SWS', 'S')
CNEST
C             BOUNDARY VALUE FILE.
C
C     IF (IBOUNC.EQ.1)
C    1   CALL GSFILE (IU06, IU19, 0, CDTPRO, CDATEF, 'CBO', 'S')
CNEST

C     IF (ITEST.GE.2)
C    1   WRITE(IU06,*) '   SUB. INITMDL: OUTPUT FILES OPENED'
CCMFE
CNEST
C
C        WRITE BOUNDARY VALUE FILE HEADER.
C
      IF (IBOUNC.EQ.1) THEN
            CALL HEADBC (NBOUNC, IDELPRO, TH(1), FR(1), IU19, IU06)
         IF (ITEST.GE.2)
     1      WRITE(IU06,'(''    SUB. INITMDL: HEADER FOR  '',
     2      ''COARSE GRID WAS WRITTEN OF UNIT IU19 = '',A8)') IU19
      ENDIF
CNEST

      RETURN
      END
