      SUBROUTINE USERIN (IU05, IBOUNC, IBOUNF, IFORCA)

C ----------------------------------------------------------------------
C
C**** *USERIN* - ROUTINE TO READ AND WRITE NAMELIST INPUT.
C
C     H. GUNTHER     GKSS/ECMWF     NOVEMBER 1989
C     J. BIDLOT      ECMWF          JUNE  1996
C     J. BIDLOT      ECMWF          JUNE  1997: coupled/uncoupled
C                                               restarted/normal
C     B. HANSEN      ECMWF          JULY  1997  COUPLED AND ANALYSIS
C                                               (IE. PREPAN)
C     B. HANSEN      ECMWF          JANUARY 1998 NAMELIST INPUT.
C
C*    PURPOSE.
C     --------
C
C       READ USER INPUT CONCERNING PERIOD OF INTEREST,TIMESTEPS AND
C       MODEL OPTIONS TO INITIALIZE COMMON USERD. A CONSISTENCY CHECK
C       IS DONE TOO.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *USERIN (IU05, IBOUNC, IBOUNF, IFORCA)*
C          *IU05*   -  INPUT UNIT USER INPUT FILE.
C          *IBOUNC* -  COARSE GRID BOUNDARY FLAG.
C          *IBOUNF* -  FINE GRID BOUNDARY FLAG.
C          *IFORCA* -  FORCAST START OPTION.
C
C     METHOD.
C     -------
C
C        USER INFORMATION IS BEING READ WITH THE PRESUMPTIONS THAT:
C         1. EVERY LINE STARTING WITH 'C' IS A COMMENT LINE
C         2. VALUES ARE PUT IN BELOW POSITIONS INDICATED WITH '-'
C            (RIGHT-JUSTIFIED)
C
C     EXTERNALS.
C     ----------
C
C       *ABORT1*    - TERMINATES PROCESSING.
C       *DIFDATE*   - COMPUTES A TIME DIFFERENCE.
C       *IINITFDB*  - INITIALISE DATA BASE.

C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comcoup.h"
C
#include "comcout.h"
!
#include "commesspass.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comstat.h"
C
#include "comtest.h"
C
#include "comtext.h"
C
#include "comunit.h"
C
#include "comwami.h"
C
C ----------------------------------------------------------------------
C
      LOGICAL LERROR
      CHARACTER (LEN=  10), PARAMETER :: ZERO = ' '
      CHARACTER (LEN=  70)            :: CLHEADER
      CHARACTER (LEN= 256)            :: CLFORM
C
C ----------------------------------------------------------------------
C
      CHARACTER*1 CLMTSU(4), CLOTSU(7)
      NAMELIST /NALINE/ CLHEADER,
     1   CBPLTDT, CEPLTDT, CDATEF,
     2   IDELPRO, IDELT, IDELWO, IDELWI, CLMTSU, IDELALT,
     3   IDELINT, IDELINS, IDELAGI, IDELSPT, IDELSPS, IDELAGS, IDELRES,
     4   CLOTSU, CDATER, CDATES,
     5   PFLAG,  FFLAG,  GFLAG,
     6   LFDB,  LPBIOIN, LGRIBIN, LPBIOOUT, LFDBIOIN, LFDBIOOUT,
     7   ICASE, ISHALLO, IREFRA, ITEST, ITESTB, IREST, IASSI, IFORCA,
     8   IBOUNC, IBOUNF,
     9   USERID, RUNID,  PATH,   YCLASS, YEXPVER, CPATH,
     x   NENSFNB, NTOTENS,
     1   LMESSPASS, LWCOU, NTASKS, NSIZE, CFDBSF, CFDB2DSP, LNOCDIN
C
      CHARACTER* 10 CLOUT
      NAMELIST /NAOT/ CLOUT
C
C ----------------------------------------------------------------------
C
C*    0. SET DEFAULT VALUES FOR THE NAMELIST ELEMENTS.
C        ---------------------------------------------
C
      CLMTSU    = 'S'
      CLOTSU    = 'S'
      CLHEADER  = ZERO
      CBPLTDT   =  ZERO
      CEPLTDT   =  ZERO
      CDATEF    =  ZERO
      IDELPRO   =  0
      IDELT     =  0
      IDELWO    =  0
      IDELWI    =  0
      IDELALT   =  0
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
      GFLAG     = .FALSE. 
      LFDB      = .FALSE. 
      LPBIOIN   = .FALSE. 
      LGRIBIN   = .FALSE. 
      LPBIOOUT  = .FALSE. 
      LFDBIOIN  = .FALSE. 
      LFDBIOOUT = .FALSE. 
      ICASE     = 0 
      ISHALLO   = 0 
      IREFRA    = 0 
      ITEST     = 0 
      ITESTB    = 0 
      IREST     = 0 
      IASSI     = 0 
      IFORCA    = 0 
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
      LMESSPASS = .TRUE.
      NTASKS    = 0
      NSIZE     = 0
C
      COUTT     = ZERO
      NOUTT     = 0
C
      CFDBSF    = "/vpp700/fdb18"
      CFDB2DSP  = ""
      LNOCDIN   = .FALSE. 
C
C ----------------------------------------------------------------------
C
C*    1. READ NAMELIST NALINE.
C        ---------------------
C
      REWIND(IU05)
C
      READ (IU05, NALINE)
      IF(.NOT.LWCOU) THEN
        CDATEA = CBPLTDT
        CDATEE = CEPLTDT
      ELSE
        CBPLTDT = CDATEA
        CEPLTDT = CDATEF
      ENDIF
      CALL u2l1cr( YCLASS )
      CFDB2DSP  = ""
C
C           **** MODEL TIME STEPS ****
      IF (CLMTSU(1) .EQ. 'H') IDELPRO = IDELPRO*3600
      IF (CLMTSU(2) .EQ. 'H') IDELT   = IDELT*3600
      IF (CLMTSU(3) .EQ. 'H') IDELWO  = IDELWO*3600
      IF (CLMTSU(4) .EQ. 'H') IDELWI  = IDELWI*3600
C           **** OUTPUT TIME IN FIXED INTERVALS ****
      IF (CLOTSU(1) .EQ. 'H') IDELINT = IDELINT*3600
      IF (CLOTSU(2) .EQ. 'H') IDELINS = IDELINS*3600
      IF (CLOTSU(3) .EQ. 'H') IDELAGI = IDELAGI*3600
      IF (CLOTSU(4) .EQ. 'H') IDELSPT = IDELSPT*3600
      IF (CLOTSU(5) .EQ. 'H') IDELSPS = IDELSPS*3600
      IF (CLOTSU(6) .EQ. 'H') IDELAGS = IDELAGS*3600
      IF (CLOTSU(7) .EQ. 'H') IDELRES = IDELRES*3600
C           **** CHECK LENGTH OF YEXPVER AND PUT IT RIGHT JUSTIFIED ****
      LYEP=IECF_LEN(YEXPVER)
      YEXPVER(5-LYEP:4)=YEXPVER(1:LYEP)
      YEXPVER(1:4-LYEP)='0000'
C           **** OUTPUT TIME AT SPECIFIED TIMES ****
      REWIND(IU05)
      SPT: DO
        READ (IU05, NAOT, END=1900)
        NOUTT = NOUTT+1
        IF (NOUTT.GT.MOUTT) EXIT SPT
        COUTT(NOUTT) = CLOUT
      ENDDO SPT

      IF(LWCOU) THEN
C     IF COUPLED MODEL
C     SET TEMPORARILY IDELWO and IDELWI TO IDELPRO
C     to AVOID ALL THE CHECKS. SEE IN INITMDL WHERE THEY ARE RESET
        IDELWO=IDELPRO
        IDELWI=IDELPRO
      ENDIF
C
C   1 FORMAT(A,A3)
C   2 FORMAT(A,A3,A,I10,A)
C
C* CHECK FLAG FOR GRIBING AS SOME OPTION ARE NOT IMPLEMENTED
C* IN ROUTINE OUTINT
C
 1900 IF (GFLAG(4)) THEN
          GFLAG(4) = .FALSE.
          WRITE(IU06,*) ' ******************* NOTE ********************'
          WRITE(IU06,*) ' NO GRIB OUTPUT POSSIBLE FOR FRICTION VELOCITY'
          WRITE(IU06,*) ' NOT IMPLEMENTED IN ROUTINE OUTINT ' 
          WRITE(IU06,*) ' '
      ENDIF
      IF (GFLAG(5)) THEN
          GFLAG(5) = .FALSE.
          WRITE(IU06,*) ' ******************* NOTE ***************** ' 
          WRITE(IU06,*) ' NO GRIB OUTPUT POSSIBLE FOR WIND DIRECTION'
          WRITE(IU06,*) ' NOT IMPLEMENTED IN ROUTINE OUTINT ' 
          WRITE(IU06,*) ' '
      ENDIF
      IF (GFLAG(8)) THEN
          GFLAG(8) = .FALSE.
          WRITE(IU06,*) ' *********** NOTE *************** ' 
          WRITE(IU06,*) ' NO GRIB OUTPUT POSSIBLE FOR' 
          WRITE(IU06,*) ' NORMALISED WAVE STRESS'
          WRITE(IU06,*) ' NOT IMPLEMENTED IN ROUTINE OUTINT ' 
          WRITE(IU06,*) ' '
      ENDIF
      IF (IREST.GT.0.AND.LGRIBIN.AND..NOT.GFLAG(7)) THEN
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
      IF (IREST.GT.0.AND.LGRIBIN.AND..NOT.GFLAG(10)) THEN
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
      FFLAG20 = FFLAG( 1).OR.FFLAG( 2).OR.FFLAG( 3).OR.FFLAG( 4).OR.
     1          FFLAG( 5).OR.FFLAG( 6).OR.FFLAG( 7).OR.FFLAG( 8).OR.
     2          FFLAG( 9).OR.FFLAG(10)
      PFLAG20 = PFLAG( 1).OR.PFLAG( 2).OR.PFLAG( 3).OR.PFLAG( 4).OR.
     1          PFLAG( 5).OR.PFLAG( 6).OR.PFLAG( 7).OR.PFLAG( 8).OR.
     2          PFLAG( 9).OR.PFLAG(10)
      GFLAG20 = GFLAG( 1).OR.GFLAG( 2).OR.GFLAG( 3).OR.GFLAG( 4).OR.
     1          GFLAG( 5).OR.GFLAG( 6).OR.GFLAG( 7).OR.GFLAG( 8).OR.
     2          GFLAG( 9).OR.GFLAG(10)
      FFLAG21 = FFLAG(11).OR.FFLAG(12).OR.FFLAG(13).OR.FFLAG(14).OR.
     1          FFLAG(15).OR.FFLAG(16) 
      PFLAG21 = PFLAG(11).OR.PFLAG(12).OR.PFLAG(13).OR.PFLAG(14).OR.
     1          PFLAG(15).OR.PFLAG(16) 
      GFLAG21 = GFLAG(11).OR.GFLAG(12).OR.GFLAG(13).OR.GFLAG(14).OR.
     1          GFLAG(15).OR.GFLAG(16) 
      FFLAG25 = FFLAG(17)
      PFLAG25 = PFLAG(17)
      GFLAG25 = GFLAG(17)
      FFLAG26 = FFLAG(18)
      PFLAG26 = PFLAG(18)
      GFLAG26 = GFLAG(18)
C
      ICPLEN=IECF_LEN(CPATH)
      CWI=CPATH(1:ICPLEN)//'/waminfo'
C
C*    1.1  READ THE WAMINFO FILE AND OVERWRITE INPUT.
C          ------------------------------------------
C
      INQUIRE(FILE=CWI,EXIST=LRESTARTED)
      IF (LRESTARTED) THEN
        OPEN  (IU04,FILE=CWI,STATUS='OLD', FORM='FORMATTED')
        CALL READSTA(IU04, CBEGDT, CENDDT, IANALPD, IFOREPD, IDELWIN,
     1               CRSTDT, CLSPDT, CBPLTDT, CEPLTDT, IASSIM, NFCST,
     2               ISTAT)
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
C
        CDATEF = CEPLTDT
        IF (.NOT. LWCOU) THEN
          CDATEA = CBEGDT
C
          CDATEE = CDATEF
C         IF ( IANALPD.LT.0 ) IFOREPD=IFOREPD+IANALPD*(-1) 
          CALL INCDATE (CDATEE,IFOREPD)
        ENDIF
C
        CDATER = CRSTDT
        CDATES = CLSPDT
        IDELWI = IDELWIN
        IASSI = IASSIM
        IFORCA = NFCST
C
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
C
        IF(LFDBIOIN) THEN
          WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
          WRITE(IU06,*)'+                                          +'
          WRITE(IU06,*)'+ LFDBIOIN HAS BEEN RESET TO FALSE         +'
          WRITE(IU06,*)'+ IN ORDER TO READ RESTART FILES FROM DISK +'
          WRITE(IU06,*)'+                                          +'
          WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
          LFDBIOIN  = .FALSE. 
        ENDIF
C
      ENDIF
C

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

      IF(CDATER .NE. '0000000000') THEN
        IF(CDATER .LT. CDATEA) CDATER = CDATEE 
        IF(CDATER .GT. CDATEE) CDATER = CDATEE 
      ENDIF
      IF(CDATES .LT. CDATEA) CDATES = CDATEE 
      IF(CDATES .GT. CDATEE) CDATES = CDATEE 
C
C*    1.2 CHECK NUMBER OF OUTPUT TIMES.
C         -----------------------------
C
      IF (NOUTT.GT.MOUTT) THEN
         WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++'
         WRITE(IU06,*) '+                                           +'
         WRITE(IU06,*) '+    WARNING ERROR IN SUB. USERIN           +'
         WRITE(IU06,*) '+    ============================           +'
         WRITE(IU06,*) '+ NUMBER OF OUTPUT TIMES IN INPUT EXCEEDS   +'
         WRITE(IU06,*) '+ DIMENSION MOUTT                = ', MOUTT
         WRITE(IU06,*) '+ NUMBER OF TIMES INPUT IS NOUTT = ', NOUTT
         WRITE(IU06,*) '+ PROGRAM WILL IGNORE THE LAST OUTPUT TIMES +'
         WRITE(IU06,*) '+                                           +'
         WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++'
        NOUTT = MOUTT
      END IF
C
C     1.3 CHECK IF IDELALT WAS SET ELSE SET IT TO IDELWO
C         ----------------------------------------------
C
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
C
C*    1.4  PRINT INITIAL CONDITIONS.
C          -------------------------
        WRITE(IU06,*) '  '
        IV = INDEX(CLHEADER,'4.0')
        IF (IV .GT. 30 ) THEN
          CLHEADER(IV:IV+2) = '4.7'
        ENDIF
        WRITE(IU06,*) ' ',CLHEADER
        WRITE(IU06,*) ' COUPLING WITH ATMOS. MODEL (LWCOU) : ',LWCOU
        IF(.NOT.LMESSPASS) THEN
        WRITE(IU06,*) ' USING ',NTASKS,' PROCESSORS'
        WRITE(IU06,*) ' WITH A WORKLOAD OF ',NSIZE,' PER PROCESSOR'
        ENDIF
        WRITE(IU06,*) '  '
        WRITE(IU06,*) ' STARTING  DATE (FORMAT:YYMMDDHHMM) : ',CDATEA,
     1                ' END DATE :',CDATEE
        WRITE(IU06,*) '  '
        WRITE(IU06,*) ' MODEL TIME STEPS:'
        WRITE(IU06,*) ' SOURCE TERM INTEGRATION TIME STEP : ',
     1                  IDELT,' SECS'
        WRITE(IU06,*) ' PROPAGATION TIME STEP ............: ',
     1                  IDELPRO,' SECS'
        IF(.NOT.LWCOU) THEN
          WRITE(IU06,*) ' MODEL WIND INPUT TIME STEP .......: ',
     1                    IDELWI,' SECS'
          WRITE(IU06,*) ' MODEL WIND OUTPUT TIME STEP.......: ',
     1                  IDELWO,' SECS'
        ENDIF
        WRITE(IU06,*) '  '
        WRITE(IU06,*) ' MODEL OPTIONS:'
        WRITE(IU06,*) '  '
        IF (ICASE.EQ.1) THEN
         WRITE(IU06,*) ' PROPAGATION GRID SPHERICAL LAT/LON COORDINATES'
        ELSE
          WRITE(IU06,*) ' PROPAGATION GRID CARTESIAN COORDINATES'
        ENDIF
        IF (ISHALLO.EQ.1) THEN
           WRITE(IU06,*) ' THIS IS A DEEP WATER RUN '
        ELSE
           WRITE(IU06,*) ' THIS IS A SHALLOW WATER RUN '
        ENDIF
        IF (IREFRA.EQ.0) THEN
           WRITE(IU06,*) ' MODEL RUNS WITHOUT REFRACTION'
        ELSE IF (IREFRA.EQ.1) THEN
           WRITE(IU06,*) ' MODEL RUNS WITH DEPTH REFRACTION'
        ELSE IF (IREFRA.EQ.2) THEN
         WRITE(IU06,*) ' MODEL RUNS WITH DEPTH AND CURRENT REFRACTION'
        ENDIF
        IF (IDELWO.GE.IDELWI) THEN
           WRITE(IU06,*) ' WIND FIELDS ARE NOT INTERPOLATED IN TIME'
        ELSE
           WRITE(IU06,*) ' WIND FIELDS ARE INTERPOLATED IN TIME'
        ENDIF
        IF (IBOUNC.EQ.1) THEN
         WRITE(IU06,*) ' MODEL RUNS WITH BOUNDARY POINTS (COURSE GRID)'
        ELSE
           WRITE(IU06,*) ' MODEL RUNS WITHOUT BOUNDARY POINTS',
     1                   ' (COURSE GRID)'
        ENDIF
        IF (IBOUNF.EQ.1) THEN
           WRITE(IU06,*) ' MODEL RUNS WITH BOUNDARY POINTS (FINE GRID)'
        ELSE
           WRITE(IU06,*) ' MODEL RUNS WITHOUT BOUNDARY POINTS',
     1                   ' (FINE GRID)'
        ENDIF
        IF (IFORCA.EQ.1) THEN
           WRITE(IU06,*) ' MODEL STARTS FROM ANALYSIS FIELDS'
        ELSE
           WRITE(IU06,*) ' MODEL STARTS FROM FORECAST FIELDS'
        ENDIF
        IF (IASSI.EQ.1) THEN
           WRITE(IU06,*) ' WAVE DATA ASSIMILATION IS CARRIED OUT'
           IF(LWCOU) THEN
           WRITE(IU06,*) ' AT DATE ', CDATEF
           ELSE
           WRITE(IU06,*) ' UNTIL THE END OF THE ANALYSIS PERIOD'
           WRITE(IU06,*) ' AT DATE ', CDATEF
           ENDIF
           WRITE(IU06,*) ' WITH ALTIMETER DATA IN TIME WINDOW(S) OF ' 
           WRITE(IU06,*) ' IDELALT = ', IDELALT,' SECONDS' 
           WRITE(IU06,*) ' CENTERED AROUND THE ASSIMILATION TIME(S) '

        ELSE
           WRITE(IU06,*) ' DATA ASSIMILATION IS NOT CARRIED OUT'
        ENDIF
        WRITE(IU06,*) '  '
        WRITE(IU06,*) ' MODEL OUTPUT SELECTION:'
        WRITE(IU06,*) '  '
        IF (NOUTT.GT.0) THEN
           WRITE(IU06,*) ' NUMBER OF OUTPUT TIMES IS NOUTT = ', NOUTT
           WRITE(IU06,*) ' OUTPUT WILL BE PROCESSED AT:'
           WRITE(IU06,'(6(2X,A10))') (COUTT(I),I=1,NOUTT)
           WRITE(IU06,*) '  '
        ENDIF

         IF(FFLAG20) CALL OUTXT(IU06,IU20,IU30,IU20,IU32,IDELINT,IASSI,
     1                 NOUTT,'IU20','IU30','IU20','IU32')

        WRITE(IU06,*) '    F = FALSE   T = TRUE '
        WRITE(IU06,*) '                        PRINTER   UNIT    GRIB'
        WRITE(IU06,*) ' SIGNIFICANT WAVE HEIGHT ..', PFLAG( 1),
     1                '......', FFLAG( 1),'.....',GFLAG( 1)
        WRITE(IU06,*) ' MEAN WAVE DIRECTION ......', PFLAG( 2),
     1                '......', FFLAG( 2),'.....',GFLAG( 2)
        WRITE(IU06,*) ' WAVE MEAN FREQUENCY ......', PFLAG( 3),
     1                '......', FFLAG( 3),'.....',GFLAG( 3)
        WRITE(IU06,*) ' FRICTION VELOCITY ........', PFLAG( 4),
     1                '......', FFLAG( 4),'.....',GFLAG( 4)
        WRITE(IU06,*) ' WIND DIRECTION ...........', PFLAG( 5),
     1                '......', FFLAG( 5),'.....',GFLAG( 5)
        WRITE(IU06,*) ' WAVE PEAK FREQUENCY ......', PFLAG( 6),
     1                '......', FFLAG( 6),'.....',GFLAG( 6)
        WRITE(IU06,*) ' DRAG COEFFICIENT .........', PFLAG( 7),
     1                '......', FFLAG( 7),'.....',GFLAG( 7)
        WRITE(IU06,*) ' NORMALISED WAVE STRESS ...', PFLAG( 8),
     1                '......', FFLAG( 8),'.....',GFLAG( 8)
        WRITE(IU06,*) ' MEAN SQUARE SLOPE.........', PFLAG( 9),
     1                '......', FFLAG( 9),'.....',GFLAG( 9)
        WRITE(IU06,*) ' WAVE MODEL WIND SPEED.....', PFLAG(10),
     1                '......', FFLAG(10),'.....',GFLAG(10)
        WRITE(IU06,*) '  '

        IF(FFLAG21) CALL OUTXT(IU06,IU21,IU31,IU23,IU33,IDELINS,IASSI,
     1                 NOUTT,'IU21','IU31','IU23','IU33')
        WRITE(IU06,*) '                        PRINTER   UNIT    GRIB'
          WRITE(IU06,*) ' WIND SEA WAVE HEIGHT......', PFLAG(11),
     1                '......', FFLAG(11),'.....',GFLAG(11)
        WRITE(IU06,*) ' WIND SEA DIRECTION........', PFLAG(13),
     1                '......', FFLAG(12),'.....',GFLAG(13)
        WRITE(IU06,*) ' WIND SEA MEAN FREQUENCY...', PFLAG(15),
     1                '......', FFLAG(15),'.....',GFLAG(15)
        WRITE(IU06,*) ' SWELL WAVE HEIGHT.........', PFLAG(12),
     1                '......', FFLAG(12),'.....',GFLAG(12)
        WRITE(IU06,*) ' SWELL WAVE DIRECTION......', PFLAG(14),
     1                '......', FFLAG(14),'.....',GFLAG(14)
        WRITE(IU06,*) ' SWELL MEAN FREQUENCY .....', PFLAG(16),
     1                '......', FFLAG(16),'.....',GFLAG(16)
        WRITE(IU06,*) ' PACKING UNDER TESTING !!!'
        WRITE(IU06,*) '  '

C
        IF(FFLAG25) CALL OUTXT(IU06,IU25,IU35,IU27,IU37,IDELSPT,IASSI,
     1                 NOUTT,'IU25','IU35','IU27','IU37')
        WRITE(IU06,*) '                        PRINTER   UNIT    GRIB'
        WRITE(IU06,*) ' SPECTRA OF TOTAL SEA .....', PFLAG(17),
     1              '......', FFLAG(17),'.....',GFLAG(17)
        WRITE(IU06,*) ' GRIBBING NOT YET IMPLEMENTED !!!'
        WRITE(IU06,*) '  '
C
        IF(FFLAG26) CALL OUTXT(IU06,IU26,IU36,IU28,IU38,IDELSPS,IASSI,
     1                 NOUTT,'IU26','IU36','IU28','IU38')
        WRITE(IU06,*) '                        PRINTER   UNIT    GRIB'
        WRITE(IU06,*) ' SWELL SPECTRA ............', PFLAG(18),
     1                '......', FFLAG(18),'.....',GFLAG(18)
        WRITE(IU06,*) ' GRIBBING NOT YET IMPLEMENTED !!!'
        WRITE(IU06,*) '  '

        WRITE(IU06,*) '  '

        WRITE(IU06,*) ' ACCESS TO THE FIELD DATA BASE: '
     1             ,'    F = DISABLED   T = ENABLED    ', LFDB
        IF (LFDB) THEN
          ILCFDBSF = LEN_TRIM(CFDBSF)
          CLFORM='("  IT''s ROOT DIRECTORY IS ", A256)'
          WRITE (CLFORM(32:34), '(I3.3)') ILCFDBSF
          WRITE(IU06,CLFORM) CFDBSF(1:ILCFDBSF)
        ENDIF
        WRITE(IU06,*) '  '


      IF ( LFDB .AND. .NOT. (GFLAG20.OR.GFLAG21))
     1  THEN
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
            WRITE(IU06,*) ' OUTPUT OF GRIBBED PARAMETERS REDIRECTED '
     .                   ,' TO FIELD DATA BASE'
            WRITE(IU06,*) '                    '
      ENDIF

        WRITE(IU06,*) '  '
        IF (IREST.EQ.1.OR.IREST.EQ.3) THEN
        WRITE(IU06,*) ' OUTPUT AND RESTART FILES WILL BE DISPOSED EVERY'
     1               ,'...', IDELRES, ' SECONDS'
          IF(CDATER.LT.CDATEE) WRITE(IU06,*)
     1              ' !! HOWEVER BOTH RESTART FILES WILL ONLY BE SAVED',
     2              ' AT ...', CDATER 
          IF(CDATES.LT.CDATEE) WRITE(IU06,*)
     1                ' BUT SPECTRA FILES ALONE  WILL BE SAVED UNTIL '
     2               ,'...', CDATES 
          IF(IREST.EQ.3) THEN
             WRITE(IU06,*) ' RESTART FILES WILL BE PREALLOCATED '
             WRITE(IU06,*) ' MAKE SURE WE ARE USING A VFL FILE SYSTEM'
          ENDIF
        ELSEIF (IREST.EQ.2.OR.IREST.EQ.4) THEN
        WRITE(IU06,*) ' OUTPUT AND RESTART FILES WILL BE DISPOSED EVERY'
     1               ,'...', IDELRES, ' SECONDS'
        WRITE(IU06,*) ' !!!! RESTART FILES ARE SAVED IN SPLIT MODE'
          IF(CDATER.LT.CDATEE) WRITE(IU06,*)
     1              '!!! HOWEVER BOTH RESTART FILES WILL ONLY BE SAVED',
     2              ' AT ...', CDATER 
          IF(CDATES.LT.CDATEE) WRITE(IU06,*)
     1                ' BUT SPECTRA FILES ALONE  WILL BE SAVED UNTIL '
     2               ,'...', CDATES 
          IF(IREST.EQ.4) THEN
             WRITE(IU06,*) ' RESTART FILES WILL BE PREALLOCATED '
             WRITE(IU06,*) ' MAKE SURE WE ARE USING A VFL FILE SYSTEM'
          ENDIF
        ELSE
        WRITE(IU06,*) ' OUTPUT FILES WILL BE DISPOSED EVERY ...........'
     1               ,'...', IDELRES, ' SECONDS'
        ENDIF
        WRITE(IU06,*) '  '
        IF(LPBIOIN) THEN
           WRITE (IU06,*)
     1     ' PBIO SOFTWARE IS USED TO READ SELECTED INPUT FILES' 
        ENDIF
        IF(LNOCDIN) THEN
           WRITE(IU06,*) '  '
           WRITE (IU06,*) 
     1     ' NO DRAG COEFFICIENT FIELD IS PROVIDED AS INPUT'
           WRITE (IU06,*) ' THE FIELD WILL BE INITIALISED BY TAKING' 
           WRITE (IU06,*) ' ZERO WAVE STRESS (TAUW)' 
        ENDIF
        IF(LGRIBIN) THEN
           WRITE(IU06,*) '  '
           WRITE (IU06,*)
     1     ' SPLITTED GRIB SPECTRA FIELD ARE USED AS INPUT' 
        ENDIF
        WRITE(IU06,*) '  '
        IF(LPBIOOUT) THEN
           WRITE (IU06,*)
     1     ' PBIO SOFTWARE IS USED TO WRITE TO SELECTED OUTPUT FILES' 
        ENDIF
        WRITE(IU06,*) '  '
        IF(LFDBIOIN) THEN
           WRITE (IU06,*)
     1     ' FDB SOFTWARE IS USED TO READ SELECTED INPUT FILES' 
        ENDIF
        IF(LFDBIOOUT) THEN
           WRITE (IU06,*)
     1     ' FDB SOFTWARE IS USED TO WRITE TO SELECTED OUTPUT FILES' 
        ENDIF
        IF (LFDBIOIN.OR.LFDBIOOUT) THEN
          CLFORM='("  IT''s ROOT DIRECTORY IS ", A256)'
          ILCFDB2DSP = LEN_TRIM(CFDB2DSP)
          IF ( ILCFDB2DSP .EQ. 0 ) THEN
            ILCFDBSF = LEN_TRIM(CFDBSF)
            WRITE (CLFORM(32:34), '(I3.3)') ILCFDBSF
            WRITE(IU06,CLFORM) CFDBSF(1:ILCFDBSF)
          ELSE
            WRITE (CLFORM(32:34), '(I3.3)') ILCFDB2DSP
            WRITE(IU06,CLFORM) CFDB2DSP(1:ILCFDB2DSP)
          ENDIF
        ENDIF
        WRITE(IU06,*) '  '

        WRITE(IU06,'("  HARD DRIVE PATH NAME : ",/,5X,A70)') CPATH
        WRITE(IU06,'("  CURRENT RUN:")')
        WRITE(IU06,'("  CLASS......: ", A4)') YCLASS
        WRITE(IU06,'("  EXPERIMENT.: ", A4)') YEXPVER
        IF (NENSFNB.NE.0.OR.NTOTENS.NE.0) THEN
          WRITE(IU06,'("  ENSEMBLE FORECAST RUN : ")')
          WRITE(IU06,'("  *********************** ")')
          WRITE(IU06,'("  ENSEMBLE NUMBER:         ", I4  )') NENSFNB 
          WRITE(IU06,'("  TOTAL NUMBER OF ENSEMBLE:", I4,/)') NTOTENS 
        ENDIF
        WRITE(IU06,*) '  '
        WRITE(IU06,*) ' TEST OUTPUT LEVEL IS .............. ITEST = ',
     1                  ITEST
        WRITE(IU06,*) ' TEST OUTPUT OF IN BLOCK LOOPS UPTO ITESTB = ',
     1                  ITESTB
        WRITE(IU06,*) '  '
C
        IF (.NOT. LWCOU .AND.
     .    (  (LFDB .OR. LGRIBIN .OR. LFDBIOIN .OR. LFDBIOOUT ))) THEN
          WRITE(IU06,*) ' ACCESS TO FIELD DATA BASE REQUESTED '
          WRITE(IU06,*) ' DATA BASE WILL BE INITIALIZED '
          INITFDB_STAT = IINITFDB()
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
     .       INITFDB_STAT
             WRITE(IU06,'("*                                      *")')
             WRITE(IU06,'("****************************************")')
             CALL ABORT1
          ENDIF
        ENDIF
C
C ----------------------------------------------------------------------
C
C*    2. CHECK INTEGER RATIOS BETWEEN TIMESTEPS.
C        ---------------------------------------
C
      LERROR = .FALSE.
C
C*    2.1 WIND OUTPUT AND PROPAGATION TIME STEP.
C         --------------------------------------
C
      IF ((IDELWO.LT.IDELPRO.AND.MOD(IDELPRO,IDELWO).NE.0) .OR.
     1    (IDELWO.GE.IDELPRO.AND.MOD(IDELWO,IDELPRO).NE.0)) THEN
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
      END IF
C
C*    2.2 SOURCE FUNCTION AND PROPAGATION TIMESTEP.
C         -----------------------------------------
C
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
      END IF
C
C*    2.3 SOURCE FUNCTION AND WIND OUTPUT TIMESTEP.
C         -----------------------------------------
C
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
      END IF
C
C*    2.4 WIND INPUT AND WIND OUTPUT TIMESTEP.
C         ------------------------------------
C
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
      END IF
C
C*    2.5 FILE DISPOSE TIMESTEP.
C         ----------------------
C
      IF (MOD(IDELRES,IDELPRO).NE.0 .AND.
     1    MOD(IDELRES,IDELWI).NE.0) THEN
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
      END IF
C
C*    2.5 OUTPUT OPTION.
C         --------------
C
      IF (NOUTT.GT.0) THEN
         DO 2501 J=1,NOUTT
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
 2501    CONTINUE
      ELSE
         IF ((FFLAG20.OR.PFLAG20.OR.GFLAG20) .AND.
     1        MOD(IDELINT,IDELPRO).NE.0) THEN
            WRITE(IU06,*) '*******************************************'
            WRITE(IU06,*) '*                                         *'
            WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
            WRITE(IU06,*) '*    ==========================           *'
            WRITE(IU06,*) '* OUTPUT OF INTEGRATED DATA (TOTAL SEA)   *'
            WRITE(IU06,*) '* IS REQUESTED.                           *'
            WRITE(IU06,*) '* OUTPUT TIME STEP HAS TO BE A MULTIPLE   *'
            WRITE(IU06,*) '* OF THE PROPAGATION TIME STEP.           *'
            WRITE(IU06,*) '* OUTPUT TIME STEP IS      IDELINT = ',
     1                    IDELINT
            WRITE(IU06,*) '* PROPAGATION TIME STEP IS IDELPRO = ',
     1                    IDELPRO
            WRITE(IU06,*) '*                                         *'
            WRITE(IU06,*) '*******************************************'
            LERROR = .TRUE.
         END IF
         IF ((FFLAG21.OR.PFLAG21.OR.GFLAG21) .AND.
     1        MOD(IDELINS,IDELPRO).NE.0) THEN
            WRITE(IU06,*) '*******************************************'
            WRITE(IU06,*) '*                                         *'
            WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
            WRITE(IU06,*) '*    ==========================           *'
            WRITE(IU06,*) '* OUTPUT OF INTEGRATED DATA (SEA + SWELL) *'
            WRITE(IU06,*) '* IS REQUESTED.                           *'
            WRITE(IU06,*) '* OUTPUT TIME STEP HAS TO BE A MULTIPLE   *'
            WRITE(IU06,*) '* OF THE PROPAGATION TIME STEP.           *'
            WRITE(IU06,*) '* OUTPUT TIME STEP IS      IDELINS = ',
     1                    IDELINS
            WRITE(IU06,*) '* PROPAGATION TIME STEP IS IDELPRO = ',
     1                    IDELPRO
            WRITE(IU06,*) '*                                         *'
            WRITE(IU06,*) '*******************************************'
            LERROR = .TRUE.
         END IF
         IF ((FFLAG25.OR.PFLAG25.OR.GFLAG25) .AND.
     2        MOD(IDELSPT,IDELPRO).NE.0) THEN
            WRITE(IU06,*) '*******************************************'
            WRITE(IU06,*) '*                                         *'
            WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
            WRITE(IU06,*) '*    ==========================           *'
            WRITE(IU06,*) '* OUTPUT OF SPECTRA (TOTAL SEA)           *'
            WRITE(IU06,*) '* IS REQUESTED.                           *'
            WRITE(IU06,*) '* OUTPUT TIME STEP HAS TO BE A MULTIPLE   *'
            WRITE(IU06,*) '* OF THE PROPAGATION TIME STEP.           *'
            WRITE(IU06,*) '* OUTPUT TIME STEP IS      IDELSPT = ',
     1                     IDELSPT
            WRITE(IU06,*) '* PROPAGATION TIME STEP IS IDELPRO = ',
     1                     IDELPRO
            WRITE(IU06,*) '*                                         *'
            WRITE(IU06,*) '*******************************************'
            LERROR = .TRUE.
         END IF
         IF ((FFLAG26.OR.PFLAG26.OR.GFLAG26) .AND.
     1        MOD(IDELSPS,IDELPRO).NE.0) THEN
            WRITE(IU06,*) '*******************************************'
            WRITE(IU06,*) '*                                         *'
            WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
            WRITE(IU06,*) '*    ==========================           *'
            WRITE(IU06,*) '* OUTPUT OF SPECTRA (SEA + SWELL)         *'
            WRITE(IU06,*) '* IS REQUESTED.                           *'
            WRITE(IU06,*) '* OUTPUT TIME STEP HAS TO BE A MULTIPLE   *'
            WRITE(IU06,*) '* OF THE PROPAGATION TIME STEP.           *'
            WRITE(IU06,*) '* OUTPUT TIME STEP IS      IDELSPS = ',
     1                     IDELSPS
            WRITE(IU06,*) '* PROPAGATION TIME STEP IS IDELPRO = ',
     1                     IDELPRO
            WRITE(IU06,*) '*                                         *'
            WRITE(IU06,*) '*******************************************'
            LERROR = .TRUE.
         ENDIF
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    3. ERROR CHECK.
C        ------------
C
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
C
C ----------------------------------------------------------------------
C
      END
