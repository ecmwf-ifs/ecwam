      SUBROUTINE USERIN (IU05, IBOUNC, IBOUNF, IFORCA)

C ----------------------------------------------------------------------
C
C**** *USERIN* - ROUTINE TO READ AND WRITE USER INPUT.
C
C     H. GUNTHER     GKSS/ECMWF     NOVEMBER 1989
C
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
C       *ABORT*     - TERMINATES PROCESSING.
C       *DIFDATE*   - COMPUTES A TIME DIFFERENCE.
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
#include "comcout.h"
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
#include "comgrib.h"
C
C ----------------------------------------------------------------------
C
      LOGICAL LERROR
      CHARACTER LINE*72, HEADER*70
C
C ----------------------------------------------------------------------
C
C*    1. READ USER INPUT.
C        ----------------
C
      GFLAG=.F.
      LPACK=.F.
      LFDB =.F.
      ICOUNT = 0
      NOUTT = 0
 1000 CONTINUE
      READ(IU05, '(A72)',ERR=4000,END=4000,IOSTAT=IOS) LINE
      IF (LINE(1:1).EQ.'C') GOTO 1000
      ICOUNT = ICOUNT + 1
      IF (ICOUNT.EQ.1) THEN
         HEADER = LINE(1:70)
         GOTO 1000
      ELSE IF (ICOUNT.EQ.2) THEN
         READ(LINE( 2:11),'(I10)',ERR=4100,IOSTAT=IOS) IDATEA
         READ(LINE(14:23),'(I10)',ERR=4100,IOSTAT=IOS) IDATEE
         READ(LINE(26:35),'(I10)',ERR=4100,IOSTAT=IOS) IDATEF
         GOTO 1000
      ELSE IF (ICOUNT.EQ.3) THEN
         READ(LINE( 2: 8),'(I7)',ERR=4100,IOSTAT=IOS) IDELPRO
         IF (LINE(10:10).EQ.'H') IDELPRO = IDELPRO*3600
         READ(LINE(13:19),'(I7)',ERR=4100,IOSTAT=IOS) IDELT
         IF (LINE(21:21).EQ.'H') IDELT = IDELT*3600
         READ(LINE(24:30),'(I7)',ERR=4100,IOSTAT=IOS) IDELWO
         IF (LINE(32:32).EQ.'H') IDELWO = IDELWO*3600
         READ(LINE(35:41),'(I7)',ERR=4100,IOSTAT=IOS) IDELWI
         IF (LINE(43:43).EQ.'H') IDELWI = IDELWI*3600
         GOTO 1000
      ELSE IF (ICOUNT.EQ.4) THEN
         READ(LINE( 2: 8),'(I7)',ERR=4100,IOSTAT=IOS) IDELINT
         IF (LINE(10:10).EQ.'H') IDELINT = IDELINT*3600
         READ(LINE(13:19),'(I7)',ERR=4100,IOSTAT=IOS) IDELINS
         IF (LINE(21:21).EQ.'H') IDELINS = IDELINS*3600
         READ(LINE(24:30),'(I7)',ERR=4100,IOSTAT=IOS) IDELSPT
         IF (LINE(32:32).EQ.'H') IDELSPT = IDELSPT*3600
         READ(LINE(35:41),'(I7)',ERR=4100,IOSTAT=IOS) IDELSPS
         IF (LINE(43:43).EQ.'H') IDELSPS = IDELSPS*3600
         READ(LINE(46:52),'(I7)',ERR=4100,IOSTAT=IOS) IDELRES
         IF (LINE(54:54).EQ.'H') IDELRES = IDELRES*3600
         GOTO 1000
      ELSE IF (ICOUNT.EQ.5) THEN
         IF (LINE( 2: 4).NE.'END') THEN
            NOUTT = NOUTT+1
            IF (NOUTT.LE.MOUTT)
     1        READ(LINE( 2:11),'(I10)',ERR=4100,IOSTAT=IOS) IOUTT(NOUTT)
            ICOUNT = ICOUNT-1
         ELSE
            GOTO 1000
         ENDIF
         IF (LINE(14:16).NE.'END') THEN
            NOUTT = NOUTT+1
            IF (NOUTT.LE.MOUTT)
     1        READ(LINE(14:23),'(I10)',ERR=4100,IOSTAT=IOS) IOUTT(NOUTT)
         ELSE
            ICOUNT = ICOUNT + 1
            GOTO 1000
         ENDIF
         IF (LINE(26:28).NE.'END') THEN
            NOUTT = NOUTT+1
            IF (NOUTT.LE.MOUTT)
     1        READ(LINE(26:35),'(I10)',ERR=4100,IOSTAT=IOS) IOUTT(NOUTT)
         ELSE
            ICOUNT = ICOUNT + 1
            GOTO 1000
         ENDIF
         IF (LINE(38:40).NE.'END') THEN
            NOUTT = NOUTT+1
            IF (NOUTT.LE.MOUTT)
     1        READ(LINE(38:47),'(I10)',ERR=4100,IOSTAT=IOS) IOUTT(NOUTT)
         ELSE
            ICOUNT = ICOUNT + 1
            GOTO 1000
         ENDIF
         IF (LINE(50:52).NE.'END') THEN
            NOUTT = NOUTT+1
            IF (NOUTT.LE.MOUTT)
     1        READ(LINE(50:59),'(I10)',ERR=4100,IOSTAT=IOS) IOUTT(NOUTT)
         ELSE
            ICOUNT = ICOUNT + 1
            GOTO 1000
         ENDIF
         IF (LINE(62:64).NE.'END') THEN
            NOUTT = NOUTT+1
            IF (NOUTT.LE.MOUTT)
     1        READ(LINE(62:71),'(I10)',ERR=4100,IOSTAT=IOS) IOUTT(NOUTT)
         ELSE
            ICOUNT = ICOUNT + 1
         ENDIF
         GOTO 1000
      ELSE IF (ICOUNT.EQ.6) THEN
         IF (LINE( 2: 2).EQ.'Y') PFLAG( 1) = .TRUE.
         IF (LINE( 4: 4).EQ.'Y') FFLAG( 1) = .TRUE.
         IF (LINE( 6: 6).EQ.'Y') GFLAG( 1) = .TRUE.
         IF (LINE(40:40).EQ.'Y') PFLAG( 2) = .TRUE.
         IF (LINE(42:42).EQ.'Y') FFLAG( 2) = .TRUE.
         IF (LINE(44:44).EQ.'Y') GFLAG( 2) = .TRUE.
         READ(IU05, '(A72)',ERR=4000,END=4000,IOSTAT=IOS) LINE
         IF (LINE( 2: 2).EQ.'Y') PFLAG( 3) = .TRUE.
         IF (LINE( 4: 4).EQ.'Y') FFLAG( 3) = .TRUE.
         IF (LINE( 6: 6).EQ.'Y') GFLAG( 3) = .TRUE.
         IF (LINE(40:40).EQ.'Y') PFLAG( 4) = .TRUE.
         IF (LINE(42:42).EQ.'Y') FFLAG( 4) = .TRUE.
         IF (LINE(44:44).EQ.'Y') GFLAG( 4) = .TRUE.
         READ(IU05, '(A72)',ERR=4000,END=4000,IOSTAT=IOS) LINE
         IF (LINE( 2: 2).EQ.'Y') PFLAG( 5) = .TRUE.
         IF (LINE( 4: 4).EQ.'Y') FFLAG( 5) = .TRUE.
         IF (LINE( 6: 6).EQ.'Y') GFLAG( 5) = .TRUE.
         IF (LINE(40:40).EQ.'Y') PFLAG( 6) = .TRUE.
         IF (LINE(42:42).EQ.'Y') FFLAG( 6) = .TRUE.
         IF (LINE(44:44).EQ.'Y') GFLAG( 6) = .TRUE.
         READ(IU05, '(A72)',ERR=4000,END=4000,IOSTAT=IOS) LINE
         IF (LINE( 2: 2).EQ.'Y') PFLAG( 7) = .TRUE.
         IF (LINE( 4: 4).EQ.'Y') FFLAG( 7) = .TRUE.
         IF (LINE( 6: 6).EQ.'Y') GFLAG( 7) = .TRUE.
         IF (LINE(40:40).EQ.'Y') PFLAG( 8) = .TRUE.
         IF (LINE(42:42).EQ.'Y') FFLAG( 8) = .TRUE.
         IF (LINE(44:44).EQ.'Y') GFLAG( 8) = .TRUE.
         READ(IU05, '(A72)',ERR=4000,END=4000,IOSTAT=IOS) LINE
         IF (LINE( 2: 2).EQ.'Y') PFLAG( 9) = .TRUE.
         IF (LINE( 4: 4).EQ.'Y') FFLAG( 9) = .TRUE.
         IF (LINE( 6: 6).EQ.'Y') GFLAG( 9) = .TRUE.
         IF (LINE(40:40).EQ.'Y') PFLAG(10) = .TRUE.
         IF (LINE(42:42).EQ.'Y') FFLAG(10) = .TRUE.
         IF (LINE(44:44).EQ.'Y') GFLAG(10) = .TRUE.
         READ(IU05, '(A72)',ERR=4000,END=4000,IOSTAT=IOS) LINE
         IF (LINE( 2: 2).EQ.'Y') PFLAG(11) = .TRUE.
         IF (LINE( 4: 4).EQ.'Y') FFLAG(11) = .TRUE.
         IF (LINE( 6: 6).EQ.'Y') GFLAG(11) = .TRUE.
         IF (LINE(40:40).EQ.'Y') PFLAG(12) = .TRUE.
         IF (LINE(42:42).EQ.'Y') FFLAG(12) = .TRUE.
         IF (LINE(44:44).EQ.'Y') GFLAG(12) = .TRUE.
         READ(IU05, '(A72)',ERR=4000,END=4000,IOSTAT=IOS) LINE
         IF (LINE( 2: 2).EQ.'Y') PFLAG(13) = .TRUE.
         IF (LINE( 4: 4).EQ.'Y') FFLAG(13) = .TRUE.
         IF (LINE( 6: 6).EQ.'Y') GFLAG(13) = .TRUE.
         IF (LINE(40:40).EQ.'Y') PFLAG(14) = .TRUE.
         IF (LINE(42:42).EQ.'Y') FFLAG(14) = .TRUE.
         IF (LINE(44:44).EQ.'Y') GFLAG(14) = .TRUE.
         GOTO 1000
      ELSE IF (ICOUNT.EQ.7) THEN
         IF (LINE( 5: 5).EQ.'Y') LFDB      = .TRUE.
         GOTO 1000
      ELSE IF (ICOUNT.EQ.8) THEN
         READ(LINE( 2: 8),'(I7)',ERR=4100,IOSTAT=IOS) ICASE
         READ(LINE(11:17),'(I7)',ERR=4100,IOSTAT=IOS) ISHALLO
         READ(LINE(20:26),'(I7)',ERR=4100,IOSTAT=IOS) IREFRA
         READ(LINE(29:35),'(I7)',ERR=4100,IOSTAT=IOS) ITEST
         READ(LINE(38:44),'(I7)',ERR=4100,IOSTAT=IOS) ITESTB
         READ(LINE(47:53),'(I7)',ERR=4100,IOSTAT=IOS) IREST
         READ(LINE(56:62),'(I7)',ERR=4100,IOSTAT=IOS) IASSI
         READ(LINE(65:71),'(I7)',ERR=4100,IOSTAT=IOS) IFORCA
         GOTO 1000
      ELSE IF (ICOUNT.EQ.9) THEN
         READ(LINE( 2: 8),'(I7)',ERR=4100,IOSTAT=IOS) IBOUNC
         READ(LINE(11:17),'(I7)',ERR=4100,IOSTAT=IOS) IBOUNF
         GOTO 1000
      ELSE IF (ICOUNT.EQ.10) THEN
         USERID = LINE(2:4)
         RUNID  = LINE(7:9)
         PATH   = LINE(12:71)
      ENDIF
C
C*    1.1  READ THE WAMINFO FILE AND OVERWRITE INPUT.
C          ------------------------------------------
C
      CALL READSTA(IU04, IBEGDT, IENDDT, IANALPD, IFOREPD, IDELWIN,
     1             IBPLTDT, IEPLTDT, IASSIM, NFCST, ISTAT)
      IDATEA = IBEGDT
      IDATEF = IBEGDT
      CALL INCDATE (IDATEF,IANALPD*3600)
      IDATEE = IDATEF
      CALL INCDATE (IDATEE,IFOREPD*3600)
      IF (IDATEE.NE.IENDDT) THEN
         WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++'
         WRITE(IU06,*) '+                                           +'
         WRITE(IU06,*) '+    WARNING ERROR IN SUB. USERIN           +'
         WRITE(IU06,*) '+    ============================           +'
         WRITE(IU06,*) '+ END DATE OF RUN AND RUN LENGTH DO NOT     +'
         WRITE(IU06,*) '+ MATCH IN THE WAM INFO FILE.               +'
         WRITE(IU06,*) '+ START DATE OF RUN       IBEGDT  = ', IBEGDT
         WRITE(IU06,*) '+ ANALYSIS PERIOD (HOURS) IANALPD = ', IANALPD
         WRITE(IU06,*) '+ FORECAST PERIOD (HOURS) IFOREPD = ', IFOREPD
         WRITE(IU06,*) '+ END DATE FROM WAMINFO   IENDDT  = ', IENDDT
         WRITE(IU06,*) '+ END DATE COMPUTED       IDATEE  = ', IDATEE
         WRITE(IU06,*) '+                                           +'
         WRITE(IU06,*) '+ MODEL CONTINUES USING COMPUTED END DATE.  +'
         WRITE(IU06,*) '+                                           +'
         WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++'
         IENDDT = IDATEE
      ENDIF
      IDELWI = 3600*IDELWIN
      IASSI = IASSIM
      IFORCA = NFCST
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
C*    1.3  PRINT INITIAL CONDITIONS.
C          -------------------------

      WRITE(IU06,'(1H1)')
      WRITE(IU06,'(1H0,A70)') HEADER
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' STARTING  DATE (FORMAT:YYMMDDHHMM) : ',IDATEA,
     1              ' END DATE :',IDATEE
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' MODEL TIME STEPS:'
      WRITE(IU06,*) ' SOURCE TERM INTEGRATION TIME STEP : ',
     1                IDELT,' SECS'
      WRITE(IU06,*) ' PROPAGATION TIME STEP ............: ',
     1                IDELPRO,' SECS'
      WRITE(IU06,*) ' MODEL WIND INPUT TIME STEP .......: ',
     1                IDELWI,' SECS'
      WRITE(IU06,*) ' MODEL WIND OUTPUT TIME STEP.......: ',
     1                IDELWO,' SECS'
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' MODEL OPTIONS:'
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
         WRITE(IU06,*) ' WIND FIELD ARE NOT INTERPOLATED IN TIME'
      ELSE
         WRITE(IU06,*) ' WIND FIELD ARE INTERPOLATED IN TIME'
      ENDIF
      IF (IBOUNC.EQ.1) THEN
         WRITE(IU06,*) ' MODEL RUNS WITH BOUNDARY POINTS (COURSE GRID)'
      ELSE
         WRITE(IU06,*) ' MODEL RUNS WITHOUT BOUNDARY POINTS',
     1                 ' (COURSE GRID)'
      ENDIF
      IF (IBOUNF.EQ.1) THEN
         WRITE(IU06,*) ' MODEL RUNS WITH BOUNDARY POINTS (FINE GRID)'
      ELSE
         WRITE(IU06,*) ' MODEL RUNS WITHOUT BOUNDARY POINTS',
     1                 ' (FINE GRID)'
      ENDIF
      IF (IFORCA.EQ.1) THEN
         WRITE(IU06,*) ' MODEL STARTS FROM ANALYSIS FIELDS'
      ELSE
         WRITE(IU06,*) ' MODEL STARTS FROM FORECAST FIELDS'
      ENDIF
      IF (IASSI.EQ.1) THEN
         WRITE(IU06,*) ' DATA ASSIMILATION IS CARRIED OUT'
      ELSE
         WRITE(IU06,*) ' DATA ASSIMILATION IS NOT CARRIED OUT'
      ENDIF
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' MODEL OUTPUT SELECTION:'
      WRITE(IU06,*) '  '
      IF (NOUTT.GT.0) THEN
         WRITE(IU06,*) ' NUMBER OF OUTPUT TIMES IS NOUTT = ', NOUTT
         WRITE(IU06,*) ' OUTPUT WILL BE PROCESSED AT:'
         WRITE(IU06,'(6(2X,I10.10))') (IOUTT(I),I=1,NOUTT)
         WRITE(IU06,*) '  '
         WRITE(IU06,*) ' TO PRINTER AND/OR UNIT IU20 = ',IU20
      ELSE
         WRITE(IU06,*) ' TO PRINTER AND/OR UNIT IU20 = ',IU20,
     1                 ' EVERY ', IDELINT,' SECONDS'
      ENDIF
      WRITE(IU06,*) '    F = FALSE   T = TRUE '
      WRITE(IU06,*) '                        PRINTER  UNIT'
      WRITE(IU06,*) ' SIGNIFICANT WAVE HEIGHT ..', PFLAG( 1),
     1              '......', FFLAG( 1)
      WRITE(IU06,*) ' MEAN WAVE DIRECTION ......', PFLAG( 2),
     1              '......', FFLAG( 2)
      WRITE(IU06,*) ' WAVE MEAN FREQUENCY ......', PFLAG( 3),
     1              '......', FFLAG( 3)
      WRITE(IU06,*) ' FRICTION VELOCITY ........', PFLAG( 4),
     1              '......', FFLAG( 4)
      WRITE(IU06,*) ' WIND DIRECTION ...........', PFLAG( 5),
     1              '......', FFLAG( 5)
      WRITE(IU06,*) ' WAVE PEAK FREQUENCY ......', PFLAG( 6),
     1              '......', FFLAG( 6)
      WRITE(IU06,*) ' DRAG COEFFICIENT .........', PFLAG( 7),
     1              '......', FFLAG( 7)
      WRITE(IU06,*) ' NORMALISED WAVE STRESS ...', PFLAG( 8),
     1              '......', FFLAG( 8)
      WRITE(IU06,*) '  '
      IF (NOUTT.GT.0) THEN
         WRITE(IU06,*) ' TO PRINTER AND/OR UNIT IU21 = ',IU21
      ELSE
         WRITE(IU06,*) ' TO PRINTER AND/OR UNIT IU21 = ',IU21,
     1                 ' EVERY ', IDELINS,' SECONDS'
      ENDIF
      WRITE(IU06,*) '                        PRINTER  UNIT'
      WRITE(IU06,*) ' SWELL WAVE HEIGHT ........', PFLAG( 9),
     1              '......', FFLAG( 9)
      WRITE(IU06,*) ' SWELL DIRECTION ..........', PFLAG(10),
     1              '......', FFLAG(10)
      WRITE(IU06,*) ' WIND WAVE DIRECTION ......', PFLAG(11),
     1              '......', FFLAG(11)
      WRITE(IU06,*) ' SWELL MEAN FREQUENCY .....', PFLAG(12),
     1              '......', FFLAG(12)
      WRITE(IU06,*) '  '
      IF (NOUTT.GT.0) THEN
         WRITE(IU06,*) ' TO PRINTER AND/OR UNIT IU25 = ',IU25
      ELSE
         WRITE(IU06,*) ' TO PRINTER AND/OR UNIT IU25 = ',IU25,
     1                 ' EVERY ', IDELSPT,' SECONDS'
      ENDIF
      WRITE(IU06,*) '                        PRINTER  UNIT'
      WRITE(IU06,*) ' SPECTRA OF TOTAL SEA .....', PFLAG(13),
     1              '......', FFLAG(13)
      WRITE(IU06,*) '  '
      IF (NOUTT.GT.0) THEN
         WRITE(IU06,*) ' TO PRINTER AND/OR UNIT IU26 = ',IU26
      ELSE
         WRITE(IU06,*) ' TO PRINTER AND/OR UNIT IU26 = ',IU26,
     1              ' EVERY ', IDELSPS,' SECONDS'
      ENDIF
      WRITE(IU06,*) '                        PRINTER  UNIT'
      WRITE(IU06,*) ' SWELL SPECTRA ............', PFLAG(14),
     1              '......', FFLAG(14)
      WRITE(IU06,*) '  '

      IF (IREST.EQ.1) THEN
        WRITE(IU06,*) ' OUTPUT AND RESTART FILES WILL BE DISPOSED EVERY'
     1               ,'...', IDELRES, ' SECONDS'
      ELSE
        WRITE(IU06,*) ' OUTPUT FILES WILL BE DISPOSED EVERY ...........'
     1               ,'...', IDELRES, ' SECONDS'
      ENDIF
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' MODEL FILE NAMES:'
      WRITE(IU06,*) ' USERID IS ........: ', USERID
      WRITE(IU06,*) ' RUN IDENTIFIER IS : ', RUNID
      WRITE(IU06,*) ' PATH NAME IS .....: ', PATH
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' TEST OUTPUT LEVEL IS .............. ITEST = ',
     1                ITEST
      WRITE(IU06,*) ' TEST OUTPUT OF IN BLOCK LOOPS UPTO ITESTB = ',
     1                ITESTB
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
      FFLAG20 = FFLAG( 1).OR.FFLAG( 2).OR.FFLAG( 3).OR.FFLAG( 4).OR.
     1          FFLAG( 5).OR.FFLAG( 6).OR.FFLAG( 7).OR.FFLAG( 8)
      PFLAG20 = PFLAG( 1).OR.PFLAG( 2).OR.PFLAG( 3).OR.PFLAG( 4).OR.
     1          PFLAG( 5).OR.PFLAG( 6).OR.PFLAG( 7).OR.PFLAG( 8)
      FFLAG21 = FFLAG( 9).OR.FFLAG(10).OR.FFLAG(11).OR.FFLAG(12)
      PFLAG21 = PFLAG( 9).OR.PFLAG(10).OR.PFLAG(11).OR.PFLAG(12)
      FFLAG25 = FFLAG(13)
      PFLAG25 = PFLAG(13)
      FFLAG26 = FFLAG(14)
      PFLAG26 = PFLAG(14)
      IF (NOUTT.GT.0) THEN
         DO 2501 J=1,NOUTT
            CALL DIFDATE (IDATEA, IOUTT(J), ISHIFT)
            IF (ISHIFT.LE.0 .OR. MOD(ISHIFT,IDELPRO).NE.0) THEN
               WRITE(IU06,*) '++++++++++++++++++++++++++++++++++++++++'
               WRITE(IU06,*) '+                                      +'
               WRITE(IU06,*) '+    WARNING ERROR IN SUB. USERIN      +'
               WRITE(IU06,*) '+    ============================      +'
               WRITE(IU06,*) '+ OUTPUT DATE IS NOT AT THE END OF A   +'
               WRITE(IU06,*) '+ PROPAGATION TIMESTEP.                +'
               WRITE(IU06,*) '+ DATE IS : ', IOUTT(J)
               WRITE(IU06,*) '+ PROGRAM WILL IGNORE THIS OUTPUT TIME +'
               WRITE(IU06,*) '+                                      +'
               WRITE(IU06,*) '++++++++++++++++++++++++++++++++++++++++'
            ENDIF
 2501    CONTINUE
      ELSE
         IF ((FFLAG20.OR.PFLAG20) .AND. MOD(IDELINT,IDELPRO).NE.0) THEN
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
         IF ((FFLAG21.OR.PFLAG21) .AND. MOD(IDELINS,IDELPRO).NE.0) THEN
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
         IF ((FFLAG25.OR.PFLAG25) .AND. MOD(IDELSPT,IDELPRO).NE.0) THEN
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
         IF ((FFLAG26.OR.PFLAG26) .AND. MOD(IDELSPS,IDELPRO).NE.0) THEN
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
         CALL ABORT
      ELSE
         RETURN
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    4. READ ERROR MESSAGES.
C        --------------------
C
 4000 CONTINUE
         WRITE(IU06,*) ' ********************************************'
         WRITE(IU06,*) ' *                                          *'
         WRITE(IU06,*) ' *     FATAL ERROR IN SUB. USERIN           *'
         WRITE(IU06,*) ' *     ==========================           *'
         WRITE(IU06,*) ' * READ ERROR ON INPUT FILE:                *'
         WRITE(IU06,*) ' * ERROR IS LATER THAN ICOUNT = ', ICOUNT
         WRITE(IU06,*) ' * LAST LINE READ IS     LINE = ', LINE
         WRITE(IU06,*) ' *                                          *'
         WRITE(IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS         *'
         WRITE(IU06,*) ' *                                          *'
         WRITE(IU06,*) ' ********************************************'
         CALL ABORT
 4100 CONTINUE
         WRITE(IU06,*) ' ********************************************'
         WRITE(IU06,*) ' *                                          *'
         WRITE(IU06,*) ' *     FATAL ERROR IN SUB. USERIN           *'
         WRITE(IU06,*) ' *     ==========================           *'
         WRITE(IU06,*) ' * READ ERROR ON CHARACTER STRING           *'
         WRITE(IU06,*) ' * ERROR IS IN DATA LINE ICOUNT = ', ICOUNT
         WRITE(IU06,*) ' * CHARACTER STRING IS   LINE = ', LINE
         WRITE(IU06,*) ' *                                          *'
         WRITE(IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS         *'
         WRITE(IU06,*) ' *                                          *'
         WRITE(IU06,*) ' ********************************************'
         CALL ABORT

      END
