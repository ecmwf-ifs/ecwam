C #include "wave_module.F"
C
C ----------------------------------------------------------------------
C
      PROGRAM RFL4WAM

C ----------------------------------------------------------------------
C
C**** *RFL4WAM* - PREPARES QUALITY CHECKED SATELLITE DATA FOR WAMASSI.
C
C     B. HANSEN      ECMWF         DECEMBER 1997
C
C     MODIFIED BY:    NONE.
C
C*    PURPOSE.
C     --------
C
C       SPEED UP ASSIMILATION BY TUNING THE I/O OF OBSERVATIONS.
C
C**   INTERFACE.
C     ----------
C
C          ---- FORMAL PARAMETERS ----
C
C    *rfl4wam*  -d date [ -d date ... ] [ -i itest ] [ -w idelalt ]
C      *date*      STRING   DATE/TIME GROUP YYMMDDHHMM.
C      *itest*     INTEGER  TEST OUTPUT LEVEL (*COMMON* *TESTO*).
C      *idelalt*   INTEGER  LENGTH OF ASSIMILATION PERIOD (SECONDS).
C
C          ---- INPUT/OUTPUT UNITS ---
C
C           *IU00*   - ERROR STANDARD OUTPUT.
C           *IU05*   - USER INPUT UNIT.
C           *IU06*   - PRINTER OUTPUT.
C           *IU07*   - INPUT  UNIT OF PRECOMPUTED COMMON BLOCKS.
C                      (OUTPUT OF PREPROC).
C           *IU08*   - INPUT  UNIT OF COMMON UBUF.
C                      (OUTPUT OF PREPROC).
C           *IU90*   - INPUT UNIT OF SEA SURFACE TEMPERATURE
C           *IUME*   - ALTIMETER DATA INPUT UNIT (FILES: RFL$DATE).
C           *IO"     - GRIDDED ALTIMETER DATA OUTPUT UNIT
C                       (FILES: rfl$DATE).
C
C          MESSAGE PASSING NOT CATERED FOR YET.
C          IF MESSAGE PASSING ENABLED THE UNITS WILL BE REASSIGNED.
C
C     METHOD.
C     -------
C          READ SATELLITE DATA AND DISTRIBUTE ONTO MODEL GRID.
C          FOR EACH COMMAND LINE OPTION '-d date' A CORRESPONDING FILE
C          NAMED RFL$DATE IS EXPECTED TO BE ON CURRENT WORKING
C          DIRECTORY. IF ICE INFORMATION IS SUPPLIED TO THE WAVE MODEL 
C          IT HAS ALSO TO BE PRESENT HERE IN JUST THE SAME WAY!  THE TWO
C          OUTPUT FILES FROM PREPROC WILL BE ATTACHED TO IU07 AND IU08.
C          NO FURTHER USER INPUT IS REQUIRED. FOR EACH FILE RFL$DATE ONE
C          CORRESPONDING FILE rfl$DATE WILL BE GENERATED CONTAINING ONLY
C          TWO RECORDS:
C          A) THE CHARACTER STRING "CDTPRO"
C          B) THE TWO ARRAYS WHME AND USME 
C
C     EXTERNALS.
C     ----------
C          "wave_module.F"
C
C          "abort1.F"
C          "adjust.F"
C          "confile.F"
C          "cormars.F"
C          "difdate.F"
C          "dummy_mp.F"
C          "expand_string.F"
C          "getclo.F"
C          "grdata.F"
C          "i_get_unit.F"
C          "iecf_len.F"
C          "incdate.F"
C          "inmars.F"
C          "makegrid.F"
C          "micep.F"
C          "mpbcastgrdfld.F"
C          "mpbcastgrid.F"
C          "mpclose_unit.F"
C          "mpdecomp.F"
C          "outpp.F"
C          "readice.F"
C          "readpre.F"
C          "readsat.F"
C          "readt.F"
C
C     REFERENCE
C     ---------
C
C          SUBROUTINE WAMASSI.
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
#include "comcoup.h"
C
#include "txtcoup.h"
C
#include "comcurr.h"
C
#include "txtcurr.h"
C
#include "comfred.h"
C
#include "txtfred.h"
C
#include "comgrid.h"
C
#include "txtgrid.h"
C
#include "comice.h"
C
#include "txtice.h"
C
#include "comindn.h"
C
#include "txtindn.h"
C
#include "commap.h"
C
#include "txtmap.h"
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
#include "commesspass.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comintp.h"
C
#include "txtintp.h"
C
      INTEGER GETCLO, GETCLA, IOPTVAL
      CHARACTER*  3 CLL1
      CHARACTER*  6 CLOPTS
      CHARACTER* 12 CLFMT
      CHARACTER*128 CLARG
      CHARACTER     CLOPTLET
      DATA CLOPTS/'d;i;a;'/
C
      CHARACTER*10 ZERO, CDUM
      CHARACTER*24 FILNM
      CHARACTER*80 LOGFILENAME 
C
      REAL Z4(2)
      INTEGER I4(2)
C
      CHARACTER (LEN=10) , ALLOCATABLE :: CLDATES(:), CLMD(:)
C
      REAL,ALLOCATABLE,DIMENSION(:,:) :: USME, WHME
      REAL, ALLOCATABLE :: CD(:,:)
C
C ----------------------------------------------------------------------
C
C*    1. INITIAL VALUES SET AND CRACK COMMAND LINE.
C        -----------------------------------------
      ITEST   =     0
      IDATES  =     0
      IMO     =    10
      IDELALT = 21600
      CDG     =   0.0
C
      ALLOCATE (CLDATES(IMO))
C
      CMDLINE: DO
        IOPTVAL=GETCLO(clopts,clarg)
        IF (IOPTVAL .LE. 0 )  THEN
          IF (ITEST.GE.1)WRITE (*,'("IOPTVAL = ", I3)') IOPTVAL
          EXIT CMDLINE
        ENDIF

        CLOPTLET=CHAR(IOPTVAL)
        IF (ITEST.GE.1) WRITE(*,*) ' OPTION = ', CLOPTLET

C       GETS VARIABLE ARGUMENT FOR OPTION
        MORARG=GETCLA(clarg)
        IF (ITEST.GE.1) WRITE(*,*) ' MORARG : ', MORARG
        IF (MORARG.NE.0) THEN
          IF ( CLOPTLET .EQ. 'i' ) THEN
C
C ITEST.
C ------
            I1=LEN_TRIM(CLARG)
            WRITE (CLL1,'(I3)') I1
            CLFMT = '(I'//CLL1//')'
            READ (CLARG(1:I1),FMT=CLFMT) ITEST
          ELSEIF ( CLOPTLET .EQ. 'a' ) THEN
C
C IDELALT.
C --------
            I1=LEN_TRIM(CLARG)
            WRITE (CLL1,'(I3)') I1
            CLFMT = '(I'//CLL1//')'
            READ (CLARG(1:I1),FMT=CLFMT) IDELALT
          ELSEIF ( CLOPTLET .EQ. 'd' ) THEN
C
C DATES ARRAY.
C ------------
            IDATES = IDATES + 1
            IF (IDATES.GT. IMO) THEN
              ALLOCATE(CLMD(IMO))
              CLMD = CLDATES
              DEALLOCATE(CLDATES)
              IMO = IMO + 10
              ALLOCATE(CLDATES(IMO))
              DO JL = 1, IMO - 10
                CLDATES(JL) = CLMD(JL)
              ENDDO
              DEALLOCATE(CLMD)
            ENDIF
            CLDATES(IDATES) = CLARG
          ENDIF
        ENDIF
      ENDDO CMDLINE
      IF (ITEST.GE.1) THEN
        WRITE(*,'("NUMBER OF FILES TO PROCESS: ",I2)') IDATES
        WRITE(*,'(5a11)') (CLDATES(JL),JL=1,IDATES)
      ENDIF
      IF (IDATES .EQ. 0 ) THEN
        WRITE(*,*) 'NOTHING TO DO?'
        WRITE(*,*) 'YOU MAY HAVE GIVEN THE WRONG COMMAND LINE OPTIONS:'
        WRITE(*,*) 'USAGE: rfl4wam -d date [ -d date [ ... ] ] ',
     1             '[ -i itest ] [ -w idelalt ] '
        STOP
      ENDIF
C
C ---------------------------------------------------------------------
C
C     2. GLOBAL SETTINGS.
C        ----------------
C
C DETERMINE BYTE STORAGE REPRESENTATION OF REAL NUMBERS.
C ------------------------------------------------------
      Z4=1.
      NPRECR = N_PRECISION(Z4)
      I4=1
      NPRECI = N_PRECISION(I4)

      IF(LMESSPASS) THEN
       CALL MPE_INIT(LDTRMSG,KERROR)
       IF(KERROR.LT.0) THEN
         CALL ABORT1
       ELSE
         IRANK = MPE_MYRANK()
         NPROC = MPE_NPROC()
       ENDIF
C
      ELSE
         IRANK = 1
         NPROC = 1
      ENDIF
C
      NPREVIOUS=IRANK-1
      IF(IRANK.EQ.NPROC) THEN
        NNEXT=0
      ELSE
        NNEXT=IRANK+1
      ENDIF
C
      ALLOCATE (NSTART(NPROC),NEND(NPROC),KLENBOT(NPROC),KLENTOP(NPROC))
C
      NPR=NPROC
      CALL MPDECOMP(NPR,NSTART,NEND,KLENTOP,KLENBOT,MAXLEN)
      MPMAXLENGTH=MAXLEN
      NINF=NSTART(IRANK)-KLENBOT(IRANK)
      NSUP=NEND(IRANK)+KLENTOP(IRANK)
      IF(.NOT.LMESSPASS) NSUP=NIBLO
C
C SET IREAD TO THE PROCESSOR WHICH WILL READ THE INPUT FILES.
C -----------------------------------------------------------
      IREAD=1
C
C DEFINITION OF MODEL PARAMETERS.
C -------------------------------
C
      ZERO = ' '
      LMESSPASS = .FALSE.
      LWCOU = .FALSE.
C
C STANDARD OUTPUT UNIT.
C ---------------------
      IF(LMESSPASS.AND.LWCOU) THEN
        IU06=20
      ELSE IF (LMESSPASS.AND..NOT.LWCOU) THEN
        IU06=66
        LOGFILENAME='logfile.%p'
        CALL EXPAND_STRING(IRANK,NPROC,0,0,LOGFILENAME,1) 
        OPEN(IU06,FILE=LOGFILENAME,STATUS='UNKNOWN')
      ELSE
        IU06=6
      ENDIF
C
C GRID AND UBUF FILES UNITS.
C --------------------------
      IF(LWCOU) THEN
        IU07 = 29
        IU08 = 31
      ELSE
        IU07 = 7
        IU08 = 8
      ENDIF
C
C SEA SURFACE TEMPERATURE FILE UNIT.
C ----------------------------------
      IF(LWCOU) THEN
        IU90 = 38
      ELSE
        IU90 = 90
      ENDIF
C
      IF (ITEST.GE.2) THEN
         WRITE(IU06,*) '   PRG. RFL4WAM: UNITS DEFINED'
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    3. INPUT FROM PREPROCESSING PROGRAMS.
C        ----------------------------------
C
      IF((LMESSPASS.AND.IRANK.EQ.IREAD).OR.(.NOT.LMESSPASS)) THEN
         CALL READPRE (IU07, IREFRA)
         CLOSE (UNIT=IU07, STATUS='KEEP')
      ENDIF
      IF(LMESSPASS) CALL MPBCASTGRID(IU06,IREAD,KTAG,NSTART,NEND,IREFRA)

      IF (ITEST.GE.2) THEN
         IF (IU07 .LT. 111) THEN
           WRITE(IU06,*) '    PRG. RFL4WAM: PREPROC OUTPUT READ ',
     1                   'AND FILE CLOSE UNIT IU07 = ', IU07
         ELSE
           WRITE(IU06,'(2a,a3)') '    PRG. RFL4WAM: PREPROC OUTPUT ',
     1                   'READ AND FILE CLOSE UNIT IU07 = ', IU07
         ENDIF
      ENDIF
      ITESTB = MIN(ITESTB,IGL)
C
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
C SPECIFY THE LIMITS OF THE SEGMENT OVER WHICH THE PROCESS HAS DIRECT
C ACCESS
           IJS(1)=NSTART(IRANK)
           IJL(1)=NEND(IRANK)
           IJLT(1)=NSUP
         ENDIF
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    4. READ COMMON ICE.
C        ----------------
C
      LICERUN=.FALSE.
C
      FILNM = 'sfctempin'
      LIU   = IECF_LEN(FILNM)
      FILNM=FILNM(1:LIU)
      INQUIRE(FILE=FILNM,EXIST=LICERUN)

      IF (LICERUN) THEN
        CALL READICE (IU06,FILNM,ITEST,NSTART,NEND,IREAD,KTAG)
      ELSE
        ICEMASK=1
        WRITE(IU06,*)'    PRG. RFL4WAM: FILE ',FILNM,' NOT FOUND '
        WRITE(IU06,*)'    PRG. RFL4WAM: NO ICE SHEET USED '
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    5. PRINT INITIAL CONDITIONS AS READ FROM PERPROCESSING.
C        ----------------------------------------------------
C
      IF (ITEST .GE. 1) THEN
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
        WRITE(IU06,*) '  '
 3002   FORMAT(3x,a,f8.2,a)
 3003   FORMAT(3x,a,i6,  a)
      ENDIF

      ALLOCATE(WHME(NGX,NGY), USME(NGX,NGY))

C     DEFINE CDG AS A LAND SEA MASK AS USED BY GRDATA

      IF(.NOT.ALLOCATED(CD)) ALLOCATE(CD(NIBLO,NBLO))
      CD = 0.0
      IF(LMESSPASS) THEN
        IG=1
        CALL MAKEGRID (CD, CDG, 1, NEND(NPROC), IG)
      ELSE
        DO IG=1,IGL
          CALL MAKEGRID (CD, CDG, IJS(IG), IJL(IG), IG)
        ENDDO
      ENDIF
      DEALLOCATE(CD)

      DO JLD = 1, IDATES
        CDTPRO = CLDATES(JLD)
        io = i_get_unit(IU06, 'rfl'//CDTPRO, 'w', 'u', 0)
        WRITE(io)CDTPRO, 2
        CALL GRDATA (WHME, USME, NSTART, NEND)
        WRITE (io) WHME, USME
        CLOSE (io)
      ENDDO
C
      END
C 
C EXTERNALS. 
C #include "abort1.F"
C #include "adjust.F"
C #include "confile.F"
C #include "cormars.F"
C #include "difdate.F"
C #include "dummy_mp.F"
C #include "expand_string.F"
C #include "getclo.F"
C #include "grdata.F"
C #include "i_get_unit.F"
C #include "iecf_len.F"
C #include "incdate.F"
C #include "inmars.F"
C #include "makegrid.F"
C #include "micep.F"
C #include "mpbcastgrdfld.F"
C #include "mpbcastgrid.F"
C #include "mpclose_unit.F"
C #include "mpdecomp.F"
C #include "outpp.F"
C #include "readice.F"
C #include "readpre.F"
C #include "readsat.F"
C #include "readt.F"
