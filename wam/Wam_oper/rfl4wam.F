! ----------------------------------------------------------------------

      PROGRAM RFL4WAM

! ----------------------------------------------------------------------

!**** *RFL4WAM* - PREPARES QUALITY CHECKED SATELLITE DATA FOR WAMASSI.

!     B. HANSEN      ECMWF         DECEMBER 1997

!*    PURPOSE.
!     --------

!       SPEED UP ASSIMILATION BY TUNING THE I/O OF OBSERVATIONS.
!       AND PEFORMS BLACKLISTING AND CORRECTION OF THE DATA.

!**   INTERFACE.
!     ----------

!          ---- FORMAL PARAMETERS ----

!    *rfl4wam*  -d date [ -d date ... ] [ -i itest ] [ -w idelalt ]
!      *date*      STRING   DATE/TIME GROUP YYYYMMDDHHMM.
!      *itest*     INTEGER  TEST OUTPUT LEVEL (*COMMON* *TESTO*).
!      *idelalt*   INTEGER  LENGTH OF ASSIMILATION PERIOD (SECONDS).
!      if the file wave_data_blacklist is present, information from it
!      will be used to blacklist data.

!          ---- INPUT/OUTPUT UNITS ---

!           *IU00*   - ERROR STANDARD OUTPUT.
!           *IU05*   - USER INPUT UNIT.
!           *IU06*   - PRINTER OUTPUT.
!           *IU07*   - INPUT  UNIT OF PRECOMPUTED COMMON BLOCKS.
!                      (OUTPUT OF PREPROC).
!           *IU08*   - INPUT  UNIT OF COMMON UBUF.
!                      (OUTPUT OF PREPROC).
!           *IU90*   - INPUT UNIT OF SEA SURFACE TEMPERATURE
!           *IUME*   - ALTIMETER DATA INPUT UNIT (FILES: RFL$DATE).
!           *IO"     - GRIDDED ALTIMETER DATA OUTPUT UNIT
!                       (FILES: rfl$DATE).

!          MESSAGE PASSING NOT CATERED FOR YET.
!          IF MESSAGE PASSING ENABLED THE UNITS WILL BE REASSIGNED.

!     METHOD.
!     -------
!          READ SATELLITE DATA AND DISTRIBUTE ONTO MODEL GRID.
!          FOR EACH COMMAND LINE OPTION '-d date' A CORRESPONDING FILE
!          NAMED RFL$DATE IS EXPECTED TO BE ON CURRENT WORKING
!          DIRECTORY. THE TWO
!          OUTPUT FILES FROM PREPROC WILL BE ATTACHED TO IU07 AND IU08.
!          NO FURTHER USER INPUT IS REQUIRED. FOR EACH FILE RFL$DATE ONE
!          CORRESPONDING FILE rfl$DATE WILL BE GENERATED CONTAINING ONLY
!          RECORDS:
!          A) THE CHARACTER STRING "CDTPRO"
!          B) THE ARRAYS IJALT and ALTDATA

!     EXTERNALS.
!     ----------
!          "yowspec.F"

!          "abort1.F"
!          "adjust.F"
!          "confile.F"
!          "cormars.F"
!          "difdate.F"
!          "dummy_mp.F"
!          "expand_string.F"
!          "getclo.F"
!          "grdata.F"
!          "i_get_unit.F"
!          "incdate.F"
!          "inmars.F"
!          "iniwcst.F"
!          "makegrid.F"
!          "mpbcastgrdfld.F"
!          "mpbcastgrid.F"
!          "mpclose_unit.F"
!          "mpdecomp.F"
!          "outpp.F"
!          "readpre.F"
!          "readsat.F"
!          "readt.F"

!     REFERENCE
!     ---------

!          SUBROUTINE WAMASSI.

! ----------------------------------------------------------------------

      USE YOWALTAS , ONLY : IJALT    ,ALTDATA
      USE YOWCOUP  , ONLY : LWCOU
      USE YOWCURR  , ONLY : U
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL      ,IJLT
      USE YOWMAP   , ONLY : AMOWEP   ,AMOSOP   ,AMOEAP   ,AMONOP   ,
     &            XDELLA   ,XDELLO
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPREVIOUS,NNEXT    ,
     &            MPMAXLENGTH,KTAG     ,NPRECR   ,
     &            NPRECI
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NBLO     ,NIBLO
      USE YOWSTAT  , ONLY : CDTPRO   ,IDELALT  ,YCLASS
      USE YOWSPEC  , ONLY : NSTART   ,NEND     ,KLENTOP  ,KLENBOT
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB

! -------------------------------------------------------------------

      INTEGER GETCLO, GETCLA, IOPTVAL
      CHARACTER*  3 CLL1
      CHARACTER*  6 CLOPTS
      CHARACTER* 12 CLFMT
      CHARACTER*128 CLARG
      CHARACTER     CLOPTLET
      DATA CLOPTS/'d;i;a;'/

      REAL Z4(2)
      INTEGER I4(2)

      REAL,ALLOCATABLE,DIMENSION(:,:) :: WHME
      CHARACTER (LEN=12),ALLOCATABLE :: CLDATES(:), CLMD(:)

! ----------------------------------------------------------------------

!*    1. INITIAL VALUES SET AND CRACK COMMAND LINE.
!        -----------------------------------------

      CALL INIWCST

      ITEST   =     0
      IDATES  =     0
      IMO     =    10
      IDELALT = 21600
      IU06 = 6

      ALLOCATE (CLDATES(IMO))

      CMDLINE: DO
        IOPTVAL=GETCLO(clopts,clarg)
        IF (IOPTVAL .LE. 0 )  THEN
          IF (ITEST.GE.1)WRITE (*,'("IOPTVAL = ", I3)') IOPTVAL
          EXIT CMDLINE
        ENDIF

        CLOPTLET=CHAR(IOPTVAL)
        IF (ITEST.GE.1) WRITE(*,*) ' OPTION = ', CLOPTLET

!       GETS VARIABLE ARGUMENT FOR OPTION
        MORARG=GETCLA(clarg)
        IF (ITEST.GE.1) WRITE(*,*) ' MORARG : ', MORARG
        IF (MORARG.NE.0) THEN
          IF ( CLOPTLET .EQ. 'i' ) THEN

! ITEST.
! ------
            I1=LEN_TRIM(CLARG)
            WRITE (CLL1,'(I3)') I1
            CLFMT = '(I'//CLL1//')'
            READ (CLARG(1:I1),FMT=CLFMT) ITEST
          ELSEIF ( CLOPTLET .EQ. 'a' ) THEN

! IDELALT.
! --------
            I1=LEN_TRIM(CLARG)
            WRITE (CLL1,'(I3)') I1
            CLFMT = '(I'//CLL1//')'
            READ (CLARG(1:I1),FMT=CLFMT) IDELALT
          ELSEIF ( CLOPTLET .EQ. 'd' ) THEN

! DATES ARRAY.
! ------------
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
        WRITE(*,'(5a13)') (CLDATES(JL),JL=1,IDATES)
      ENDIF
      IF (IDATES .EQ. 0 ) THEN
        WRITE(*,*) 'NOTHING TO DO?'
        WRITE(*,*) 'YOU MAY HAVE GIVEN THE WRONG COMMAND LINE OPTIONS:'
        WRITE(*,*) 'USAGE: rfl4wam -d date [ -d date [ ... ] ] ',
     &             '[ -i itest ] [ -w idelalt ] '
        STOP
      ENDIF

! ---------------------------------------------------------------------

!     2. GLOBAL SETTINGS.
!        ----------------

      LMESSPASS = .FALSE.

! DETERMINE BYTE STORAGE REPRESENTATION OF REAL NUMBERS.
! ------------------------------------------------------
      Z4=1.
      NPRECR = KIND(Z4)
      I4=1
      NPRECI = KIND(I4)

      IRANK = 1
      NPROC = 1

      NPREVIOUS=IRANK-1
      IF(IRANK.EQ.NPROC) THEN
        NNEXT=0
      ELSE
        NNEXT=IRANK+1
      ENDIF

      ALLOCATE (NSTART(NPROC),NEND(NPROC),KLENBOT(NPROC),KLENTOP(NPROC))

      NPR=NPROC

      CALL MPDECOMP(NPR,NSTART,NEND,KLENTOP,KLENBOT,MAXLEN)

      MPMAXLENGTH=MAXLEN

! SET IREAD TO THE PROCESSOR WHICH WILL READ THE INPUT FILES.
! -----------------------------------------------------------
      IREAD=1

! DEFINITION OF MODEL PARAMETERS.
! -------------------------------

      LWCOU = .FALSE.

! SEA SURFACE TEMPERATURE FILE UNIT.
! ----------------------------------
      IF(LWCOU) THEN
        IU90 = 38
      ELSE
        IU90 = 90
      ENDIF


! ----------------------------------------------------------------------

!*    3. INPUT FROM PREPROCESSING PROGRAMS.
!        ----------------------------------

!     IT IS NOW DONE IN CALL TO MPDECOMP

      ITESTB = MIN(ITESTB,IGL)

! ----------------------------------------------------------------------

!*    5. PRINT INITIAL CONDITIONS AS READ FROM PREPROCESSING.
!        ----------------------------------------------------

      IF (ITEST .GE. 1) THEN
        WRITE(IU06,*) '  '
        WRITE(IU06,*) ' WAVE MODEL GRID ORGANISATION:'
        WRITE(IU06,3002) ' SOUTHERNMOST LATITUDE IN GRID IS .......: ',
     &                 AMOSOP, ' DEGREE'
        WRITE(IU06,3002) ' NORTHERNMOST LATITUDE IN GRID IS .......: ',
     &                 AMONOP, ' DEGREE'
        WRITE(IU06,3002) ' WESTERNMOST LONGITUDE IN GRID IS .......: ',
     &               AMOWEP, ' DEGREE'
        WRITE(IU06,3002) ' EASTERNMOST LONGITUDE IN GRID IS .......: ',
     &               AMOEAP, ' DEGREE'
        WRITE(IU06,3002) ' LATITUDE INCREMENT IS ..................: ',
     &                 XDELLA, ' DEGREE'
        WRITE(IU06,3002) ' LONGITUDE INCREMENT IS .................: ',
     &               XDELLO, ' DEGREE'
        WRITE(IU06,*) '  '
        WRITE(IU06,3003) ' TOTAL NUMBER OF BLOCKS IS...............: ',
     &                   IGL
        WRITE(IU06,3003) ' TOTAL LENGTH OF EACH BLOCK .............: ',
     &                    NIBLO
        WRITE(IU06,*) '  '
 3002   FORMAT(3x,a,f8.2,a)
 3003   FORMAT(3x,a,i6,  a)
      ENDIF

      ALLOCATE(WHME(NGX,NGY))

      DO JLD = 1, IDATES
        CDTPRO = CLDATES(JLD)
        io = i_get_unit(IU06, 'rfl'//CDTPRO, 'w', 'u', 0)

        CALL GRDATA (NOBS, NSTART, NEND)

        WRITE(io) CDTPRO, NOBS

        IF(NOBS.GT.0) THEN
        
!         DEFINE COEFFICIENT FOR THE ALT DATA ERROR
!!! it will need to be adapted for other satellite !!!
          SIGOBS=0.5
          A=SIGOBS
          B=-LOG(0.001)
          HSMIN=1.


          DO IOBS=1,NOBS
            ALTDATA(IOBS,2) = SIGOBS+A*EXP(-B*(ALTDATA(IOBS,1)-HSMIN))
          ENDDO
!!! we will have to adapt for the linear regression correctio
!!!

          WRITE (io) IJALT 
          WRITE (io) ALTDATA

          DEALLOCATE(IJALT)
          DEALLOCATE(ALTDATA)
        ENDIF
        CLOSE (io)
      ENDDO

      END

! EXTERNALS. 
! #include "abort1.F"
! #include "adjust.F"
! #include "confile.F"
! #include "cormars.F"
! #include "difdate.F"
! #include "dummy_mp.F"
! #include "expand_string.F"
! #include "getclo.F"
! #include "grdata.F"
! #include "i_get_unit.F"
! #include "incdate.F"
! #include "inmars.F"
! #include "iniwcst.F"
! #include "makegrid.F"
! #include "mpbcastgrdfld.F"
! #include "mpbcastgrid.F"
! #include "mpclose_unit.F"
! #include "mpdecomp.F"
! #include "outpp.F"
! #include "readpre.F"
! #include "readsat.F"
! #include "setwgribhd.F"
