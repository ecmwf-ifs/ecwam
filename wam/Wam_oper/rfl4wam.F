! ----------------------------------------------------------------------

      PROGRAM RFL4WAM

! ----------------------------------------------------------------------

!**** *RFL4WAM* - PREPARES QUALITY CHECKED SATELLITE DATA FOR WAMASSI.

!     B. HANSEN      ECMWF         DECEMBER 1997

!     MODIFIED BY:    NONE.

!*    PURPOSE.
!     --------

!       SPEED UP ASSIMILATION BY TUNING THE I/O OF OBSERVATIONS.

!**   INTERFACE.
!     ----------

!          ---- FORMAL PARAMETERS ----

!    *rfl4wam*  -d date [ -d date ... ] [ -i itest ] [ -w idelalt ]
!      *date*      STRING   DATE/TIME GROUP YYYYMMDDHHMM.
!      *itest*     INTEGER  TEST OUTPUT LEVEL (*COMMON* *TESTO*).
!      *idelalt*   INTEGER  LENGTH OF ASSIMILATION PERIOD (SECONDS).

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
!          DIRECTORY. IF ICE INFORMATION IS SUPPLIED TO THE WAVE MODEL 
!          IT HAS ALSO TO BE PRESENT HERE IN JUST THE SAME WAY!  THE TWO
!          OUTPUT FILES FROM PREPROC WILL BE ATTACHED TO IU07 AND IU08.
!          NO FURTHER USER INPUT IS REQUIRED. FOR EACH FILE RFL$DATE ONE
!          CORRESPONDING FILE rfl$DATE WILL BE GENERATED CONTAINING ONLY
!          TWO RECORDS:
!          A) THE CHARACTER STRING "CDTPRO"
!          B) THE TWO ARRAYS WHME AND USME 

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
!          "makegrid.F"
!          "micep.F"
!          "mpbcastgrdfld.F"
!          "mpbcastgrid.F"
!          "mpclose_unit.F"
!          "mpdecomp.F"
!          "outpp.F"
!          "readice.F"
!          "readpre.F"
!          "readsat.F"
!          "readt.F"

!     REFERENCE
!     ---------

!          SUBROUTINE WAMASSI.

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWCURR  , ONLY : U
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL      ,IJLT
      USE YOWICE   , ONLY : LICERUN  ,ICEMASK
      USE YOWINTP  , ONLY : CDG
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

      CHARACTER*24 FILNM
      CHARACTER*80 LOGFILENAME 

      REAL Z4(2)
      INTEGER I4(2)

      CHARACTER (LEN=12) , ALLOCATABLE :: CLDATES(:), CLMD(:)

      REAL,ALLOCATABLE,DIMENSION(:,:) :: USME, WHME
      REAL, ALLOCATABLE :: CD(:,:)

! ----------------------------------------------------------------------

!*    1. INITIAL VALUES SET AND CRACK COMMAND LINE.
!        -----------------------------------------
      ITEST   =     0
      IDATES  =     0
      IMO     =    10
      IDELALT = 21600

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
      NPRECR = N_PRECISION(Z4)
      I4=1
      NPRECI = N_PRECISION(I4)

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

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   PRG. RFL4WAM: UNITS DEFINED'
      ENDIF

! ----------------------------------------------------------------------

!*    3. INPUT FROM PREPROCESSING PROGRAMS.
!        ----------------------------------

!     IT IS NOW DONE IN CALL TO MPDECOMP

      ITESTB = MIN(ITESTB,IGL)

! ----------------------------------------------------------------------

!*    4. READ COMMON ICE.
!        ----------------

      LICERUN=.FALSE.

      FILNM = 'sfctempin'
      LIU   = LEN_TRIM(FILNM)
      FILNM=FILNM(1:LIU)
      INQUIRE(FILE=FILNM,EXIST=LICERUN)

      ALLOCATE(ICEMASK(NIBLO,NBLO))

      IF (LICERUN) THEN
        YCLASS = 'od'
        CALL SETWGRIBHD
        CALL READICE (IU06,FILNM,ITEST,NSTART,NEND,IREAD,KTAG)
      ELSE
        ICEMASK=1
        WRITE(IU06,*)'    PRG. RFL4WAM: FILE ',FILNM,' NOT FOUND '
        WRITE(IU06,*)'    PRG. RFL4WAM: NO ICE SHEET USED '
      ENDIF

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

      ALLOCATE(WHME(NGX,NGY), USME(NGX,NGY))
      ALLOCATE(CDG(NGX,NGY))
      CDG     =   0.0


!     DEFINE CDG AS A LAND SEA MASK AS USED BY GRDATA

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
! #include "makegrid.F"
! #include "micep.F"
! #include "mpbcastgrdfld.F"
! #include "mpbcastgrid.F"
! #include "mpclose_unit.F"
! #include "mpdecomp.F"
! #include "outpp.F"
! #include "readice.F"
! #include "readpre.F"
! #include "readsat.F"
! #include "readt.F"
! #include "setwgribhd.F"
