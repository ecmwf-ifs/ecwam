! ----------------------------------------------------------------------

      PROGRAM COMBINE_ODB_TXT

! ----------------------------------------------------------------------

!**** *COMBINE_ODB_TXT* - COMBINES TEXT FILES FOR USE WITH "simulobs2odb"

!     S. ABDALLA     ECMWF         JUNE 2011

!*    PURPOSE.
!     --------

!       READS ONE OR MORE TEXT FILES (OF SAME TYPE) THAT ARE
!       PREPARED TO GENERATE ODB RALT USING TOOL "simulobs2odb"
!       THEN PRODUCES A SINGLE FILE BY REMOVING THE IRRELEVANT LINES

!**   INTERFACE.
!     ----------

!          ---- FORMAL PARAMETERS ----

!    *combine_odb_txt*  infile1 [ infile2 [ infile3 [ ...] ]] 
!                   
!      *infile*    STRING   INPUT FILE NAME

!          ---- INPUT/OUTPUT UNITS ---

!           *IU00*   - STANDARD ERROR UNIT.
!           *IU06*   - STANDARD OUTPUT UNIT (NEEDS TO BE REDIRECTED TO A FILE).
!           *IU12*   - UNIT FOR INPUT FILES.
!
!     EXIT ERRORS:
!     ------------
!           * 0*     - NORMAL TERMINATION (SUCCESSFUL RUN).
!           * 1*     - ERROR: Too few number of arguments!
!           * 2*     - ERROR: Number of input files exceeds maximum permitted!
!           * 3*     - ERROR: No input files were provided!
!           * 4*     - ERROR: No output file name was provided!
!           * 5*     - ERROR: Input file does not exist!
!           * 6*     - ERROR while opening input file!
!           * 7*     - Unknown input ERROR while reading file!
!
!     METHOD.
!     -------
!          NONE
!
!     EXTERNALS.
!     ----------
!          NONE
!  

      PARAMETER (MAXINFILES=20, MAXARG=MAXINFILES)
      CHARACTER*1024 CLINE
      CHARACTER*512 CARG(0:MAXARG)
      CHARACTER*512 infile(MAXINFILES), outfile
      LOGICAL LOUT, LEXIST, LERR

! ----------------------------------------------------------------------

!*    1. INITIAL VALUES SET AND CRACK COMMAND LINE.
!        -----------------------------------------

      IU00 = 0
      IU06 = 6
      IU12 =12

      LOUT=.FALSE.
      infile=' '
      outfile=' '
      NumInFiles=0

      NARG=IARGC()
      IF(NARG.lt.1) THEN
        CALL GETARG(0,CARG(0))
        WRITE(0,*)'ERROR: Too few number of arguments!'
        WRITE(0,*)'Proper use: '//CARG(0)(1:LEN_TRIM(CARG(0)))//
     &            ' inputfile1 inputfile2 ...'
        CALL EXIT(1)
      ELSEIF(NARG.GT.MAXARG) THEN
        WRITE(0,*)'WARNING: Too many arguments! Some will be ignored!'
      END IF

      DO J=0,NARG
        CALL GETARG(J,CARG(J))
      ENDDO

      DO J=1,NARG
         NumInFiles=NumInFiles+1
         IF (NumInFiles.GT.MAXINFILES) THEN
          WRITE(0,*)'ERROR: '//
     &             'Number of input files exceeds maximum permitted'
          WRITE(0,*)'       Maximum number of input files'
     &                   //' is: ',MAXINFILES
          CALL EXIT(2)
        ENDIF
        infile(NumInFiles)=CARG(J)
      ENDDO

      IF (NumInFiles.GE.1) THEN
        WRITE(0,*)'Number of input files: ', NumInFiles
      ELSE
        WRITE(0,*)'ERROR: No input files were provided!'
        WRITE(0,*)'       Proper use: '//CARG(0)(1:LEN_TRIM(CARG(0)))//
     &            ' inputfile1 inputfile2 ...'
        WRITE(0,*)'       Aborting without doing anything!'
        CALL EXIT(3)
      ENDIF

!*    2. CHECK FOR THE EXISTENCE OF FILES
!        --------------------------------

      LERR=.FALSE.
      DO I=1, NumInFiles
        INQUIRE (FILE=infile(I), EXIST=LEXIST)
        IF (.NOT. LEXIST) THEN
          LERR=.TRUE.
          WRITE(0,*)'ERROR: Input file does not exist: '//
     &                                infile(I)(1:LEN_TRIM(infile(I)))
        ENDIF
      ENDDO
      IF (LERR) THEN
        WRITE(0,*)'       Aborting without doing anything!'
        CALL EXIT(5)
      ENDIF

!*    4. OPEN FIRST INPUT FILE AND SENDS ITS CONTENT TO OUTPUT
!        -----------------------------------------------------
      I=1
      OPEN (IU12 , FILE=infile(I), STATUS='OLD', IOSTAT=IOS)
      IF (IOS.NE.0) THEN
        WRITE(0,*)'ERROR while opening input file: '//
     &                                infile(I)(1:LEN_TRIM(infile(I)))
        CALL EXIT(6)
      ENDIF
      
      DO
        READ(IU12,'(A)',END=100, ERR=900, IOSTAT=IOS) CLINE
        WRITE(IU06,'(A)') CLINE(1:LEN_TRIM(CLINE))
      ENDDO

  100 CLOSE (IU12)

!*    5. LOOP OVER REMAINING INPUT FILES AND COPY RELEVANT CONTENT TO OUTPUT FILE
!        ------------------------------------------------------------------------
      IF (NumInFiles.GT.1) THEN
        DO I=2, NumInFiles
          OPEN (IU12 , FILE=infile(I), STATUS='OLD', IOSTAT=IOS)
          IF (IOS.NE.0) THEN
            WRITE(0,*)'ERROR while opening input file: '//
     &                                infile(I)(1:LEN_TRIM(infile(I)))
            CALL EXIT(6)
          ENDIF
      
          DO
            READ(IU12,'(A)',END=200, ERR=900, IOSTAT=IOS) CLINE
            IF (CLINE.NE.' ') THEN
              READ(CLINE,*,IOSTAT=IOS) A
              IF (IOS.EQ.0) WRITE(IU06,'(A)') CLINE(1:LEN_TRIM(CLINE))
            ENDIF
          ENDDO

  200     CLOSE (IU12)
        ENDDO
      ENDIF

!*    6. TERMINATE SUCCESSFULLY
!        ----------------------

      WRITE(0,*)'Output written successfully to standard output unit.'
      CALL EXIT(0)

!*    7. ERRORS WHILE READING INPUT FILES BRANCHES TO HERE TO TERMINATE
!        --------------------------------------------------------------
  900 WRITE(0,*)'Unknown input ERROR while reading file: '//
     &                                infile(I)(1:LEN_TRIM(infile(I)))
      CALL EXIT(7)
      END

