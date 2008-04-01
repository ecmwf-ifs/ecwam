! ----------------------------------------------------------------------

      SUBROUTINE GSFILE (IU06, IUNIT, ICH, CDATE, CDATEF, FILEID, OPT)

! ----------------------------------------------------------------------

!**** *GSFILE* - GETS/SAVES FILES FROM/TO MASS STORAGE.

!     H. GUNTHER    GKSS/ECMWF    OCTOBER 1989
!     P. JANSSEN    KNMI          OCTOBER 1990   YMP-MODIFICATION
!     H. GUNTHER    GKSS/ECMWF    OCTOBER 1990   NEW FILE NAMES.
!     J. BIDLOT     ECMWF         MAY 1996       MIGRATION TO FUJITSU 
!     J. BIDLOT     ECMWF         JUNE 1996      PARALLEL RUN 
!     J. BIDLOT     ECMWF         FEBRUARY 1997  USE OF ECFS ON FUJITSU

!*    PURPOSE.
!     --------

!       GETS OR SAVES FILES FROM / TO MASS STORAGE.

!**   INTERFACE.
!     ----------

!       *CALL* *GSFILE (IU06, IUNIT, CDATE, CDATEF, FILEID, OPT)*

!         *IU06*   INTEGER    UNIT FOR PRINTER MESSAGES.
!         *IUNIT*  INTEGER    FORTRAN UNIT.
!         *ICH*    INTEGER    IF ICH = 0 THEN FILE WITH DTG IS SAVED,
!                             ELSE FILE FROM A LIST IS SAVED.
!         *CDATE*  CHAR*14    DATE (YYYYMMDDHHMM) OF FILE TO BE SAVED.
!         *CDATEF* CHAR*14    DATE TIME GROUP OF FORECAST START.
!         *FILEID* CHARACTER  FILE ID.
!         *OPT*    CHARACTER  OPTION
!                             = 'S' SAVE FILE CONNECTED TO IUNIT
!                             = 'G' GET FILE AND CONNECT FILE TO IUNIT

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.
!       *DIFDATE*   - COMPUTE TIME DIFFERENCES.

!     METHOD.
!     -------

!        THE METHOD USED IN THIS SUB DEPENDS UPON THE COMPUTER
!        ENVIROMENT. IN ITS PRESENT FORM THE ROUTINE IS ADOPTED TO
!        THE ECMWF FUJITSU SYSTEM AND USES THE MASS STORAGE ECFS.

!        IF ICH = 0 THEN THE
!        FILE NAME IS BUILT FROM THE 3 CHARACTER FILEID FOLLOWED
!        BY YYMMDDHHFFFF WHERE YYMMDDHH IS THE LAST ANALYSIS TIME AND
!        FFFF IS THE FORCECAST PERIOD IN HOURS.
!        THESE NUMBERS ARE CONSTRUCTED FROM CDATE AND CDATEF.

!        A POSITIVE ICH INDICATES THAT THE FILE NAME HAS TO BE
!        TAKE OUT OF ARRAY NAME(ICH)

!        IN PARTICULAR THESE SUB DOES THE FOLLOWING STEPS:
!               - CLOSES THE UNIT.
!               - CONSTRUCTS THE FILE NAME.
!               - CONSTRUCTS THE UNIT NAME.
!               - COPIES !!!!!!!! THE FILE.

!     REFERENCES.
!      -----------

!       NONE.

! ----------------------------------------------------------------------


      USE YOWTEXT  , ONLY : LRESTARTED,USERID   ,PATH


      INTEGER IU06, IUNIT
      CHARACTER FILEID*3, OPT*1

      INTEGER QQNAME
      PARAMETER (QQNAME =  7)
      CHARACTER NAME(QQNAME)*17,
     &          FILENA*17, IU*8, PLIST*110,
     &          CDATE*14, CDATEF*14, CDATEH*10

! ----------------------------------------------------------------------

!*    1. DEFINE RESTART FILENAMES AND OPTIONS.
!        -------------------------------------

!*    1.1 FILE NAMES FOR RESTART FIELDS FROM AN ANALYSIS RUN.
!         ---------------------------------------------------

      DATA NAME( 1) / 'blspanal' / ,
     &     NAME( 2) / 'slatanal' / ,
     &     NAME( 3) / 'lawianal' / 

!*    1.2 FILE NAMES FOR RESTART FIELDS FROM A FORECAST RUN.
!         --------------------------------------------------

      DATA NAME( 4) / 'blspforc' / ,
     &     NAME( 5) / 'slatforc' / ,
     &     NAME( 6) / 'lawiforc' /

!*    1.3 MANAGEMENT FILE.
!         ----------------

      DATA NAME( 7) / 'waminfo' / 

! ----------------------------------------------------------------------

!*    2. CLOSE UNIT.
!        -----------

      CLOSE (UNIT=IUNIT, STATUS='KEEP')

! ----------------------------------------------------------------------

!*    3.0  CONSTRUCT FILE NAME.
!          --------------------

      INDEX = ICH
      IF (ICH.EQ.0) THEN
        FILENA = FILEID
        IF (CDATE.LE.CDATEF) THEN
          CDATEH = CDATE(1:10)
          ISHIFT = 0
        ELSE
          CDATEH = CDATEF(1:10)
          CALL DIFDATE (CDATEF, CDATE, ISHIFT)
          ISHIFT = ISHIFT/3600
        ENDIF
        FILENA(4:13) = CDATEH
        WRITE (FILENA(14:17),'(I4.4)') ISHIFT
      ELSE
        IF (INDEX.GT.QQNAME) THEN
          WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++'
          WRITE(IU06,*) ' +                                         +'
          WRITE(IU06,*) ' +      WARNING ERROR IN --GSFILE--        +'
          WRITE(IU06,*) ' +      ===========================        +'
          WRITE(IU06,*) ' + FILE NAME INDEX REQUESTED IS ',INDEX
          WRITE(IU06,*) ' + MAXIMUM INDEX ALLOWED IS QQNAME =',QQNAME
          WRITE(IU06,*) ' + NO GET OR SAVE PROCESSED                +'
          WRITE(IU06,*) ' + EXECUTION IS CONTINUED                  +'
          WRITE(IU06,*) ' +                                         +'
          WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++'
          RETURN
        ENDIF
        FILENA = NAME(INDEX)
      ENDIF
      LNAME = LEN_TRIM(FILENA)

! ----------------------------------------------------------------------

!*    4.0 CONSTRUCT UNIT NAME.
!         --------------------

      IF ((IUNIT.GT.0) .AND. (IUNIT.LT.10)) THEN
        WRITE(IU,'(A,I1)') 'fort.',IUNIT
      ELSEIF ((IUNIT.GT.9) .AND. (IUNIT.LT.100)) THEN
        WRITE(IU,'(A,I2)') 'fort.',IUNIT
      ELSE
        WRITE(IU,'(A8)') IUNIT
      ENDIF
      LIU   = LEN_TRIM(IU)

! ----------------------------------------------------------------------

!*    5. CONSTRUCTED ECFS PARAMETER LIST.
!        ----------------------------------

      LIP   = LEN_TRIM(PATH)
      LUSI = LEN_TRIM(USERID)
      IF (OPT.EQ.'S' .AND. LIP .GT.0 ) THEN

!     5.1 SAVE FILE.
!         ----------

        IF (INDEX.EQ.7) THEN
          CALL SYSTEM('chmod go+rw '// IU(1:LIU))
        ELSE
          CALL SYSTEM('chmod go+r '// IU(1:LIU))
        ENDIF
        PLIST = '-o '//IU(1:LIU)//' ec:'
     &   //PATH(1:LIP)//'/'//FILENA(1:LNAME)


!*    5.2 GET FILE.
!         ---------

      ELSEIF (OPT.EQ.'G' .AND. LIP .GT. 0) THEN
        PLIST = '-n ec:'//PATH(1:LIP)//'/'//FILENA(1:LNAME)
     &   //' '//IU(1:LIU)

!*    5.3 ERROR IN OPTION.
!         ----------------

      ELSEIF ((OPT .NE. 'S' .AND. OPT .NE. 'G').AND. LIP .GT. 0) THEN
        WRITE(IU06,*) ' *******************************************'
        WRITE(IU06,*) ' *                                         *'
        WRITE(IU06,*) ' *        FATAL ERROR IN --GSFILE--        *'
        WRITE(IU06,*) ' *        =========================        *'
        WRITE(IU06,*) ' * OPTION REQUESTED IS OPT= ',OPT
        WRITE(IU06,*) ' * ONLY OPTIONS S (SAVE) AND G (GET) EXIST *'
        WRITE(IU06,*) ' * PROCESSING WILL BE ABORTED              *'
        WRITE(IU06,*) ' *                                         *'
        WRITE(IU06,*) ' *******************************************'
        CALL ABORT1
      ELSEIF ( OPT.EQ.'G' .AND. LRESTARTED ) THEN
        WRITE(IU06,*) '  '
        WRITE(IU06,*) ' * W A R N I N G in --GSFILE-- '//
     &                 ' MODEL RESTART --> USE FILE_TRANSFER TO COPY '//
     &                 ' RESTART FILES' 
      ENDIF

! ----------------------------------------------------------------------

!*    6.0 EXECUTE ECFILE COMMAND.
!         -----------------------

      IF ( LIP .GT. 0 ) THEN
        WRITE(IU06,*) '  '
        WRITE(IU06,*) ' --Ecp COMMAND IN GSFILE--'
        WRITE(IU06,*) '  ',PLIST
 
        IFAIL = 0
        CALL SYSTEM ('ksh ecp '//PLIST)

! ----------------------------------------------------------------------

!*    7.0 CHECK ECFILE RETURN CODE.
!         -------------------------

        IF(IFAIL.NE.0) THEN
          WRITE(IU06,*) ' *******************************************'
          WRITE(IU06,*) ' *                                         *'
          WRITE(IU06,*) ' *        FATAL ERROR IN --GSFILE--        *'
          WRITE(IU06,*) ' *        =========================        *'
          WRITE(IU06,*) ' * ECFS ERROR                              *'
          WRITE(IU06,*) ' * IFAIL = ',IFAIL
          WRITE(IU06,*) ' * PROCESSING WILL BE ABORTED              *'
          WRITE(IU06,*) ' *                                         *'
          WRITE(IU06,*) ' *******************************************'
          CALL ABORT1
        ENDIF
      ELSEIF ( LIP .EQ. 0 .AND. OPT .EQ. 'S' ) THEN

! ----------------------------------------------------------------------

!*    8.0 NO PATH PROVIDED ONLY SAVES ARE RECOGNISED AND A CALL TO
!         FILE_TRANSFER IS MADE TO COPY THE FILE TO ITS NEW NAME AND
!         LOCATION.
!         -----------------------------------------------------------


        CALL FILE_TRANSFER(IU06,IUNIT,FILENA,OPT)

      ELSEIF ( LIP .EQ. 0 .AND. OPT .EQ. 'G' ) THEN
        IF (LRESTARTED) THEN
          CALL FILE_TRANSFER(IU06,IUNIT,FILENA,OPT)
        ELSE
          WRITE(IU06,*) ' FILE '//IU(1:LIU)//' EXPECTED TO BE ONLINE'
        ENDIF

      ENDIF

      RETURN
      END SUBROUTINE GSFILE
