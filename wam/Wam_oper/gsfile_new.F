C ----------------------------------------------------------------------
C
      SUBROUTINE GSFILE (IU06, IUNIT, ICH, CDATE, CDATEF, FILEID, OPT)

C ----------------------------------------------------------------------
C
C**** *GSFILE* - GETS/SAVES FILES FROM/TO MASS STORAGE.
C
C     H. GUNTHER    GKSS/ECMWF    OCTOBER 1989
C     P. JANSSEN    KNMI          OCTOBER 1990   YMP-MODIFICATION
C     H. GUNTHER    GKSS/ECMWF    OCTOBER 1990   NEW FILE NAMES.
C     J. BIDLOT     ECMWF         MAY 1996       MIGRATION TO FUJITSU 
C     J. BIDLOT     ECMWF         JUNE 1996      PARALLEL RUN 
C     J. BIDLOT     ECMWF         FEBRUARY 1997  USE OF ECFS ON FUJITSU
C
C*    PURPOSE.
C     --------
C
C       GETS OR SAVES FILES FROM / TO MASS STORAGE.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *GSFILE (IU06, IUNIT, CDATE, CDATEF, FILEID, OPT)*
C
C         *IU06*   INTEGER    UNIT FOR PRINTER MESSAGES.
C         *IUNIT*  INTEGER    FORTRAN UNIT.
C         *ICH*    INTEGER    IF ICH = 0 THEN FILE WITH DTG IS SAVED,
C                             ELSE FILE FROM A LIST IS SAVED.
C         *CDATE*  CHAR*10    DATE (YYMMDDHHMM) OF FILE TO BE SAVED.
C         *CDATEF* CHAR*10    DATE TIME GROUP OF FORECAST START.
C         *FILEID* CHARACTER  FILE ID.
C         *OPT*    CHARACTER  OPTION
C                             = 'S' SAVE FILE CONNECTED TO IUNIT
C                             = 'G' GET FILE AND CONNECT FILE TO IUNIT
C
C     EXTERNALS.
C     ----------
C
C       *ABORT1*     - TERMINATES PROCESSING.
C       *DIFDATE*   - COMPUTE TIME DIFFERENCES.
C       *IECF_LEN*  -  LENGTH OF A CHARATER VARIABLE.
C
C     METHOD.
C     -------
C
C        THE METHOD USED IN THIS SUB DEPENDS UPON THE COMPUTER
C        ENVIROMENT. IN ITS PRESENT FORM THE ROUTINE IS ADOPTED TO
C        THE ECMWF FUJITSU SYSTEM AND USES THE MASS STORAGE ECFS.
C
C        IF ICH = 0 THEN THE
C        FILE NAME IS BUILT FROM THE 3 CHARACTER FILEID FOLLOWED
C        BY YYMMDDHHFFFF WHERE YYMMDDHH IS THE LAST ANALYSIS TIME AND
C        FFFF IS THE FORCECAST PERIOD IN HOURS.
C        THESE NUMBERS ARE CONSTRUCTED FROM CDATE AND CDATEF.
C
C        A POSITIVE ICH INDICATES THAT THE FILE NAME HAS TO BE
C        TAKE OUT OF ARRAY NAME(ICH)
C
C !!!!
C        THE PATH NAME FOR ECFS IS
C            PATH/FILENAME
C !!!!
C
C        IN PARTICULAR THESE SUB DOES THE FOLLOWING STEPS:
C               - CLOSES THE UNIT.
C               - CONSTRUCTS THE FILE NAME.
C               - CONSTRUCTS THE UNIT NAME.
C               - COPIES !!!!!!!! THE FILE.
C
C     REFERENCES.
C      -----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
#include "comtext.h"
C
#include "txttext.h"
C
      INTEGER IU06, IUNIT
      CHARACTER FILEID*3, OPT*1

      INTEGER QQNAME
      PARAMETER (QQNAME =  7)
      CHARACTER NAME(QQNAME)*16,
     1          FILENA*16, IU*8, PLIST*110,
     2          CDATE*10, CDATEF*10, CDATEH*8, FILENAME*72,
     3          MODE*1
C
C ----------------------------------------------------------------------
C
C*    1. DEFINE RESTART FILENAMES AND OPTIONS.
C        -------------------------------------
C
C*    1.1 FILE NAMES FOR RESTART FIELDS FROM AN ANALYSIS RUN.
C         ---------------------------------------------------
C
      DATA NAME( 1) / 'blspanal' / ,
     1     NAME( 2) / 'slatanal' / ,
     2     NAME( 3) / 'lawianal' / 
C
C*    1.2 FILE NAMES FOR RESTART FIELDS FROM A FORECAST RUN.
C         --------------------------------------------------
C
      DATA NAME( 4) / 'blspforc' / ,
     1     NAME( 5) / 'slatforc' / ,
     2     NAME( 6) / 'lawiforc' /
C
C*    1.3 MANAGEMENT FILE.
C         ----------------
C
      DATA NAME( 7) / 'waminfo' / 
C
C ----------------------------------------------------------------------
C
C*    2. CLOSE UNIT.
C        -----------
C
 2000 CONTINUE
      CLOSE (UNIT=IUNIT, STATUS='KEEP')
C
C ----------------------------------------------------------------------
C
C*    3.0  CONSTRUCT FILE NAME.
C          --------------------
C
 3000 CONTINUE
      INDEX = ICH
      IF (ICH.EQ.0) THEN
         FILENA = FILEID
         IF (CDATE.LE.CDATEF) THEN
            CDATEH = CDATE(1:8)
            ISHIFT = 0
         ELSE
            CDATEH = CDATEF(1:8)
            CALL DIFDATE (CDATEF, CDATE, ISHIFT)
            ISHIFT = ISHIFT/3600
         ENDIF
         FILENA(4:11) = CDATEH
         WRITE (FILENA(12:15),'(I4.4)') ISHIFT
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
      LNAME = IECF_LEN(FILENA)
C
C ----------------------------------------------------------------------
C
C*    4.0 CONSTRUCT UNIT NAME.
C         --------------------
C
 4000 CONTINUE
      IF ((IUNIT.GT.0) .AND. (IUNIT.LT.10)) THEN
         WRITE(IU,'(A,I1)') 'fort.',IUNIT
      ELSE IF ((IUNIT.GT.9) .AND. (IUNIT.LT.100)) THEN
         WRITE(IU,'(A,I2)') 'fort.',IUNIT
      ELSE
         WRITE(IU,'(A8)') IUNIT
      ENDIF
      LIU   = IECF_LEN(IU)
C
C ----------------------------------------------------------------------
C
C*    5. CONSTRUCTED ECFS PARAMETER LIST.
C        ----------------------------------
C
 5000 CONTINUE
      LIP   = IECF_LEN(PATH)
      LUSI = IECF_LEN(USERID)
      IF (OPT.EQ.'S' .AND. LIP .GT.0 ) THEN
C
C     5.1 SAVE FILE.
C         ----------
C
         IF (INDEX.EQ.7) THEN
            IRET=ISHELL('chmod go+rw '// IU(1:LIU))
C           CALL SYSTEM('chmod go+rw '// IU(1:LIU))
         ELSE
            IRET=ISHELL('chmod go+r '// IU(1:LIU))
C           CALL SYSTEM('chmod go+r '// IU(1:LIU))
         ENDIF
         PLIST = '-o '//IU(1:LIU)//' ec:'
     1      //PATH(1:LIP)//'/'//FILENA(1:LNAME)
C
C
C*    5.2 GET FILE.
C         ---------
C
      ELSEIF (OPT.EQ.'G' .AND. LIP .GT. 0) THEN
         PLIST = '-n ec:'//PATH(1:LIP)//'/'//FILENA(1:LNAME)
     1      //' '//IU(1:LIU)
C
C*    5.3 ERROR IN OPTION.
C         ----------------
C
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
     .                 ' MODEL RESTART --> USE FILE_TRANSFER TO COPY '//
     .                 ' RESTART FILES' 
      ELSEIF ( LIP .EQ. 0 ) THEN
        WRITE(IU06,*) '  '
        WRITE(IU06,*) ' * W A R N I N G in --GSFILE-- '//
     .                 ' NO PATH --> NO CALL TO ECFS !! '
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    6.0 EXECUTE ECFILE COMMAND.
C         -----------------------
C
 6000 CONTINUE
      IF ( LIP .GT. 0 ) THEN
        WRITE(IU06,*) '  '
        WRITE(IU06,*) ' --Ecp COMMAND IN GSFILE--'
        WRITE(IU06,*) '  ',PLIST
 
        IFAIL = 0
        IRET=ISHELL ('ksh ecp '//PLIST)
C       CALL SYSTEM ('ksh ecp '//PLIST)
C
C ----------------------------------------------------------------------
C
C*    7.0 CHECK ECFILE RETURN CODE.
C         -------------------------
C
 7000   CONTINUE
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
C
C ----------------------------------------------------------------------
C
C*    8.0 NO PATH PROVIDED ONLY SAVES ARE RECOGNISED AND A CALL TO
C         FILE_TRANSFER IS MADE TO COPY THE FILE TO ITS NEW NAME AND
C         LOCATION.
C         -----------------------------------------------------------
C
 8000   CONTINUE

        CALL FILE_TRANSFER(IU06,IUNIT,FILENA,OPT)

      ELSEIF ( LIP .EQ. 0 .AND. OPT .EQ. 'G' ) THEN
        IF (LRESTARTED) THEN
          CALL FILE_TRANSFER(IU06,IUNIT,FILENA,OPT)
        ELSE
          WRITE(IU06,*) ' FILE '//IU(1:LIU)//' EXPECTED TO BE ONLINE'
        ENDIF

      ENDIF

      RETURN
      END
