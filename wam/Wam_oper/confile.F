C=======================================================================
C
      SUBROUTINE CONFILE (kuso, kunit, cdate, yafileid, kfail)
C
C-----------------------------------------------------------------------
C
C****  *confile*  *data aquisition*
C
C     b. hansen           ECMWF       july    1991
C
C     PURPOSE.
C     --------
C     Opens data file unformatted for the use of a program.
C     handling of files following the naming convention
C     XXXYYMMDDHHMM
C     where XXX is the yafileid
C           YY  is the year
C           MM  is the month
C           DD  is the day
C           HH  is the hour
C           MM  is the minute
C
C**   INTERFACE.
C     ----------
C        CALL  confile (kuso, kunit, cdate, yafileid, kfail)
C
C        *kuso*       INTEGER       FORTRAN UNIT FOR THE PRINTER OUTPUT.
C        *kunit*      INTEGER       FORTRAN UNIT FOR THE REQUESTED FILE.
C        *cdate*      CHAR*10       DATE TIME GROUP (YYMMDDHHMM).
C        *yafileid*   CHARACTER*3   FILE ID.
C        *kfail*      INTEGER       = 0 NO ERROR.
C                                   > 0 ERROR.
C
C     EXTERNALS.
C     ----------
C        NONE.
C
C     METHOD.
C     -------
C        "confile" checks if the required file resides in the working catalog.
C        If not "confile" will return with kfail > 0.
C        If the requested file is not already opened this will be done
C        on unit kunit. IF the unit is already in use it will be released.
C
C     REFERENCES.
C     -----------
C        NONE
C
C-------------------------------------------------------------------------------
C
C* 0.    DEFINITIONS.
C  ------------------
C
      INTEGER kunit
      CHARACTER* 3 yafileid
C
      CHARACTER*130 younitnam
      CHARACTER* 13 yoname
      CHARACTER*10 cdate
      LOGICAL   ex, od
C
C-------------------------------------------------------------------------------
C
      LOGICAL lotest, losuvi
      CHARACTER*7 csubna
      DATA csubna /"confile"/
C
C-------------------------------------------------------------------------------
C
C* 1.0   INITIALIZATION.
C  ---------------------
C
 1000 CONTINUE
      kfail = 0
      losuvi = .TRUE.
      lotest = .FALSE.
C
C-------------------------------------------------------------------------------
C
C* 2.    CONSTRUCT FILE NAME.
C  --------------------------
C
      IF (losuvi) THEN
CCC        WRITE(kuso,'(1h1,a2)')   ' '
        WRITE(kuso,*) ' FILE INQUIRY FOR ',yafileid,cdate
      ENDIF
      yoname = yafileid
      yoname(4:13) = cdate
C
C-------------------------------------------------------------------------------
C
C* 3.0   CHECK STATUS OF UNIT.
C  ---------------------------
C
      INQUIRE ( UNIT=kunit, OPENED=od, NAME=younitnam, ERR=8200)
      IF (od) THEN
        IF (iecf_len(younitnam) .EQ. 13 ) THEN
          IF ( younitnam(1:13) .EQ. yoname ) THEN
            REWIND kunit
            RETURN
          ENDIF
        ELSE
          CLOSE (UNIT=kunit, ERR=8300)
          IF (losuvi) THEN
            WRITE(kuso,*) ' UNIT ',kunit,' WAS OPENED '
            WRITE(kuso,*) ' AND CONNECTED WITH FILE ',younitnam(1:50)
            WRITE(kuso,*) ' IS NOW CLOSED '
          ENDIF
        ENDIF
      ENDIF
C
C-------------------------------------------------------------------------------
C
C* 4.0   CHECK STATUS OF FILE.
C  ---------------------------
C
      INQUIRE (FILE=yoname, EXIST=ex, OPENED=od, NUMBER=knum, ERR =8100)
      IF (.NOT.ex) THEN
         kfail = 1
         RETURN
      ELSE
        OPEN (FILE=yoname, UNIT=kunit, FORM='unformatted', ERR=8400)
        IF (losuvi) THEN
          WRITE(kuso,*) ' FILE ',yoname,
     .                  ' IS NOW OPENED ON UNIT ', kunit
          WRITE(kuso,*) '  '
        ENDIF
      ENDIF
      RETURN
C
C-------------------------------------------------------------------------------
C
C* 8.0   RETURN IN CASE OF FAILURE.
C  --------------------------------
C
 8100 CONTINUE
      kfail=20
      IF (losuvi) THEN
        WRITE(kuso,*) ' INQUIRY ON FILE NAMED ',yoname,' FAILED'
        WRITE(kuso,*) '  '
      ENDIF
      RETURN
 8200 CONTINUE
      kfail=30
      IF (losuvi) THEN
        WRITE(kuso,*) ' INQUIRY ON UNIT ',kunit,' FAILED'
        WRITE(kuso,*) '  '
      ENDIF
      RETURN
 8300 CONTINUE
      kfail=40
      IF (losuvi) THEN
        WRITE(kuso,*) ' CLOSE ON UNIT ',kunit,' FAILED'
        WRITE(kuso,*) '  '
      ENDIF
      RETURN
 8400 CONTINUE
      kfail=50
      IF (losuvi) THEN
        WRITE(kuso,*) ' OPEN ON FILE NAMED ',yoname,' FAILED'
        WRITE(kuso,*) '  '
      ENDIF
      RETURN
      END
