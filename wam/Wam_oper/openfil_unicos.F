      SUBROUTINE OPENFIL (IU06, IUNIT, CDATE, FILEID, STAT,
     1                    USERID, RUNID, PATH, IFAIL)

C -------------------------------------------------------------------
C
C**** *OPENFIL* - OPENS A FILE ( UNICOS VERSION ).
C
C     RENATE PORTZ     MPI HAMBURG   JULY 1990
C     HEINZ GUNTHER    GKSS/ECMWF    MARCH 1991  MODIFIED FOR CYCLE_4.
C
C*    PURPOSE.
C     --------
C
C       OPENS A FILE WITH VARIABLE FILE NAME.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *OPENFIL (IU06, IUNIT, CDATE, FILEID, STAT,
C                        USERID, RUNID, PATH, IFAIL)
C         *IU06*    -   PRINTER OUTPUT UNIT.
C         *IUNIT*   -   LOGCIAL UNIT FILE IS ASSIGNED TO.
C         *CDATE*   -   DATE IN FILE NAME.
C         *FILEID*  -   FILE IDENTIFIER.
C         *STAT*    -   FILE STATUS 'OLD' OR 'NEW'.
C         *USERID*  -   USER IDENTIFIER.
C         *RUNID*   -   RUN IDENTIFIER.
C         *PATH*    -   PATH NAME (WITHOUT FILE NAME).
C         *IFAIL*   -   ERROR FLAG = 0 FILE OPENED
C                                    = 1 FILE NOT OPENED.
C
C     EXTERNALS.
C     ----------
C
C       *ABORT*     - TERMINATES PROCESSING.
C       *DIFDATE*   - COMPUTE TIME DIFFERENCES.
C
C     METHOD.
C     -------
C
C       THE METHOD USED IN THIS SUB DEPENDS UPON THE COMPUTER
C       ENVIROMENT. IN ITS PRESENT FORM THE ROUTINE IS ADOPTED TO
C       THE UNICOS CRAY SYSTEM.
C
C       THE FILE NAME IS BUILD FORM THE 3 CHARACTER FILEID FOLLOWED
C       BY YYMMDDHHMM0 WHERE YYMMDDHHMM IS THE DATE OF THE LAST
C       OUTPUT IN THIS FILE.
C       THESE NUMBERS ARE CONSTRUCTED FROM CDATE.
C
C       THE FULL PATH NAME IS PATH/FILE NAME
C
C     REFERENCES.
C     -----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
      INTEGER IU06, IUNIT, IFAIL
      CHARACTER FILEID*3, STAT*3, USERID*3, RUNID*3, PATH*60, CDATE*12

      CHARACTER FILENA*15, IU*8, FNAME*128
C
C        *FILENA*  CHARACTER   GENERATED FILE NAME.
C        *IU*      CHARACTER   UNIT ALIAS FILE NAME.
C        *FNAME*   CHARACTER   FULL PATH NAME.
C
C ----------------------------------------------------------------------
C
C*    1. BUILD FILE NAME.
C        ----------------
C
 1000 CONTINUE
      FILENA = FILEID
      FILENA(4:15) = CDATE(1:12)
C
C ----------------------------------------------------------------------
C
C*    2. CONSTRUCTED FULL PATH.
C        ----------------------
C
 2000 CONTINUE
      ILEN  = LEN_TRIM(PATH)
      LNAME = LEN_TRIM(FILENA)
      FNAME = ' '
      FNAME = PATH(1:ILEN)//'/'//FILENA(1:LNAME)
C
C*    3. OPEN FILE AND RETURN.
C        ---------------------
C
 3000 CONTINUE
      ILEN = LEN_TRIM(FNAME)
      OPEN (UNIT=IUNIT, FILE=FNAME, STATUS=STAT, FORM='UNFORMATTED',
     1      ERR=4000)

      WRITE (IU06,*) ' FILE [ ', FNAME(1:ILEN),' ] OPENED',
     1               ' AND ASSIGNED TO UNIT ', IUNIT

      IFAIL = 0
      RETURN
C
C*    4. ERROR MESSAGE.
C        --------------
C
 4000 CONTINUE
      IFAIL = 1
      WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++'
      WRITE (IU06,*) ' +                                 +'
      WRITE (IU06,*) ' +    WARNING ERROR SUB. OPENFIL   +'
      WRITE (IU06,*) ' +    ==========================   +'
      WRITE (IU06,*) ' +                                 +'
      WRITE (IU06,*) ' + COULD NOT OPEN FILE             +'
      WRITE (IU06,*) ' + FULL PATH NAME IS : ', FNAME(1:ILEN)
      WRITE (IU06,*) ' +                                 +'
      WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++'
      RETURN

      END
