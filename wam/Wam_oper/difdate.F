C ======================================================================
C
      SUBROUTINE DIFDATE (CDATE1, CDATE2, KSHIFT)
C
C ======================================================================
C
C****     *DIFDATE* - TO COMPUTE TIME DIFFERENCE.
C
C     B. HANSEN           ECMWF       JULY    1998
C
C
C**   PURPOSE.
C     --------
C
C          COMPUTE THE SECONDS BETWEEN THE INPUT DATES.
C          DATES HAVE TO BE IN CONSECUTIVE YEARS.
C
C**   INTERFACE.
C     ----------
C
C         *CALL* *DIFDATE (CDATE1, CDATE2, KSHIFT)*
C
C     I/   *CDATE1* CHAR*12 - DATE TIME GROUP (YYYYMMDDHHMM) OR
C                   CHAR*14 - DATE TIME GROUP (YYYYMMDDHHMMSS).
C     I/   *CDATE2* CHAR*12 - DATE TIME GROUP (YYYYMMDDHHMM) OR
C                   CHAR*14 - DATE TIME GROUP (YYYYMMDDHHMMSS).
C      /O  *KSHIFT* INTEGER - DIFFERENCE IN SECONDS (CDATE2-CDATE1).
C
C     METHOD.
C     -------
C
C       THE DATE AND TIME CAN BE SUPPLIED IN THE Y2K COMPLIANT
C       12 or 14 CHARACTER FORMAT. IF THE 10 CHARACTER FORMAT IS
C       USED AN ERROR MESSAGE WILL BE ISSUED. 
C       THE CHARACTER STRING IS BROKEN INTO INTEGERS TO OBTAIN
C       THE INPUT ARGUMENT FOR SECDIFF. THE DIFFERENCE IN SECONDES BETWEEN
C       THE TWO INPUT DATES IS RETURNED BACK AS KSHIFT.
C
C
C     EXTERNALS.
C     ----------
C
C       SECINCR (ECLIB).
C
C
C     REFERENCES.
C     -----------
C
C       OD MEMO (6 July 1998) File: O/PRO/90 FROM GROUP 2000.
C
C
C      MODIFICATIONS
C      -----------
C
C       NONE.
C
C
C
C ---------------------------------------------------------------------
C
      CHARACTER*(*) CDATE1, CDATE2
      INTEGER       KSHIFT
C
      CHARACTER* 80 CLFMT
      LOGICAL LLND
C
C ----------------------------------------------------------------------
C
C*    1.0 SPLITE DATE TIME GROUPS INTO SECOND, MINUTE, HOUR, DAY, MONTH,
C*        YEAR DEPENDENT ON THE FORMAT OF CDATE.
C         -----------------------------------------------------------
C
C*        SPLIT CDATE1.
C
      IL = LEN_TRIM(CDATE1)
      IF (IL==10) THEN
        WRITE(CLFMT,'(A39, I2.2, A2)')
     .     '(" DIFDATE: NON-Y2K COMPLIANT  DATE IN USE ", A', IL, ' )'
        WRITE(6, CLFMT ) CDATE1
        WRITE(*, CLFMT ) CDATE1
        LLND = .FALSE.
        CALL ABORT
ccc        READ(CDATE1,'(5i2)')IYEAR1,IMON1,IDAY1,IHOUR1,IMIN1
ccc        IYEAR1 = IYEAR1 + 1900
ccc        ISEC1 = 0
      ELSEIF (IL==12) THEN
        LLND = .TRUE.
        READ(CDATE1,'(i4, 5i2)')IYEAR1,IMON1,IDAY1,IHOUR1,IMIN1
        ISEC1 = 0
      ELSEIF (IL==14) THEN
        LLND = .TRUE.
        READ(CDATE1,'(i4, 5i2)')IYEAR1,IMON1,IDAY1,IHOUR1,IMIN1,ISEC1
C       with satellite data it might happen that the seconds were rounded
C       to 60. One needs to increment the date properly.
        IF(ISEC1==60) THEN
           IYEAR=IYEAR1
           IMON=IMON1
           IDAY=IDAY1
           IHOUR=IHOUR1
           IMIN=IMIN1
           ISEC=ISEC1-1
           ISHIFT=1
           CALL SECINCR  (IYEAR, IMON, IDAY, IHOUR, IMIN, ISEC, ISHIFT,
     .                IYEAR1, IMON1, IDAY1, IHOUR1, IMIN1, ISEC1, IRET)
           IF (IRET/=0) THEN
             WRITE(6, '(" DIFDATE:  ERROR CALLING SECINCR FOR ISEC1" )')
             WRITE(6, '(" DIFDATE:  FATAL@! IRET= ", I3 )') IRET
             WRITE(6, '("           IYEAR ", I4 )') IYEAR
             WRITE(6, '("           IMON ", I2 )') IMON
             WRITE(6, '("           IDAY ", I2 )') IDAY
             WRITE(6, '("           IHOUR ", I2 )') IHOUR
             WRITE(6, '("           IMIN ", I2 )') IMIN
             WRITE(6, '("           ISEC ", I2 )') ISEC
             WRITE(6, '("           ISHIFT ", I12 )') ISHIFT
             CALL ABORT
           ENDIF
        ENDIF
      ELSE
        WRITE(CLFMT,'(A58, I2.2, A2)')
     .  '(" DIFDATE:  FATAL@! DATE IS ",I2," CHARACTERS LONG!!= ",A',
     .    IL, ' )'
        WRITE(6, CLFMT) IL, CDATE1
        STOP
      ENDIF
C
C*        SPLIT CDATE2.
C
      IL = LEN_TRIM(CDATE2)
      IF (IL==10) THEN
        WRITE(CLFMT,'(A39, I2.2, A2)')
     .     '(" DIFDATE: NON-Y2K COMPLIANT  DATE IN USE ", A', IL, ' )'
        WRITE(6, CLFMT ) CDATE2
        WRITE(*, CLFMT ) CDATE2
        LLND = .FALSE.
        CALL ABORT
ccc        READ(CDATE2,'(5i2)')IYEAR2,IMON2,IDAY2,IHOUR2,IMIN2
ccc        IYEAR2 = IYEAR2 + 1900
ccc        ISEC2 = 0
      ELSEIF (IL==12) THEN
        LLND = .TRUE.
        READ(CDATE2,'(i4, 5i2)')IYEAR2,IMON2,IDAY2,IHOUR2,IMIN2
        ISEC2 = 0
      ELSEIF (IL==14) THEN
        LLND = .TRUE.
        READ(CDATE2,'(i4, 5i2)')IYEAR2,IMON2,IDAY2,IHOUR2,IMIN2,ISEC2
C       with satellite data it might happen that the seconds were rounded
C       to 60. One needs to increment the date properly.
        IF(ISEC2==60) THEN
           IYEAR=IYEAR2
           IMON=IMON2
           IDAY=IDAY2
           IHOUR=IHOUR2
           IMIN=IMIN2
           ISEC=ISEC2-1
           ISHIFT=1
           CALL SECINCR  (IYEAR, IMON, IDAY, IHOUR, IMIN, ISEC, ISHIFT,
     .                IYEAR2, IMON2, IDAY2, IHOUR2, IMIN2, ISEC2, IRET)
           IF (IRET/=0) THEN
             WRITE(6, '(" DIFDATE:  ERROR CALLING SECINCR FOR ISEC2" )')
             WRITE(6, '(" DIFDATE:  FATAL@! IRET= ", I3 )') IRET
             WRITE(6, '("           IYEAR ", I4 )') IYEAR
             WRITE(6, '("           IMON ", I2 )') IMON
             WRITE(6, '("           IDAY ", I2 )') IDAY
             WRITE(6, '("           IHOUR ", I2 )') IHOUR
             WRITE(6, '("           IMIN ", I2 )') IMIN
             WRITE(6, '("           ISEC ", I2 )') ISEC
             WRITE(6, '("           ISHIFT ", I12 )') ISHIFT
             CALL ABORT
           ENDIF
        ENDIF
      ELSE
        WRITE(CLFMT,'(A58, I2.2, A2)')
     .  '(" DIFDATE:  FATAL@! DATE IS ",I2," CHARACTERS LONG!!= ",A',
     .    IL, ' )'
        WRITE(6, CLFMT) IL, CDATE2
        STOP
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    2.0 DETERMINE DIFFERENCE IN SECONDES.
C         ---------------------------------
C
      CALL SECDIFF  (IYEAR2,IMON2,IDAY2,IHOUR2,IMIN2,ISEC2,
     .               IYEAR1,IMON1,IDAY1,IHOUR1,IMIN1,ISEC1,
     .               KSHIFT,IRET)

      IF (IRET/=0) THEN
        WRITE(6, '(" DIFDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(6, '("           IYEAR1 ", I4 )') IYEAR1 
        WRITE(6, '("           IMON1 ", I2 )') IMON1 
        WRITE(6, '("           IDAY1 ", I2 )') IDAY1 
        WRITE(6, '("           IHOUR1 ", I2 )') IHOUR1 
        WRITE(6, '("           IMIN1 ", I2 )') IMIN1 
        WRITE(6, '("           ISEC1 ", I2 )') ISEC1 
        WRITE(6, '("           IYEAR2 ", I4 )') IYEAR2 
        WRITE(6, '("           IMON2 ", I2 )') IMON2 
        WRITE(6, '("           IDAY2 ", I2 )') IDAY2 
        WRITE(6, '("           IHOUR2 ", I2 )') IHOUR2 
        WRITE(6, '("           IMIN2 ", I2 )') IMIN2 
        WRITE(6, '("           ISEC2 ", I2 )') ISEC2 
        WRITE(6, '("           KSHIFT ", I12 )') KSHIFT 
        CALL ABORT
      ENDIF
C
      RETURN
      END
