C ======================================================================
C
      SUBROUTINE DIFDATES (CDATE1, CDATE2, KSHIFT)
C
C ======================================================================
C
C****     *DIFDATES* - TO COMPUTE TIME DIFFERENCE.
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
C         *CALL* *DIFDATES (CDATE1, CDATE2, KSHIFT)*
C
C     I/   *CDATE1* CHAR*10 - DATE TIME GROUP (YYMMDDHHMM) OR
C                   CHAR*12 - DATE TIME GROUP (YYYYMMDDHHMM) OR
C                   CHAR*14 - DATE TIME GROUP (YYYYMMDDHHMMSS).
C     I/   *CDATE2* CHAR*10 - DATE TIME GROUP (YYMMDDHHMM) OR
C                   CHAR*12 - DATE TIME GROUP (YYYYMMDDHHMM) OR
C                   CHAR*14 - DATE TIME GROUP (YYYYMMDDHHMMSS).
C      /O  *KSHIFT* INTEGER - DIFFERENCE IN SECONDS (CDATE2-CDATE1).
C
C     METHOD.
C     -------
C
C       THE DATE AND TIME CAN BE SUPPLIED IN THE OLD 10 CHARACTER
C       WAM FORMAT OR IN THE Y2K COMPLIANT 14 CHARACTER FORMAT.
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
C      MODIFIKATIONS
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
     .     '(" NON-Y2K COMPLIANT  DATE IN USE ", A', IL, ' )'
        WRITE(6, CLFMT ) CDATE1
        LLND = .FALSE.
        READ(CDATE1,'(5i2)')IYEAR1,IMON1,IDAY1,IHOUR1,IMIN1
        IYEAR1 = IYEAR1 + 1900
        ISEC1 = 0
      ELSEIF (IL==12) THEN
        LLND = .TRUE.
        READ(CDATE1,'(i4, 5i2)')IYEAR1,IMON1,IDAY1,IHOUR1,IMIN1
        ISEC1 = 0
      ELSEIF (IL==14) THEN
        LLND = .TRUE.
        READ(CDATE1,'(i4, 5i2)')IYEAR1,IMON1,IDAY1,IHOUR1,IMIN1,ISEC1
      ELSE
        WRITE(CLFMT,'(A58, I2.2, A2)')
     .  '(" DIFDATES:  FATAL@! DATE IS ",I2," CHARACTERS LONG!!= ",A',
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
     .     '(" NON-Y2K COMPLIANT  DATE IN USE ", A', IL, ' )'
        WRITE(6, CLFMT ) CDATE2
        LLND = .FALSE.
        READ(CDATE2,'(5i2)')IYEAR2,IMON2,IDAY2,IHOUR2,IMIN2
        IYEAR2 = IYEAR2 + 1900
        ISEC2 = 0
      ELSEIF (IL==12) THEN
        LLND = .TRUE.
        READ(CDATE2,'(i4, 5i2)')IYEAR2,IMON2,IDAY2,IHOUR2,IMIN2
        ISEC2 = 0
      ELSEIF (IL==14) THEN
        LLND = .TRUE.
        READ(CDATE2,'(i4, 5i2)')IYEAR2,IMON2,IDAY2,IHOUR2,IMIN2,ISEC2
      ELSE
        WRITE(CLFMT,'(A58, I2.2, A2)')
     .  '(" DIFDATES:  FATAL@! DATE IS ",I2," CHARACTERS LONG!!= ",A',
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
        WRITE(6, '(" DIFFDATES:  FATAL@! IRET= ", I3 )') IRET
        CALL ABORT
      ENDIF
C
      RETURN
      END
