! ======================================================================

      SUBROUTINE DIFDATE (CDATE1, CDATE2, KSHIFT)

! ======================================================================

!****     *DIFDATE* - TO COMPUTE TIME DIFFERENCE.

!     B. HANSEN           ECMWF       JULY    1998


!**   PURPOSE.
!     --------

!          COMPUTE THE SECONDS BETWEEN THE INPUT DATES.
!          DATES HAVE TO BE IN CONSECUTIVE YEARS.

!**   INTERFACE.
!     ----------

!         *CALL* *DIFDATE (CDATE1, CDATE2, KSHIFT)*

!     I/   *CDATE1* CHAR*12 - DATE TIME GROUP (YYYYMMDDHHMM) OR
!                   CHAR*14 - DATE TIME GROUP (YYYYMMDDHHMMSS).
!     I/   *CDATE2* CHAR*12 - DATE TIME GROUP (YYYYMMDDHHMM) OR
!                   CHAR*14 - DATE TIME GROUP (YYYYMMDDHHMMSS).
!      /O  *KSHIFT* INTEGER - DIFFERENCE IN SECONDS (CDATE2-CDATE1).

!     METHOD.
!     -------

!       THE DATE AND TIME CAN BE SUPPLIED IN THE Y2K COMPLIANT
!       12 or 14 CHARACTER FORMAT. IF THE 10 CHARACTER FORMAT IS
!       USED AN ERROR MESSAGE WILL BE ISSUED. 
!       THE CHARACTER STRING IS BROKEN INTO INTEGERS TO OBTAIN
!       THE INPUT ARGUMENT FOR SECDIFF. THE DIFFERENCE IN SECONDES
!       BETWEEN THE TWO INPUT DATES IS RETURNED BACK AS KSHIFT.


!     EXTERNALS.
!     ----------

!       SECINCR (ECLIB).


!     REFERENCES.
!     -----------

!       OD MEMO (6 July 1998) File: O/PRO/90 FROM GROUP 2000.


!      MODIFICATIONS
!      -----------

!       NONE.



! ---------------------------------------------------------------------

      CHARACTER*(*) CDATE1, CDATE2
      INTEGER       KSHIFT

      CHARACTER* 80 CLFMT
      LOGICAL LLND

! ----------------------------------------------------------------------

!*    1.0 SPLITE DATE TIME GROUPS INTO SECOND, MINUTE, HOUR, DAY, MONTH,
!*        YEAR DEPENDENT ON THE FORMAT OF CDATE.
!         -----------------------------------------------------------

!*        SPLIT CDATE1.

      IL = LEN_TRIM(CDATE1)
      IF (IL==10) THEN
        WRITE(CLFMT,'(A49, I2.2, A2)')
     &     '(" DIFDATE : NON-Y2K COMPLIANT DATE1 IN USE ", A', IL, ' )'
        WRITE(6, CLFMT ) CDATE1
        WRITE(*, CLFMT ) CDATE1
        LLND = .FALSE.
        CALL ABORT
      ELSEIF (IL==12) THEN
        LLND = .TRUE.
        READ(CDATE1,'(i4, 5i2)')IYEAR1,IMON1,IDAY1,IHOUR1,IMIN1
        ISEC1 = 0
      ELSEIF (IL==14) THEN
        LLND = .TRUE.
        READ(CDATE1,'(i4, 5i2)')IYEAR1,IMON1,IDAY1,IHOUR1,IMIN1,ISEC1
!       with satellite data it might happen that the seconds were
!       rounded to 60. One needs to increment the date properly.
        IF(ISEC1==60) THEN
          IYEAR=IYEAR1
          IMON=IMON1
          IDAY=IDAY1
          IHOUR=IHOUR1
          IMIN=IMIN1
          ISEC=ISEC1-1
          ISHIFT=1
          CALL SECINCR  (IYEAR, IMON, IDAY, IHOUR, IMIN, ISEC, ISHIFT,
     &     IYEAR1, IMON1, IDAY1, IHOUR1, IMIN1, ISEC1, IRET)
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
      ELSEIF (IL==0) THEN
        WRITE(6, *) ' '
        WRITE(6, *) ' DIFDATE:  FATAL@! DATE1 IS 0 CHARACTER LONG'
        WRITE(6, *) ' '
        CALL ABORT
      ELSE
        WRITE(6, *) ' '
        WRITE(CLFMT,'(A61, I2.2, A2)')
     &  '(" DIFDATE :  FATAL@! DATE1 IS ",I2," CHARACTERS LONG!!= ",A',
     &    IL, ' )'
        WRITE(6, CLFMT) IL, CDATE1
        WRITE(6, *) ' '
        CALL ABORT
      ENDIF

!*        SPLIT CDATE2.

      IL = LEN_TRIM(CDATE2)
      IF (IL==10) THEN
        WRITE(CLFMT,'(A49, I2.2, A2)')
     &     '(" DIFDATE : NON-Y2K COMPLIANT DATE2 IN USE ", A', IL, ' )'
        WRITE(6, CLFMT ) CDATE2
        WRITE(*, CLFMT ) CDATE2
        LLND = .FALSE.
        CALL ABORT
      ELSEIF (IL==12) THEN
        LLND = .TRUE.
        READ(CDATE2,'(i4, 5i2)')IYEAR2,IMON2,IDAY2,IHOUR2,IMIN2
        ISEC2 = 0
      ELSEIF (IL==14) THEN
        LLND = .TRUE.
        READ(CDATE2,'(i4, 5i2)')IYEAR2,IMON2,IDAY2,IHOUR2,IMIN2,ISEC2
!       with satellite data it might happen that the seconds were
!       rounded to 60. One needs to increment the date properly.
        IF(ISEC2==60) THEN
          IYEAR=IYEAR2
          IMON=IMON2
          IDAY=IDAY2
          IHOUR=IHOUR2
          IMIN=IMIN2
          ISEC=ISEC2-1
          ISHIFT=1
          CALL SECINCR  (IYEAR, IMON, IDAY, IHOUR, IMIN, ISEC, ISHIFT,
     &     IYEAR2, IMON2, IDAY2, IHOUR2, IMIN2, ISEC2, IRET)
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
      ELSEIF (IL==0) THEN
        WRITE(6, *) ' '
        WRITE(6, *) ' DIFDATE:  FATAL@! DATE2 IS 0 CHARACTER LONG'
        WRITE(6, *) ' '
        CALL ABORT
      ELSE
        WRITE(6, *) ' '
        WRITE(CLFMT,'(A61, I2.2, A2)')
     &  '(" DIFDATE :  FATAL@! DATE2 IS ",I2," CHARACTERS LONG!!= ",A',
     &    IL, ' )'
        WRITE(6, CLFMT) IL, CDATE2
        WRITE(6, *) ' '
        CALL ABORT
      ENDIF

! ----------------------------------------------------------------------

!*    2.0 DETERMINE DIFFERENCE IN SECONDES.
!         ---------------------------------

      CALL SECDIFF  (IYEAR2,IMON2,IDAY2,IHOUR2,IMIN2,ISEC2,
     &               IYEAR1,IMON1,IDAY1,IHOUR1,IMIN1,ISEC1,
     &               KSHIFT,IRET)

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

      RETURN
      END SUBROUTINE DIFDATE


