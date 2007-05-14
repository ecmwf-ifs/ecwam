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

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ---------------------------------------------------------------------

      CHARACTER*(*) CDATE1, CDATE2
      INTEGER       KSHIFT

      CHARACTER* 80 CLFMT
      LOGICAL LLND

      REAL ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('DIFDATE',0,ZHOOK_HANDLE)

! ----------------------------------------------------------------------

!*    1.0 SPLITE DATE TIME GROUPS INTO SECOND, MINUTE, HOUR, DAY, MONTH,
!*        YEAR DEPENDENT ON THE FORMAT OF CDATE.
!         -----------------------------------------------------------

!*        SPLIT CDATE1.

      IL = LEN_TRIM(CDATE1)
      IF (IL==10) THEN
        WRITE(CLFMT,'(A49, I2.2, A2)')
     &     '(" DIFDATE : NON-Y2K COMPLIANT DATE1 IN USE ", A', IL, ' )'
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
            WRITE(*, '(" DIFDATE:  ERROR CALLING SECINCR FOR ISEC1" )')
            WRITE(*, '(" DIFDATE:  FATAL@! IRET= ", I3 )') IRET
            WRITE(*, '("           IYEAR ", I4 )') IYEAR
            WRITE(*, '("           IMON ", I2 )') IMON
            WRITE(*, '("           IDAY ", I2 )') IDAY
            WRITE(*, '("           IHOUR ", I2 )') IHOUR
            WRITE(*, '("           IMIN ", I2 )') IMIN
            WRITE(*, '("           ISEC ", I2 )') ISEC
            WRITE(*, '("           ISHIFT ", I12 )') ISHIFT
            CALL ABORT
          ENDIF
        ENDIF
      ELSEIF (IL==0) THEN
        WRITE(*, *) ' '
        WRITE(*, *) ' DIFDATE:  FATAL@! DATE1 IS 0 CHARACTER LONG'
        WRITE(*, *) ' '
        CALL ABORT
      ELSE
        WRITE(*, *) ' '
        WRITE(CLFMT,'(A61, I2.2, A2)')
     &  '(" DIFDATE :  FATAL@! DATE1 IS ",I2," CHARACTERS LONG!!= ",A',
     &    IL, ' )'
        WRITE(*, CLFMT) IL, CDATE1
        WRITE(*, *) ' '
        CALL ABORT
      ENDIF

!*        SPLIT CDATE2.

      IL = LEN_TRIM(CDATE2)
      IF (IL==10) THEN
        WRITE(CLFMT,'(A49, I2.2, A2)')
     &     '(" DIFDATE : NON-Y2K COMPLIANT DATE2 IN USE ", A', IL, ' )'
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
            WRITE(*, '(" DIFDATE:  ERROR CALLING SECINCR FOR ISEC2" )')
            WRITE(*, '(" DIFDATE:  FATAL@! IRET= ", I3 )') IRET
            WRITE(*, '("           IYEAR ", I4 )') IYEAR
            WRITE(*, '("           IMON ", I2 )') IMON
            WRITE(*, '("           IDAY ", I2 )') IDAY
            WRITE(*, '("           IHOUR ", I2 )') IHOUR
            WRITE(*, '("           IMIN ", I2 )') IMIN
            WRITE(*, '("           ISEC ", I2 )') ISEC
            WRITE(*, '("           ISHIFT ", I12 )') ISHIFT
            CALL ABORT
          ENDIF
        ENDIF
      ELSEIF (IL==0) THEN
        WRITE(*, *) ' '
        WRITE(*, *) ' DIFDATE:  FATAL@! DATE2 IS 0 CHARACTER LONG'
        WRITE(*, *) ' '
        CALL ABORT
      ELSE
        WRITE(*, *) ' '
        WRITE(CLFMT,'(A61, I2.2, A2)')
     &  '(" DIFDATE :  FATAL@! DATE2 IS ",I2," CHARACTERS LONG!!= ",A',
     &    IL, ' )'
        WRITE(*, CLFMT) IL, CDATE2
        WRITE(*, *) ' '
        CALL ABORT
      ENDIF

! ----------------------------------------------------------------------

!*    2.0 DETERMINE DIFFERENCE IN SECONDES.
!         ---------------------------------

      CALL SECDIFF  (IYEAR2,IMON2,IDAY2,IHOUR2,IMIN2,ISEC2,
     &               IYEAR1,IMON1,IDAY1,IHOUR1,IMIN1,ISEC1,
     &               KSHIFT,IRET)

      IF (IRET/=0) THEN
        WRITE(*, '(" DIFDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(*, '("           IYEAR1 ", I4 )') IYEAR1 
        WRITE(*, '("           IMON1 ", I2 )') IMON1 
        WRITE(*, '("           IDAY1 ", I2 )') IDAY1 
        WRITE(*, '("           IHOUR1 ", I2 )') IHOUR1 
        WRITE(*, '("           IMIN1 ", I2 )') IMIN1 
        WRITE(*, '("           ISEC1 ", I2 )') ISEC1 
        WRITE(*, '("           IYEAR2 ", I4 )') IYEAR2 
        WRITE(*, '("           IMON2 ", I2 )') IMON2 
        WRITE(*, '("           IDAY2 ", I2 )') IDAY2 
        WRITE(*, '("           IHOUR2 ", I2 )') IHOUR2 
        WRITE(*, '("           IMIN2 ", I2 )') IMIN2 
        WRITE(*, '("           ISEC2 ", I2 )') ISEC2 
        WRITE(*, '("           KSHIFT ", I12 )') KSHIFT 
        CALL ABORT
      ENDIF

      IF (LHOOK) CALL DR_HOOK('DIFDATE',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE DIFDATE


