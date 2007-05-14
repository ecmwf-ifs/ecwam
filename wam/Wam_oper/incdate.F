! ======================================================================

      SUBROUTINE INCDATE (CDATE, ISHIFT)

! ======================================================================

!**** *INCDATE* - TO UPDATE DATE TIME GROUP

!     B. HANSEN           ECMWF       JULY    1998

!**   PURPOSE.
!     --------
!       UPDATING DATE TIME GROUP.

!**   INTERFACE.
!     ----------

!       *CALL* *INCDATE (CDATE,ISHIFT)*
!         *CDATE*  CHAR*12 - DATE TIME GROUP (YYYYMMDDHHMM) OR
!                  CHAR*14 - DATE TIME GROUP (YYYYMMDDHHMMSS)
!         *ISHIFT* INTEGER - TIME INCREMENT IN SECONDS, RESTRICTIONS:
!                            -2147483648 <= ISHIFT <=  +2147483647

!     METHOD.
!     -------

!       THE DATE AND TIME CAN BE SUPPLIED IN THE Y2K COMPLIANT
!       14 CHARACTER OR 12 CHARACTER FORMAT. IF THE 10 CHARACTER
!       FORMAT IS USED, AN ERROR MESSAGE WILL BE ISSUED. 
!       THE CHARACTER STRING IS BROKEN INTO INTEGERS TO OBTAIN 
!       THE INPUT ARGUMENT FOR SECINCR. THE UPDATED DATE AND TIME 
!       IS PUT BACK INTO AN OUTPUT STRING AS USED BY THE CALLING 
!       ROUTINE.

!     EXTERNALS.
!     ----------

!       SECINCR (ECLIB).

!     REFERENCES.
!     -----------

!       OD MEMO (6 July 1998) File: O/PRO/90 FROM GROUP 2000.

! ----------------------------------------------------------------------
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

      CHARACTER*(*) CDATE
      CHARACTER* 80 CLFMT
      INTEGER ISHIFT
      LOGICAL LLND
      REAL ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('INCDATE',0,ZHOOK_HANDLE)

! ----------------------------------------------------------------------

!*    1.0 SPLITE DATE TIME GROUP INTO SECOND, MINUTE, HOUR, DAY, MONTH,
!*        YEAR DEPENDENT ON THE FORMAT OF CDATE.
!         -----------------------------------------------------------

      IL = LEN_TRIM(CDATE)
      IF (IL==10) THEN
        WRITE(CLFMT,'(A49, I2.2, A2)')
     &     '(" INCDATE : NON-Y2K COMPLIANT  DATE IN USE ", A', IL, ' )'
        WRITE(*, CLFMT ) CDATE
        LLND = .FALSE.
        CALL ABORT
      ELSEIF (IL==12) THEN
        LLND = .TRUE.
        READ(CDATE,'(i4, 4i2)')IYEAR, IMON, IDAY, IHOUR, IMIN
        ISEC = 0
      ELSEIF (IL==14) THEN
        LLND = .TRUE.
        READ(CDATE,'(i4, 5i2)')IYEAR, IMON, IDAY, IHOUR, IMIN, ISEC
      ELSEIF (IL==0) THEN
        WRITE(*, *) ' INCDATE:  FATAL@! DATE IS 0 CHARACTER LONG'
        CALL ABORT
      ELSE
        WRITE(CLFMT,'(A59, I2.2, A2)')
     &  '(" INCDATE:  FATAL@! DATE IS ",I2," CHARACTERS LONG!!= ",A',
     &    IL, ' )'
        WRITE(*, CLFMT) IL, CDATE
        STOP
      ENDIF

! ----------------------------------------------------------------------

!*    2.0 INCREMENT OR DECREMENT THE DATE AND TIME BY ISHIFT SECONDS.
!         -----------------------------------------------------------

      CALL SECINCR  (IYEAR, IMON, IDAY, IHOUR, IMIN, ISEC, ISHIFT,
     &               INYEAR, INMON, INDAY, INHOUR, INMIN, INSEC, IRET)

      IF (IRET/=0) THEN
        WRITE(*, '(" INCDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(*, '("           IYEAR ", I4 )') IYEAR
        WRITE(*, '("           IMON ", I2 )') IMON
        WRITE(*, '("           IDAY ", I2 )') IDAY
        WRITE(*, '("           IHOUR ", I2 )') IHOUR
        WRITE(*, '("           IMIN ", I2 )') IMIN
        WRITE(*, '("           ISEC ", I2 )') ISEC
        WRITE(*, '("           ISHIFT ", I12 )') ISHIFT
        CALL ABORT
      ENDIF

! ----------------------------------------------------------------------

!*    3.0 CREATE OUTPUT STRING OF NEW DATE AND TIME USING THE INPUT
!*        FORMAT.
!         ---------------------------------------------------------

      IF (LLND) THEN
        IF (IL==12) THEN
          WRITE(CDATE,'(i4.4, 4i2.2)')INYEAR, INMON, INDAY, INHOUR,
     &                                INMIN
        ELSEIF (IL==14) THEN
          WRITE(CDATE,'(i4.4, 5i2.2)')INYEAR, INMON, INDAY, INHOUR,
     &                                INMIN, INSEC
        ELSE
          WRITE(*,'(" THE LENGTH OF THE INPUT CHARACTER CHANGED @!")')
          CALL ABORT
        ENDIF
      ELSE
        INYEAR = INYEAR - 1900
        IF (INYEAR >= 100) THEN
          WRITE(CLFMT,'(A39, I2.2, A2)')
     &     '(" NON-Y2K COMPLIANT  DATE IN USE ", A', IL, ' )'
          WRITE(*, CLFMT ) CDATE
          WRITE(*,'(" THIS CANNOT BE DONE ANY LONGER")')
          CALL ABORT
        ENDIF
        WRITE(CDATE,'(5i2.2)')INYEAR, INMON, INDAY, INHOUR, INMIN
      ENDIF
      IF (LHOOK) CALL DR_HOOK('INCDATE',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE INCDATE


