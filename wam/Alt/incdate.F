      SUBROUTINE INCDATE (CDATE, ISHIFT)

C ----------------------------------------------------------------------
C
C**** *INCDATE* - TO UPDATE DATE TIME GROUP
C
C     L. BERTOTTI, P.JANSSEN.
C
C     H. GUNTHER   ECMWF  NOVEMBER 1989    NEGATIVE INCREMENTS.
C
C*    PURPOSE.
C     --------
C
C       UPDATING DATE TIME GROUP.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *INCDATE (CDATE,ISHIFT)*
C         *CDATE*  CHAR*10 - DATE TIME GROUP (YYMMDDHHMM)
C         *ISHIFT* INTEGER - TIME INCREMENT IN SECONDS, WHERE
C                            ABS (ISHIFT) HAS TO BE LESS THEN 1 YEAR
C
C     METHOD.
C     -------
C
C       NONE.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCES.
C     -----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
      DIMENSION MON(12)
      CHARACTER*10 CDATE
C
      DATA MON /31,28,31,30,31,30,31,31,30,31,30,31/
C
C ----------------------------------------------------------------------
C
C*    1.0 SPLITE DATE TIME GROUP INTO MINUTE, HOUR, DAY, MONTH, YEAR.
C         -----------------------------------------------------------
C
 1000 CONTINUE
C 
      READ(CDATE,5) IYEAR,MONTH,IDAY,IHOUR,MINUT
      IF(MOD(IYEAR,4).EQ.0) MON(2)=29
C
C ----------------------------------------------------------------------
C
C*    2.0 ADD SECONDS AND UPDATE DATE AND TIME.
C         -------------------------------------
C
 2000 CONTINUE
      MINUT=MINUT+ISHIFT/60
C
C     2.1 POSITIVE SHIFT GREATER THAN 1 MINUTE.
C
      IF (MINUT.GE.60) THEN
         IHOUR = IHOUR + MINUT/60
         MINUT = MINUT - (MINUT/60)*60
         IF (IHOUR.GE.24) THEN
            IDAY = IDAY + IHOUR/24
            IHOUR = IHOUR - (IHOUR/24)*24
            IF (IDAY.GT.MON(MONTH)) THEN
 1300          CONTINUE
               IDAY=IDAY-MON(MONTH)
               MONTH=MONTH+1
               IF(MONTH.EQ.13) THEN
                  MONTH = 1
                  IYEAR=MOD(IYEAR+1,100)
                  IF (MOD(IYEAR,4).EQ.0) THEN
                     MON(2)=29
                  ELSE
                     MON(2)=28
                  ENDIF
               END IF
               IF(IDAY.GT.MON(MONTH)) GO TO 1300
            END IF
         END IF
      ELSE IF (MINUT.LT.0) THEN
C
C     2.2 NEGATIVE SHIFT.
C
         IHOUR = IHOUR + (MINUT-59)/60
         MINUT = MINUT - ((MINUT-59)/60)*60
         IF (IHOUR.LT.0) THEN
            IDAY = IDAY + (IHOUR-23)/24
            IHOUR = IHOUR - ((IHOUR-23)/24)*24
            IF (IDAY.LT.1) THEN
 1400          CONTINUE
               MONTH=MONTH-1
               IF(MONTH.EQ.0) THEN
                  MONTH = 12
                  IYEAR=MOD(IYEAR+99,100)
                  IF (MOD(IYEAR,4).EQ.0) THEN
                     MON(2)=29
                  ELSE
                     MON(2)=28
                  ENDIF
               END IF
               IDAY=IDAY+MON(MONTH)
               IF(IDAY.LT.1) GO TO 1400
            END IF
         END IF
      END IF
C
C ----------------------------------------------------------------------
C
C*    3.0 COMPOSE NEW DATE TIME GROUP.
C         ----------------------------
C
 3000 CONTINUE
      WRITE(CDATE,5) IYEAR,MONTH,IDAY,IHOUR,MINUT
C
    5 FORMAT(5I2.2)
C
      RETURN
      END
