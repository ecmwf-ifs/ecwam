C======================================================================
C
      SUBROUTINE incdates (cd_date, kshift)
C
C----------------------------------------------------------------------
C
C****     *INCDATES* - TO UPDATE DATE TIME GROUP
C
C          PROGRAMMER:
C          -----------
C
C          L. BERTOTTI,P.JANSSEN.
C
C          MODIFIED BY:
C          ------------
C
C          H. GUNTHER   ECMWF  NOVEMBER 1989
C          b. hansen    ECMWF  may      1991
C
C          PURPOSE.
C          ________
C
C          UPDATING DATE TIME GROUP.
C
C**   INTERFACE. (modified date time group)
C     ----------
C
C         *CALL* *INCDATES(cd_DATE,kshift)*
C
C         *cd_date*  INTEGER - DATE TIME GROUP (YYYYMMDDhhmmss)
C         *kshift*   INTEGER - TIME INCREMENT IN SECONDS, WHERE
C                              ABS (kshift) HAS TO BE LESS THEN 1 YEAR
C
C      METHOD.
C      -------
C          NONE.
C
C      EXTERNALS.
C      ----------
C          NONE.
C
C      REFERENCES.
C      -----------
C          NONE.
C
C--------------------------------------------------------------------
C
      CHARACTER*14 cd_date
      DIMENSION i_mon(12)
C
      DATA i_mon /31,28,31,30,31,30,31,31,30,31,30,31/
C
C ----------------------------------------------------------------------
C
C*    1.0 SPLITE DATE TIME GROUP INTO MINUTE, HOUR, DAY, MONTH, YEAR.
C         -----------------------------------------------------------
C
 1000 CONTINUE
      READ(cd_date,1010) i_year,i_month,i_day,i_hour,i_minut,i_second
 1010 FORMAT(i4.4,5I2.2)
C
      IF(MOD(i_year,4).EQ.0) THEN
         i_mon(2)=29
      ELSE
         i_mon(2)=28
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    2.0 ADD SECONDS AND UPDATE DATE AND TIME.
C         -------------------------------------
C
 2000 CONTINUE
      i_second = i_second + kshift
C
C     2.1 POSITIVE SHIFT GREATER THAN 1 MINUTE.
C
      IF (i_second .GE. 60) THEN
         i_minut = i_minut + i_second/60
         i_second = i_second - (i_second/60)*60
         IF (i_minut.GE.60) THEN
            i_hour = i_hour + i_minut/60
            i_minut = i_minut - (i_minut/60)*60
            IF (i_hour.GE.24) THEN
               i_day = i_day + i_hour/24
               i_hour = i_hour - (i_hour/24)*24
               IF (i_day.GT.i_mon(i_month)) THEN
 2100             CONTINUE
                  i_day=i_day-i_mon(i_month)
                  i_month=i_month+1
                  IF(i_month.EQ.13) THEN
                     i_month = 1
                     i_year=MOD(i_year+1,10000)
                     IF (MOD(i_year,4).EQ.0) THEN
                        i_mon(2)=29
                     ELSE
                        i_mon(2)=28
                     ENDIF
                  ENDIF
                  IF(i_day.GT.i_mon(i_month)) GO TO 2100
               ENDIF
            ENDIF
         ENDIF
      ELSE IF (i_second.LT.0) THEN
C
C     2.2 NEGATIVE SHIFT.
C
         i_minut = i_minut + (i_second-59)/60
         i_second = i_second - ((i_second-59)/60)*60
         IF (i_minut .LT. 0) THEN
            i_hour = i_hour + (i_minut-59)/60
            i_minut = i_minut - ((i_minut-59)/60)*60
            IF (i_hour.LT.0) THEN
               i_day = i_day + (i_hour-23)/24
               i_hour = i_hour - ((i_hour-23)/24)*24
               IF (i_day.LT.1) THEN
 2200             CONTINUE
                  i_month=i_month-1
                  IF(i_month.EQ.0) THEN
                     i_month = 12
                     i_year=MOD(i_year+9999,10000)
                     IF (MOD(i_year,4).EQ.0) THEN
                        i_mon(2)=29
                     ELSE
                        i_mon(2)=28
                     ENDIF
                  END IF
                  i_day=i_day+i_mon(i_month)
                  IF(i_day.LT.1) GO TO 2200
               END IF
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
      WRITE(cd_date,1010) i_year,i_month,i_day,i_hour,i_minut,i_second
C
      RETURN
      END
