C ======================================================================
C
      SUBROUTINE DIFDATES (cd_date1, cd_date2, kshift)
C
C ======================================================================
C
C****     *DIFDATE* - TO COMPUTE TIME DIFFERENCE.
C
C          PROGRAMMER:
C          -----------
C
C          H. GUNTHER   ECMWF  NOVEMBER 1989
C
C          PURPOSE.
C          ________
C
C          COMPUTE THE SECONDS BETWEEN THE INPUT DATES.
C          DATES HAVE TO BE IN CONSECUTIVE YEARS.
C
C**   INTERFACE.
C     ----------
C
C         *CALL* *DIFDATE (cd_date1, cd_date2, kshift)*
C
C     I/   *cd_date1* INTEGER - DATE TIME GROUP (YYMMDDHHMM).
C     I/   *cd_date2* INTEGER - DATE TIME GROUP (YYMMDDHHMM).
C      /O  *kshift* INTEGER - DIFFERENCE IN SECONDS (cd_date2-cd_date1).
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
C      MODIFIKATIONS
C      -----------
C          NEW FORMAT FOR THE INPUT VARIABLES cd_date1 AND cd_date2
C         *cd_date1* INTEGER - DATE TIME GROUP (YYYYMMDDhhhhss).
C         *cd_date2* INTEGER - DATE TIME GROUP (YYYYMMDDhhmmss).
C
C
C
C
C
C ---------------------------------------------------------------------
C
      CHARACTER*14 cd_date1, cd_date2
      CHARACTER*14 cldt1, cldt2, cldat
      INTEGER i_sec_1, i_sec_2, isi
C
      DIMENSION i_mon(12)
C
      DATA i_mon /31,28,31,30,31,30,31,31,30,31,30,31/
C
C ---------------------------------------------------------------------
C
C*    1.0 CHANGE DATE TIME GROUPS TO ENSURE THAT THE SECOND IS LARGER.
C         ------------------------------------------------------------
C
 1000 CONTINUE
C
      cldt1 = cd_date1
      cldt2 = cd_date2
      IF (cldt1 .GT. cldt2) THEN
         cldt1 = cd_date2
         cldt2 = cd_date1
         isi  = -60
      ELSE
         isi  =  60
      ENDIF
C
C ---------------------------------------------------------------------
C
C*    2.0 SPLITE DATE TIME GROUPS INTO MINUTE, HOUR, DAY, MONTH, YEAR.
C         ------------------------------------------------------------
C
 2000 CONTINUE
      READ(cldt1,2100) i_yy_1,i_mm_1,i_dd_1,i_hh_1,i_min_1,i_sec_1
      READ(cldt2,2100) i_yy_2,i_mm_2,i_dd_2,i_hh_2,i_min_2,i_sec_2
 2100 FORMAT(i4.4,5I2.2)
C
C ---------------------------------------------------------------------
C
C*    3.0 CORRECT DAYS IN FEBRUARY OF FIRST YEAR FOR LEAP-YEAR.
C         -----------------------------------------------------
C
 3000 CONTINUE
      IF (MOD(i_yy_1,4).EQ.0) THEN
         i_mon(2) = 29
      ELSE
         i_mon(2) = 28
      ENDIF
C
C ---------------------------------------------------------------------
C
C*    4.0 COMPUTE TIME DIFFERENCE IN MINUTES.
C         -----------------------------------
C
 4000 CONTINUE
C
C     4.1 DIFFERENCE BETWEEN DAY, HOUR ,MINITE.
C
      kshift = ((i_dd_2-i_dd_1)*24+i_hh_2-i_hh_1)*60+i_min_2-i_min_1
C
C      4.2 ADD DIFFERENCE FROM MONTH.
C
      IF (i_yy_2 .GT. i_yy_1) THEN
C
C      4.2.1 START AND END MONTH ARE IN DIFFERENT YEARS.
C
         DO 4211 jm=i_mm_1,12
            kshift = kshift + i_mon(jm)*24*60
 4211    CONTINUE
         IF (i_mm_2.GT.1) THEN
            IF (MOD(i_yy_2,4).EQ.0) THEN
               i_mon(2) = 29
            ELSE
               i_mon(2) = 28
            ENDIF
            DO 4212 jm=1,i_mm_2-1
               kshift = kshift + i_mon(jm)*24*60
 4212       CONTINUE
         ENDIF
      ELSE
C
C      4.2.2 START AND END MONTH ARE IN THE SAME YEAR.
C
         IF (i_mm_2.GT.i_mm_1) THEN
            DO 4221 jm=i_mm_1,i_mm_2-1
               kshift = kshift + i_mon(jm)*24*60
 4221       CONTINUE
         ENDIF
      ENDIF
C
C ---------------------------------------------------------------------
C
C*    5.0 CHANGE SIGN OF DIFFERENCE IF DATES HAVE BEEN EXCHANGED
C*        AND CONVERT TO SECONDS.
C         ------------------------------------------------------
C
 5000 CONTINUE
      kshift = kshift * isi - i_sec_1 + i_sec_2
C
      RETURN
      END
