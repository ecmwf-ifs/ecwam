C ----------------------------------------------------------------------
      SUBROUTINE DIFDATE (CDATE1, CDATE2, ISHIFT)
C ----------------------------------------------------------------------
C
C**** *DIFDATE* - TO COMPUTE TIME DIFFERENCE.
C
C     H. GUNTHER   GKSS/ECMWF  NOVEMBER 1989
C
C*    PURPOSE.
C     --------
C
C       COMPUTE THE SECONDS BETWEEN THE INPUT DATES.
C       DATES HAVE TO BE IN CONSECUTIVE YEARS.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *DIFDATE (CDATE1, CDATE2, ISHIFT)*
C
C         *CDATE1* - DATE TIME GROUP (YYMMDDHHMM).
C         *CDATE2* - DATE TIME GROUP (YYMMDDHHMM).
C         *ISHIFT* - DIFFERENCE IN SECONDS (CDATE2-CDATE1).
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
C      REFERENCES.
C      -----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
      DIMENSION MON(12)
C
      CHARACTER*10 CDATE1, CDATE2, CDT1, CDT2, CDAT
C
      DATA MON /31,28,31,30,31,30,31,31,30,31,30,31/
C
C ----------------------------------------------------------------------
C
C*    1.0 CHANGE DATE TIME GROUPS TO ENSURE THAT THE SECOND IS LARGER.
C         ------------------------------------------------------------
C
C
C
      CDT1=CDATE1
      CDT2=CDATE2
      ISI = 60
C
C         END OF CENTURY CORRECTION.
C
      READ(CDT1(1:4),'(I4)') IYM1
      READ(CDT2(1:4),'(I4)') IYM2
      IF (IYM1.LT. 100 .AND. IYM2.GT.9900)
     1    IYM1 = 10000 + IYM1
      IF (IYM2.LT. 100 .AND. IYM1.GT.9900)
     1    IYM2 = 10000 + IYM2

      IF (CDT1.GT.CDT2.OR.IYM1.GT.IYM2) THEN
         CDT1 = CDATE2
         CDT2 = CDATE1
         ISI  = -60
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    2.0 SPLIT DATE TIME GROUPS INTO MINUTE, HOUR, DAY, MONTH, YEAR.
C         -----------------------------------------------------------
C
      READ(CDT1,5) IYEAR1,MONTH1,IDAY1,IHOUR1,MINUT1
      IF (IYM1.GT.10000) IYEAR1=IYEAR1+100
C
      READ(CDT2,5) IYEAR2,MONTH2,IDAY2,IHOUR2,MINUT2
      IF (IYM2.GT.10000) IYEAR2=IYEAR2+100
C
C ----------------------------------------------------------------------
C
C*    3.0 CORRECT DAYS IN FEBRUARY OF FIRST YEAR FOR LEAP-YEAR.
C         -----------------------------------------------------
C
      IF (MOD(IYEAR1,4).EQ.0) MON(2) = 29
C
C ----------------------------------------------------------------------
C
C*    4.0 COMPUTE TIME DIFFERENCE IN MINUTES.
C         -----------------------------------
C
C     4.1 DIFFERENCE BETWEEN DAY, HOUR ,MINITE.
C
      ISHIFT = ((IDAY2-IDAY1)*24+IHOUR2-IHOUR1)*60+MINUT2-MINUT1
C
C      4.2 ADD DIFFERENCE FROM MONTH.
C
      IF (IYEAR2.GT.IYEAR1) THEN
C
C      4.2.1 START AND END MONTH ARE IN DIFFERENT YEARS.
C
         DO WHILE (IYEAR2 .GT. IYEAR1)
            DO 4211 M=MONTH1,12
               ISHIFT = ISHIFT + INT(MON(M)*24*60)
 4211       CONTINUE
            MONTH1 = 1
            IYEAR1 = IYEAR1 + 1
            IF (MOD(IYEAR1,4).EQ.0) THEN
               MON(2) = 29
            ELSE
               MON(2) = 28
            ENDIF

         ENDDO
         IF (MONTH2.GT.1) THEN
            IF (MOD(IYEAR2,4).EQ.0) THEN
               MON(2) = 29
            ELSE
               MON(2) = 28
            ENDIF
            DO 4212 M=1,MONTH2-1
               ISHIFT = ISHIFT + INT(MON(M)*24*60)
 4212       CONTINUE
         ENDIF
      ELSE
C
C      4.2.2 START AND END MONTH ARE IN THE SAME YEAR.
C
         IF (MONTH2.GT.MONTH1) THEN
            DO 4221 M=MONTH1,MONTH2-1
               ISHIFT = ISHIFT + INT(MON(M)*24*60)
 4221       CONTINUE
         ENDIF
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    5.0 CHANGE SIGN OF DIFFERENCE IF DATES HAVE BEEN EXCHANGED
C*        AND CONVERT TO SECONDS.
C         ------------------------------------------------------
C
      ISHIFT = ISHIFT*ISI
C
    5 FORMAT(5I2.2)
C
      RETURN
      END
