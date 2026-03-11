! ======================================================================
! TEMPLATE FOR INCDATE IMPLEMENTATION
! This file is included by specific typed versions of INCDATE
! The following must be defined before including:
!   - INCDATE_ROUTINE: name of the specific routine
!   - ISHIFT_KIND: the integer kind for ISHIFT argument
! ======================================================================

!*    1.0 SPLIT DATE TIME GROUP INTO SECOND, MINUTE, HOUR, DAY, MONTH,
!*        YEAR DEPENDENT ON THE FORMAT OF CDATE.
!         -----------------------------------------------------------

      IL = LEN_TRIM(CDATE)
      IF (IL==10) THEN
        WRITE(CLFMT,'(A49, I2.2, A2)')                                  &
     &     '(" INCDATE : NON-Y2K COMPLIANT  DATE IN USE ", A', IL, ' )'
        WRITE(6, CLFMT ) CDATE
        WRITE(*, CLFMT ) CDATE
        LLND = .FALSE.
        CALL ABORT1
      ELSEIF (IL==12) THEN
        LLND = .TRUE.
        READ(CDATE,'(i4, 4i2)')IYEAR, IMON, IDAY, IHOUR, IMIN
        ISEC = 0
      ELSEIF (IL==14) THEN
        LLND = .TRUE.
        READ(CDATE,'(i4, 5i2)')IYEAR, IMON, IDAY, IHOUR, IMIN, ISEC
      ELSEIF (IL==0) THEN
        WRITE(6, *) ' INCDATE:  FATAL@! DATE IS 0 CHARACTER LONG'
        WRITE(*, *) ' INCDATE:  FATAL@! DATE IS 0 CHARACTER LONG'
        CALL ABORT1
      ELSE
        WRITE(CLFMT,'(A59, I2.2, A2)')                                  &
     &  '(" INCDATE:  FATAL@! DATE IS ",I2," CHARACTERS LONG!!= ",A',   &
     &    IL, ' )'
        WRITE(6, CLFMT) IL, CDATE
        WRITE(*, CLFMT) IL, CDATE
        STOP
      ENDIF

      IRET=0
      MON(2)=MFEB_LENGTH(IYEAR)
      IF (IMON<1 .OR. IMON>12) IRET=-1
      IF (IDAY<1 .OR. IDAY>MON(MIN(MAX(IMON,1),12))) IRET=-2
      IF (IHOUR<0 .OR. IHOUR>23) IRET=-3
      IF (IMIN<0 .OR. IMIN>59) IRET=-4
      IF (ISEC<0 .OR. ISEC>59) IRET=-5
      IF (IRET/=0) THEN
        WRITE(6, '(" INCDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(6, '("           INPUT DATE INCORRECT,")')
        WRITE(6, '("           IYEAR ", I4 )') IYEAR
        WRITE(6, '("           IMON ", I2 )') IMON
        WRITE(6, '("           IDAY ", I2 )') IDAY
        WRITE(6, '("           IHOUR ", I2 )') IHOUR
        WRITE(6, '("           IMIN ", I2 )') IMIN
        WRITE(6, '("           ISEC ", I2 )') ISEC
        WRITE(*, '(" INCDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(*, '("           INPUT DATE INCORRECT,")')
        WRITE(*, '("           IYEAR ", I4 )') IYEAR
        WRITE(*, '("           IMON ", I2 )') IMON
        WRITE(*, '("           IDAY ", I2 )') IDAY
        WRITE(*, '("           IHOUR ", I2 )') IHOUR
        WRITE(*, '("           IMIN ", I2 )') IMIN
        WRITE(*, '("           ISEC ", I2 )') ISEC
        CALL ABORT1
      ENDIF
      IF (IL==12 .AND. MOD(ISHIFT,60) /= 0) THEN
        IRET=-6
        WRITE(6, '(" INCDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(6, '(" WHEN HANDELING CHARACTER*12 DATES,")')
        WRITE(6, '(" ISHIFT MUST BE A MULTIPLE OF 60 SECONDS")')
        WRITE(6, '(" ISHIFT = ", I12 )') ISHIFT
        WRITE(*, '(" INCDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(*, '(" WHEN HANDELING CHARACTER*12 DATES,")')
        WRITE(*, '(" ISHIFT MUST BE A MULTIPLE OF 60 SECONDS")')
        WRITE(*, '(" ISHIFT = ", I12 )') ISHIFT
        CALL ABORT1
      ENDIF

! ----------------------------------------------------------------------

!*    2.0 INCREMENT OR DECREMENT THE DATE AND TIME BY ISHIFT SECONDS.
!         -----------------------------------------------------------


      ISEC=ISEC+MOD(ISHIFT,INT(60,KIND=ISHIFT_KIND))
      IF (ISEC >= 60) THEN
        IMIN = IMIN+ISEC/60
        ISEC = ISEC-(ISEC/60)*60
      ELSE IF (ISEC < 0) THEN
        IMIN = IMIN +(ISEC-59)/60
        ISEC = ISEC-((ISEC-59)/60)*60
      END IF

      IMIN=IMIN+ISHIFT/INT(60,KIND=ISHIFT_KIND)

!     2.1 POSITIVE SHIFT GREATER THAN 1 MINUTE.

      IF (IMIN >= 60) THEN
         IHOUR = IHOUR + IMIN/60
         IMIN = IMIN - (IMIN/60)*60
         IF (IHOUR >= 24) THEN
            IDAY = IDAY + IHOUR/24
            IHOUR = IHOUR - (IHOUR/24)*24
            DO WHILE (IDAY > MON(IMON))
               IDAY=IDAY-MON(IMON)
               IMON=IMON+1
               IF (IMON == 13) THEN
                  IMON = 1
                  IYEAR=IYEAR+1
                  MON(2)=MFEB_LENGTH(IYEAR)
               END IF
            ENDDO
         END IF
      ELSE IF (IMIN < 0) THEN
 
!     2.2 NEGATIVE SHIFT.

         IHOUR = IHOUR + (IMIN-59)/60
         IMIN = IMIN - ((IMIN-59)/60)*60
         IF (IHOUR < 0) THEN
            IDAY = IDAY + (IHOUR-23)/24
            IHOUR = IHOUR - ((IHOUR-23)/24)*24
            DO WHILE (IDAY < 1)
               IMON=IMON-1
               IF (IMON == 0) THEN
                  IMON = 12
                  IYEAR=IYEAR-1
                  MON(2)=MFEB_LENGTH(IYEAR)
               END IF
               IDAY=IDAY+MON(IMON)
            ENDDO
         END IF
      END IF


! ----------------------------------------------------------------------

!*    3.0 CREATE OUTPUT STRING OF NEW DATE AND TIME USING THE INPUT
!*        FORMAT.
!         ---------------------------------------------------------

      IRET=0
      MON(2)=MFEB_LENGTH(IYEAR)
      IF (IMON<1 .OR. IMON>12) IRET=-1
      IF (IDAY<1 .OR. IDAY>MON(MIN(MAX(IMON,1),12))) IRET=-2
      IF (IHOUR<0 .OR. IHOUR>23) IRET=-3
      IF (IMIN<0 .OR. IMIN>59) IRET=-4
      IF (ISEC<0 .OR. ISEC>59) IRET=-5
      IF (IRET/=0) THEN
        WRITE(6, '(" INCDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(6, '("           OUTPUT DATE INCORRECT,")')
        WRITE(6, '("           IYEAR ", I4 )') IYEAR
        WRITE(6, '("           IMON ", I2 )') IMON
        WRITE(6, '("           IDAY ", I2 )') IDAY
        WRITE(6, '("           IHOUR ", I2 )') IHOUR
        WRITE(6, '("           IMIN ", I2 )') IMIN
        WRITE(6, '("           ISEC ", I2 )') ISEC
        WRITE(6, '("           ISHIFT ", I12 )') ISHIFT
        CALL ABORT1
      ENDIF

      IF (LLND) THEN
        IF (IL==12) THEN
          WRITE(CDATE,'(I4.4, 4I2.2)')IYEAR, IMON, IDAY, IHOUR,         &
     &                                IMIN
        ELSEIF (IL==14) THEN
          WRITE(CDATE,'(I4.4, 5I2.2)')IYEAR, IMON, IDAY, IHOUR,         &
     &                                IMIN, ISEC
        ELSE
          WRITE(6,'(" THE LENGTH OF THE INPUT CHARACTER CHANGED @!")')
          CALL ABORT1
        ENDIF
      ELSE
        IYEAR = IYEAR - 1900
        IF (IYEAR >= 100) THEN
          WRITE(CLFMT,'(A39, I2.2, A2)')                                &
     &     '(" NON-Y2K COMPLIANT  DATE IN USE ", A', IL, ' )'
          WRITE(6, CLFMT ) CDATE
          WRITE(6,'(" THIS CANNOT BE DONE ANY LONGER")')
          CALL ABORT1
        ENDIF
        WRITE(CDATE,'(5I2.2)')IYEAR, IMON, IDAY, IHOUR, IMIN
      ENDIF
