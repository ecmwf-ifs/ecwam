      SUBROUTINE READSTRESS(RU10OLD,RTHWOLD,RUSOLD,RTAUW,RZ0OLD,
     &                      FILENAME)

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      MARCH 1997

!*    PURPOSE.
!     --------
!     READS THE RESTART WIND AND STRESS FIELDS

!**   INTERFACE.
!     ----------
!     CALL *READSTRESS*(RU10OLD,RTHWOLD,RUSOLD,RTAUW,RZ0OLD,FILENAME)
!     *RU10OLD*   WIND SPEED.
!     *RTHWOLD*   WIND DIRECTION (RADIANS).
!     *RUSOLD*    FRICTION VELOCITY.
!     *RTAUW*     WAVE STRESS.
!     *RZ0OLD*    ROUGHNESS LENGTH IN.
!     *FILENAME* FILENAME (INCLUDING PATH) OF INPUT FILE (ONLY IF 
!                MESSAGE PASSING IS USED)

!     METHOD.
!     -------
!     READS ARRAYS FROM FILE. FORTRAN READS OR PBIO READ ARE
!     USED DEPENDING ON THE CHOICE OF OUTPUT MODE AS SELECTED
!     BY LPBIOOUT

!     IF PBIO IS USED, THEN PBREAD IS CALLED TO READ EACH ARRAY.
!     IF FORTRAN READ IS USED THEN FILNAME IS READ AS UNFORMATTED BINARY
!     ONLY THE NECESSARY CONTRIBUTION FROM FL IS KEPT ON PE IN MESSAGE.

!     EXTERNALS.
!     ----------
!     *PBOPEN*
!     *PBREAD*
!     *PBCLOSE*
!     *ABORT1*

!     REFERENCE.
!     ----------
!     NONE
! ----------------------------------------------------------------------

      USE YOWMESPAS, ONLY : LPBIOIN
      USE YOWMPP   , ONLY : NPRECR
      USE YOWPARAM , ONLY : NBLO     ,NIBLO
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWUNIT  , ONLY : IU15
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL

! ----------------------------------------------------------------------
!     MODIFICATIONS.
!     --------------
!     B. HANSEN    ECMWF      NOVEMBER 1998  KOUNT IS 12 FOR Y2K.

      CHARACTER*1 MODE
      CHARACTER*120 FILENAME

      REAL,DIMENSION(NIBLO,NBLO) :: RU10OLD,RTHWOLD,RUSOLD,RZ0OLD,RTAUW
! ----------------------------------------------------------------------

      IF (ITEST.GE.2) THEN
          WRITE(IU06,*) ' SUB. READSTRESS: ABOUT TO INPUT RESTART FILES'
          WRITE(IU06,*) ' LPBIOIN=', LPBIOIN
          CALL FLUSH (IU06)
      ENDIF

      IF(LPBIOIN) THEN

        MODE = 'r'
        CALL PBOPEN(IUNIT,FILENAME,MODE,KRET)
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN       *'
          WRITE (IU06,*) '*   IN READSTRESS                        *' 
          IF(KRET.EQ.-1)
     &      WRITE (IU06,*) 'COULD NOT OPEN FILE ',FILENAME
          IF(KRET.EQ.-2) WRITE (IU06,*) 'INVALID FILENAME ',FILENAME
          IF(KRET.EQ.-3) WRITE (IU06,*) 'INVALID OPEN MODE SPECIFIED'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          CALL ABORT1
        ENDIF


        KOUNT = 12 
        CALL PBREAD(IUNIT,CDTPRO,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        CALL PBREAD(IUNIT,CDATEWO,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        CALL PBREAD(IUNIT,CDAWIFL,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        CALL PBREAD(IUNIT,CDATEFL,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        KOUNT = NIBLO*NBLO*NPRECR
        CALL PBREAD(IUNIT,RU10OLD,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        CALL PBREAD(IUNIT,RTHWOLD,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        CALL PBREAD(IUNIT,RUSOLD,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        CALL PBREAD(IUNIT,RTAUW,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        CALL PBREAD(IUNIT,RZ0OLD,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

1099    CONTINUE
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '*****************************************'
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*    ERROR FOLLOWING CALL TO PBREAD     *'
          WRITE (IU06,*) '*    IN READSTRESS                      *'
          WRITE (IU06,*) '*    FILE ',FILENAME
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*    KRET = ', KRET
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*****************************************'
          CALL ABORT1
        ELSEIF(KRET.NE.KOUNT) THEN
          WRITE (IU06,*) '*****************************************'
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*    ERROR FOLLOWING CALL TO PBREAD     *'
          WRITE (IU06,*) '*    IN READSTRESS                      *'
          WRITE (IU06,*) '*    FILE ',FILENAME
          WRITE (IU06,*) '*    NUMBER OF BYTES READ : ',KRET
          WRITE (IU06,*) '*    IS NOT EQUAL TO THE NUMBER OF BYTES*'
          WRITE (IU06,*) '*    REQUESTED : ',KOUNT
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*****************************************'
          CALL ABORT1
        ENDIF

        CALL PBCLOSE(IUNIT,KRET)
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '************************************'
          WRITE (IU06,*) '*                                  *'
          WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
          WRITE (IU06,*) '* IN READSTRESS                    *'
          WRITE (IU06,*) '* FILE ',FILENAME
          WRITE (IU06,*) '*                                  *'
          WRITE (IU06,*) '************************************'
          CALL ABORT1
        ENDIF

      ELSE

        READ(IU15) CDTPRO,CDATEWO,CDAWIFL,CDATEFL
        READ(IU15) RU10OLD
        READ(IU15) RTHWOLD
        READ(IU15) RUSOLD
        READ(IU15) RTAUW
        READ(IU15) RZ0OLD

        CLOSE(IU15)

      ENDIF

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) ' SUB. READSTRESS: INPUT OF RESTART FILES DONE'
        CALL FLUSH(IU06)
      ENDIF

      RETURN
      END SUBROUTINE READSTRESS
