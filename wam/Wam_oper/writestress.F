      SUBROUTINE WRITESTRESS(IJINF,IJSUP,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     &                       ROAIRO,ZIDLOLD,FILENAME,LPBOPEN,LPBCLOSE) 

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      APRIL 1997

!*    PURPOSE.
!     --------
!     WRITES THE RESTART WIND AND STRESS FIELDS TO FILE.

!**   INTERFACE.
!     ----------

!     CALL *WRITESTRESS*(IJINF,IJSUP,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
!                        FILENAME,LPBOPEN,LPBCLOSE) 
!     *IJINF*    FIRST LOWER DIMENSION OF U10OLD,THWOLD,USOLD,TAUW,Z0OLD
!     *IJSUP*    FIRST UPPER DIMENSION OF U10OLD,THWOLD,USOLD,TAUW,Z0OLD
!     *U10OLD*   WIND SPEED.
!     *THWOLD*   WIND DIRECTION (RADIANS).
!     *USOLD*    FRICTION VELOCITY.
!     *TAUW*     WAVE STRESS.
!     *Z0OLD*    ROUGHNESS LENGTH IN.
!     *ROAIRO*    AIR/WATER DENSITY RATIO.
!     *ZIDLOLD*   Zi/L (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!     *FILENAME* FILENAME (INCLUDING PATH) OF TARGET FILE
!     *LPBOPEN*  LOGICAL, TRUE IF PBOPEN HAS TO BE CALLED
!     *LPBCLOSE* LOGICAL, TRUE IF PBCLOSE HAS TO BE CALLED

!     METHOD.
!     -------
!     WRITES ARRAYS TO FILE. FORTRAN WRITES OR PBIO WRITES ARE
!     USED DEPENDING ON THE CHOICE OF OUTPUT MODE AS SELECTED BY
!     LPBIOOUT

!     IF PBIO IS USED, THEN THE ARRAYS  WILL BE SAVED IN FILENAME 
!     USING PBWRITE ON CONDITION THAT FILENAME WAS NOT OPENED BEFORE
!     AND CONNECTED TO IUNIT.
!     IT WILL BE APPENDED TO THE BOTTOM OF THE CONTENT OF FILENAME IF
!     FILENAME WAS ALREADY OPENED AND CONNECTED TO UNIT IUNIT.

!     IF FORTRAN WRITES ARE USED THEN IN THE CASE OF MESSAGE PASSING
!     THE OUTPUT IS WRITTEN AS UNFORMATTED BINARY TO FILENAME CONNECTED
!     TO UNIT IU15. IN THE NON MESSAGE PASSING CASE, IT IS WRITTEN AS
!     UNFORMATTED BINARY TO FORTRAN UNIT SPECIFIED BY IU15.


!     EXTERNALS.
!     ----------
!     *PBOPEN*
!     *PBWRITE*
!     *PBCLOSE*
!     *ABORT1*

!     REFERENCE.
!     ----------
!     NONE
! ----------------------------------------------------------------------

      USE YOWMESPAS, ONLY : LMESSPASS,LPBIOOUT
      USE YOWMPP   , ONLY : NPRECR
      USE YOWPARAM , ONLY : NBLO
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWUNIT  , ONLY : IU15
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL
      USE YOMGSTATS, ONLY : LSYNCSTATS

! ----------------------------------------------------------------------

      CHARACTER*1 MODE
      CHARACTER*120 FILENAME
      LOGICAL LPBOPEN,LPBCLOSE
      REAL,DIMENSION(IJINF:IJSUP,NBLO) :: U10OLD,THWOLD,USOLD,
     &                                    Z0OLD,TAUW,ROAIRO,ZIDLOLD
! ----------------------------------------------------------------------

      IF(LPBIOOUT) THEN

        IF(LPBOPEN) THEN
          MODE = 'w'
          CALL PBOPEN(IUNIT,FILENAME,MODE,KRET)
          IF(KRET.LT.0) THEN
            WRITE (IU06,*) '******************************************'
            WRITE (IU06,*) '*                                        *'
            WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN       *'
            WRITE (IU06,*) '*   IN WRITESTRESS                       *' 
            IF(KRET.EQ.-1)
     &       WRITE (IU06,*) 'COULD NOT OPEN FILE ',FILENAME
            IF(KRET.EQ.-2) WRITE (IU06,*) 'INVALID FILENAME ',FILENAME
            IF(KRET.EQ.-3) WRITE (IU06,*) 'INVALID OPEN MODE SPECIFIED'
            WRITE (IU06,*) '*                                        *'
            WRITE (IU06,*) '******************************************'
            CALL ABORT1
          ENDIF
        ENDIF

        CALL PBWRITE(IUNIT,CDTPRO,12,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,CDATEWO,12,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,CDAWIFL,12,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,CDATEFL,12,KRET)
        IF(KRET.LT.0) GOTO 1099

        KOUNT = (IJSUP-IJINF+1)*NBLO*NPRECR
        CALL PBWRITE(IUNIT,U10OLD,KOUNT,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,THWOLD,KOUNT,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,USOLD,KOUNT,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,TAUW,KOUNT,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,Z0OLD,KOUNT,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,ROAIRO,KOUNT,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,ZIDLOLD,KOUNT,KRET)
        IF(KRET.LT.0) GOTO 1099

1099    CONTINUE 
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '*****************************************'
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*    ERROR FOLLOWING CALL TO PBWRITE    *'
          WRITE (IU06,*) '*    IN WRITESTRESS                     *' 
          WRITE (IU06,*) '*    FILE ',FILENAME 
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*****************************************'
          CALL ABORT1
        ENDIF


        IF(LPBCLOSE) THEN
          CALL PBCLOSE(IUNIT,KRET)
          IF(KRET.LT.0) THEN
            WRITE (IU06,*) '************************************'
            WRITE (IU06,*) '*                                  *'
            WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
            WRITE (IU06,*) '* IN WRITESTRESS                   *'
            WRITE (IU06,*) '* FILE ',FILENAME
            WRITE (IU06,*) '*                                  *'
            WRITE (IU06,*) '************************************'
            CALL ABORT1
          ENDIF
        ENDIF


      ELSE
        IF(LMESSPASS) THEN

          OPEN(IU15,FILE=FILENAME,STATUS='UNKNOWN',FORM='UNFORMATTED')

          WRITE(IU15) CDTPRO,CDATEWO,CDAWIFL,CDATEFL
          WRITE(IU15)((U10OLD(J1,IG),J1=IJINF,IJSUP),IG=1,NBLO)
          WRITE(IU15)((THWOLD(J1,IG),J1=IJINF,IJSUP),IG=1,NBLO)
          WRITE(IU15)((USOLD(J1,IG),J1=IJINF,IJSUP),IG=1,NBLO)
          WRITE(IU15)((TAUW(J1,IG),J1=IJINF,IJSUP),IG=1,NBLO)
          WRITE(IU15)((Z0OLD(J1,IG),J1=IJINF,IJSUP),IG=1,NBLO)
          WRITE(IU15)((ROAIRO(J1,IG),J1=IJINF,IJSUP),IG=1,NBLO)
          WRITE(IU15)((ZIDLOLD(J1,IG),J1=IJINF,IJSUP),IG=1,NBLO)

          CLOSE(IU15)

        ELSE

          REWIND IU15
          WRITE(IU15) CDTPRO,CDATEWO,CDAWIFL,CDATEFL
          WRITE(IU15) U10OLD
          WRITE(IU15) THWOLD
          WRITE(IU15) USOLD
          WRITE(IU15) TAUW
          WRITE(IU15) Z0OLD
          WRITE(IU15) ROAIRO 
          WRITE(IU15) ZIDLOLD 

        ENDIF
      ENDIF
      IF (ITEST.GE.2)
     &    WRITE(IU06,*) ' SUB. WRITESTRESS: OUTPUT DONE'

      RETURN
      END SUBROUTINE WRITESTRESS
