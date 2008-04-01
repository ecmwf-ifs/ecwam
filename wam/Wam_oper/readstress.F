      SUBROUTINE READSTRESS(RU10OLD,RTHWOLD,RUSOLD,RTAUW,RZ0OLD,
     &                      RADOLD,RZIDLOLD,KICEMASK,FILENAME)

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      MARCH 1997

!*    PURPOSE.
!     --------
!     READS THE RESTART WIND AND STRESS FIELDS

!**   INTERFACE.
!     ----------
!     CALL *READSTRESS*(RU10OLD,RTHWOLD,RUSOLD,RTAUW,RZ0OLD,
!                       RADOLD,RZIDLOLD,KICEMASK,FILENAME)
!     *RU10OLD*   WIND SPEED.
!     *RTHWOLD*   WIND DIRECTION (RADIANS).
!     *RUSOLD*    FRICTION VELOCITY.
!     *RTAUW*     WAVE STRESS.
!     *RZ0OLD*    ROUGHNESS LENGTH IN.
!     *RADOLD*    AIR/WATER DENSITY RATIO.
!     *RZIDLOLD*  Zi/L (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!     *KICEMASK*  SEA ICE MASK.
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
      USE YOWMPP   , ONLY : NPROC    ,NPRECR   ,NPRECI
      USE YOWPARAM , ONLY : NBLO     ,NIBLO    ,LL1D
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWSPEC  , ONLY : IJ2NEWIJ
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWUNIT  , ONLY : IU15
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL

! ----------------------------------------------------------------------
!     MODIFICATIONS.
!     --------------
!     B. HANSEN    ECMWF      NOVEMBER 1998  KOUNT IS 12 FOR Y2K.

      CHARACTER*1 MODE
      CHARACTER*120 FILENAME

      REAL,DIMENSION(NIBLO,NBLO) :: RU10OLD,RTHWOLD,RUSOLD,RZ0OLD,RTAUW,
     &                              RADOLD,RZIDLOLD
      REAL, ALLOCATABLE :: RDUM(:,:)
      INTEGER,DIMENSION(NIBLO,NBLO) :: KICEMASK 
      INTEGER, ALLOCATABLE :: IDUM(:,:)
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


        KOUNT = 14 
        CALL PBREAD(IUNIT,CDTPRO,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        CALL PBREAD(IUNIT,CDATEWO,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        CALL PBREAD(IUNIT,CDAWIFL,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        CALL PBREAD(IUNIT,CDATEFL,KOUNT,KRET)
        IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099

        KOUNT = NIBLO*NBLO*NPRECR
        IOUNT = NIBLO*NBLO*NPRECI
        IF(LL1D .OR. NPROC.EQ.1 ) THEN
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
          CALL PBREAD(IUNIT,RADOLD,KOUNT,KRET)
          IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099
          CALL PBREAD(IUNIT,RZIDLOLD,KOUNT,KRET)
          IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099
          CALL PBREAD(IUNIT,KICEMASK,IOUNT,KRET)
          IF(KRET.LT.0.OR.KRET.NE.IOUNT) GOTO 1098
        ELSE
          ALLOCATE(RDUM(NIBLO,NBLO))
          CALL PBREAD(IUNIT,RDUM,KOUNT,KRET)
          IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099
          DO IG=1,NBLO
            DO IJ=1,NIBLO
               RU10OLD(IJ2NEWIJ(IJ),IG)=RDUM(IJ,IG)
            ENDDO
          ENDDO
          CALL PBREAD(IUNIT,RDUM,KOUNT,KRET)
          IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099
          DO IG=1,NBLO
            DO IJ=1,NIBLO
               RTHWOLD(IJ2NEWIJ(IJ),IG)=RDUM(IJ,IG)
            ENDDO
          ENDDO
          CALL PBREAD(IUNIT,RDUM,KOUNT,KRET)
          IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099
          DO IG=1,NBLO
            DO IJ=1,NIBLO
               RUSOLD(IJ2NEWIJ(IJ),IG)=RDUM(IJ,IG)
            ENDDO
          ENDDO
          CALL PBREAD(IUNIT,RDUM,KOUNT,KRET)
          IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099
          DO IG=1,NBLO
            DO IJ=1,NIBLO
               RTAUW(IJ2NEWIJ(IJ),IG)=RDUM(IJ,IG)
            ENDDO
          ENDDO
          CALL PBREAD(IUNIT,RDUM,KOUNT,KRET)
          IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099
          DO IG=1,NBLO
            DO IJ=1,NIBLO
               RZ0OLD(IJ2NEWIJ(IJ),IG)=RDUM(IJ,IG)
            ENDDO
          ENDDO
          CALL PBREAD(IUNIT,RDUM,KOUNT,KRET)
          IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099
          DO IG=1,NBLO
            DO IJ=1,NIBLO
               RADOLD(IJ2NEWIJ(IJ),IG)=RDUM(IJ,IG)
            ENDDO
          ENDDO
          CALL PBREAD(IUNIT,RDUM,KOUNT,KRET)
          IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099
          DO IG=1,NBLO
            DO IJ=1,NIBLO
               RZIDLOLD(IJ2NEWIJ(IJ),IG)=RDUM(IJ,IG)
            ENDDO
          ENDDO

          DEALLOCATE(RDUM)

          ALLOCATE(IDUM(NIBLO,NBLO))
          CALL PBREAD(IUNIT,IDUM,IOUNT,KRET)
          IF(KRET.LT.0.OR.KRET.NE.IOUNT) GOTO 1098
          DO IG=1,NBLO
            DO IJ=1,NIBLO
               KICEMASK(IJ2NEWIJ(IJ),IG)=IDUM(IJ,IG)
            ENDDO
          ENDDO

          DEALLOCATE(IDUM)
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
        IF(LL1D .OR. NPROC.EQ.1 ) THEN
          READ(IU15) RU10OLD
          READ(IU15) RTHWOLD
          READ(IU15) RUSOLD
          READ(IU15) RTAUW
          READ(IU15) RZ0OLD
          READ(IU15) RADOLD
          READ(IU15) RZIDLOLD 
          READ(IU15) KICEMASK
        ELSE
!       WHEN 2-D DECOMPOSITION IS USED THEN THE INDEXES IJ ARE
!       RE-LABELLED BUT THE BINARY INPUT FILES ARE IN THE OLD MAPPING
          READ(IU15) ((RU10OLD(IJ2NEWIJ(IJ),IG),IJ=1,NIBLO),IG=1,NBLO)
          READ(IU15) ((RTHWOLD(IJ2NEWIJ(IJ),IG),IJ=1,NIBLO),IG=1,NBLO)
          READ(IU15) ((RUSOLD(IJ2NEWIJ(IJ),IG),IJ=1,NIBLO),IG=1,NBLO)
          READ(IU15) ((RTAUW(IJ2NEWIJ(IJ),IG),IJ=1,NIBLO),IG=1,NBLO)
          READ(IU15) ((RZ0OLD(IJ2NEWIJ(IJ),IG),IJ=1,NIBLO),IG=1,NBLO)
          READ(IU15) ((RADOLD(IJ2NEWIJ(IJ),IG),IJ=1,NIBLO),IG=1,NBLO)
          READ(IU15) ((RZIDLOLD(IJ2NEWIJ(IJ),IG),IJ=1,NIBLO),IG=1,NBLO)
          READ(IU15) ((KICEMASK(IJ2NEWIJ(IJ),IG),IJ=1,NIBLO),IG=1,NBLO)
        ENDIF

        CLOSE(IU15)

      ENDIF

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) ' SUB. READSTRESS: INPUT OF RESTART FILES DONE'
        CALL FLUSH(IU06)
      ENDIF

      RETURN

1098  CONTINUE
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
      RETURN

1099  CONTINUE
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

      RETURN

      END SUBROUTINE READSTRESS
