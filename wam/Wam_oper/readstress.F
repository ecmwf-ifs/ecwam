      SUBROUTINE READSTRESS(NREAL, RFIELD, FILENAME)

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      MARCH 1997

!*    PURPOSE.
!     --------
!     READS THE RESTART WIND AND STRESS FIELDS

!**   INTERFACE.
!     ----------
!     CALL *READSTRESS*(NREAL, RFIELD, FILENAME)  
!     *NREAL*    NUMBER OF FIELDS IN RFIELD
!     *RFIELD*     REAL RESTART FIELDS
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
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL

! ----------------------------------------------------------------------
!     MODIFICATIONS.
!     --------------
!     B. HANSEN    ECMWF      NOVEMBER 1998  KOUNT IS 12 FOR Y2K.

      CHARACTER*1 MODE
      CHARACTER*120 FILENAME

      INTEGER :: IFLD
      INTEGER :: NREAL

      REAL,DIMENSION(NIBLO,NBLO,NREAL) :: RFIELD 

      REAL, ALLOCATABLE :: RDUM(:,:)
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
          DO IFLD=1,NREAL
            CALL PBREAD(IUNIT,RFIELD(:,:,IFLD),KOUNT,KRET)
            IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099
          ENDDO
        ELSE
          ALLOCATE(RDUM(NIBLO,NBLO))
          DO IFLD=1,NREAL
            CALL PBREAD(IUNIT,RDUM,KOUNT,KRET)
            IF(KRET.LT.0.OR.KRET.NE.KOUNT) GOTO 1099
            DO IG=1,NBLO
              DO IJ=1,NIBLO
                RFIELD(IJ2NEWIJ(IJ),IG,IFLD)=RDUM(IJ,IG)
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(RDUM)

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

        LFILE=0
        IF (FILENAME.NE. ' ') LFILE=LEN_TRIM(FILENAME)
        IUNIT=I_GET_UNIT(IU06, FILENAME(1:LFILE) , 'r', 'u', 0)

        READ(IUNIT) CDTPRO,CDATEWO,CDAWIFL,CDATEFL
        IF(LL1D .OR. NPROC.EQ.1 ) THEN
          DO IFLD=1,NREAL
            READ(IUNIT) RFIELD(:,:,IFLD) 
          ENDDO
        ELSE
!       WHEN 2-D DECOMPOSITION IS USED THEN THE INDEXES IJ ARE
!       RE-LABELLED BUT THE BINARY INPUT FILES ARE IN THE OLD MAPPING
          DO IFLD=1,NREAL
            READ(IUNIT)
     &      ((RFIELD(IJ2NEWIJ(IJ),IG,IFLD),IJ=1,NIBLO),IG=1,NBLO)
          ENDDO
        ENDIF

        CLOSE(IUNIT)

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
