      SUBROUTINE WRITESTRESS(NREAL, RFIELD, NINTGR, IFIELD,
     &                       FILENAME,LPBOPEN,LPBCLOSE) 

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      APRIL 1997

!*    PURPOSE.
!     --------
!     WRITES THE RESTART WIND AND STRESS FIELDS TO FILE.

!**   INTERFACE.
!     ----------

!     CALL *WRITESTRESS*(NREAL, RFIELD, NINTGR, IFIELD,
!                        FILENAME,LPBOPEN,LPBCLOSE) 
!     *NREAL*    NUMBER OF FIELDS IN RFIELD
!     *RFIELD*     REAL RESTART FIELDS
!     *NINTGR*   NUMBER OF FIELDS IN IFIELD
!     *IFIELD*     INTEGER RESTART FIELDS
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
      USE YOWMPP   , ONLY : NPROC   ,NPRECR   ,NPRECI
      USE YOWPARAM , ONLY : NBLO    ,NIBLO    ,LL1D
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWSPEC  , ONLY : IJ2NEWIJ
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWUNIT  , ONLY : IU15
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL

! ----------------------------------------------------------------------

      INTEGER :: IFLD
      INTEGER :: NREAL, NINTGR
      CHARACTER*1 MODE
      CHARACTER*120 FILENAME
      LOGICAL LPBOPEN,LPBCLOSE
      INTEGER,DIMENSION(NIBLO,NBLO,NINTGR) :: IFIELD 
      REAL,DIMENSION(NIBLO,NBLO,NREAL) :: RFIELD 
      REAL, ALLOCATABLE :: RDUM(:,:)
      INTEGER, ALLOCATABLE :: IDUM(:,:)
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

        CALL PBWRITE(IUNIT,CDTPRO,14,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,CDATEWO,14,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,CDAWIFL,14,KRET)
        IF(KRET.LT.0) GOTO 1099
        CALL PBWRITE(IUNIT,CDATEFL,14,KRET)
        IF(KRET.LT.0) GOTO 1099

        KOUNT = NIBLO*NBLO*NPRECR
        IOUNT = NIBLO*NBLO*NPRECI
        IF(LL1D .OR. NPROC.EQ.1 ) THEN
          DO IFLD=1,NREAL
            CALL PBWRITE(IUNIT,RFIELD(:,:,IFLD),KOUNT,KRET)
            IF(KRET.LT.0) GOTO 1099
          ENDDO
          DO IFLD=1,NINTGR
            CALL PBWRITE(IUNIT,IFIELD(:,:,IFLD),IOUNT,KRET)
            IF(KRET.LT.0) GOTO 1099
          ENDDO
        ELSE
          ALLOCATE(RDUM(NIBLO,NBLO))
          DO IFLD=1,NREAL
            DO IG=1,NBLO
              DO IJ=1,NIBLO
                 RDUM(IJ,IG)=RFIELD(IJ2NEWIJ(IJ),IG,IFLD)
              ENDDO
            ENDDO
            CALL PBWRITE(IUNIT,RDUM,KOUNT,KRET)
            IF(KRET.LT.0) GOTO 1099
          ENDDO
          DEALLOCATE(RDUM)

          ALLOCATE(IDUM(NIBLO,NBLO))
          DO IFLD=1,NINTGR
            DO IG=1,NBLO
              DO IJ=1,NIBLO
                 IDUM(IJ,IG)=IFIELD(IJ2NEWIJ(IJ),IG,IFLD)
              ENDDO
            ENDDO
            CALL PBWRITE(IUNIT,IDUM,IOUNT,KRET)
            IF(KRET.LT.0) GOTO 1099
          ENDDO

          DEALLOCATE(IDUM)
        ENDIF

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

        OPEN(IU15,FILE=FILENAME,STATUS='UNKNOWN',FORM='UNFORMATTED')

        WRITE(IU15) CDTPRO,CDATEWO,CDAWIFL,CDATEFL
        IF(LL1D .OR. NPROC.EQ.1 ) THEN
          DO IFLD=1,NREAL
            WRITE(IU15)((RFIELD(J1,IG,IFLD),J1=1,NIBLO),IG=1,NBLO)
          ENDDO
          DO IFLD=1,NINTGR
            WRITE(IU15)((IFIELD(J1,IG,IFLD),J1=1,NIBLO),IG=1,NBLO)
          ENDDO
        ELSE
          DO IFLD=1,NREAL
            WRITE(IU15)
     &      ((RFIELD(IJ2NEWIJ(J1),IG,IFLD),J1=1,NIBLO),IG=1,NBLO)
          ENDDO
          DO IFLD=1,NINTGR
            WRITE(IU15)
     &      ((IFIELD(IJ2NEWIJ(J1),IG,IFLD),J1=1,NIBLO),IG=1,NBLO)
          ENDDO
        ENDIF

        IF(LMESSPASS) CLOSE(IU15)

      ENDIF
      IF (ITEST.GE.2)
     &    WRITE(IU06,*) ' SUB. WRITESTRESS: OUTPUT DONE'

      RETURN
      END SUBROUTINE WRITESTRESS
