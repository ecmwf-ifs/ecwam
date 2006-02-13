! ----------------------------------------------------------------------

      SUBROUTINE FILE_TRANSFER (IU06, IUNIT, FILENA, CDOPT) 

! ----------------------------------------------------------------------

!**** *FILE_TRANSFER* - COPY FORTRAN UNIT FILE TO A FILE 

!     J. BIDLOT     ECMWF JUNE  1996  MESSAGE PASSING
!     B. HANSEN     ECMWF APRIL 1996  MOVE CYMPPATH TO COMTEXT AND ADD 
!                                     CDOPT.

!*    PURPOSE.
!     --------

!     COPY FORTRAN UNIT FILE TO A FILE 

!**   INTERFACE.
!     ----------

!       *CALL* *FILE_TRANSFER*(IU06, IUNIT, FILENA) 

!         *IU06*   INTEGER    UNIT FOR PRINTER MESSAGES.
!         *IUNIT*  INTEGER    FORTRAN UNIT.
!         *FILENA* CHARACTER  DESTINATION FILENAME
!         *CDOPT*  CHARACTER  OPTION
!                             = 'S' SAVE FILE CONNECTED TO IUNIT
!                             = 'G' GET FILE AND CONNECT FILE TO IUNIT

!     EXTERNALS.
!     ----------

!       *PBOPEN*
!       *PBREAD2*
!       *PBWRITE*
!       *PBCLOSE* 
!       *EXPRAND_STRING - ADD EXTENSION TO FILENAME FUNCTION OF THE PE 

!     METHOD.
!     -------

!     REFERENCES.
!      -----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,
     &            NPRECR
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH


      INTEGER IU06, IUNIT
      CHARACTER*  1 MODE,
     &              CDOPT
      CHARACTER* 17 FILENA
      CHARACTER* 160 FILENAM,
     &              FILENAMEF, 
     &              FILENAMET 

      INTEGER, ALLOCATABLE :: KARRAY(:) 
      LOGICAL LMORE

! ----------------------------------------------------------------------
      LMORE=.TRUE.

      IF (IUNIT.LT.10) THEN
        WRITE(FILENAMEF,'(A,I1)') 'fort.',IUNIT
      ELSE
        WRITE(FILENAMEF,'(A,I2)') 'fort.',IUNIT
      ENDIF
      FILENAM=FILENAMEF
      LNAME = LEN_TRIM(FILENAM)
      IF(ICPLEN.GT.0) THEN
        FILENAMEF=CPATH(1:ICPLEN)//'/'//FILENAM(1:LNAME)
      ENDIF

      FILENAM=FILENA
      LNAME = LEN_TRIM(FILENAM)
      IF(ICPLEN.EQ.0) THEN
        FILENAMET=FILENAM(1:LNAME)
      ELSE
        FILENAMET=CPATH(1:ICPLEN)//'/'//FILENAM(1:LNAME)
      ENDIF

      IF (CDOPT.EQ.'G') THEN
        FILENAM = FILENAMET
        FILENAMET = FILENAMEF
        FILENAMEF = FILENAM
      ENDIF

      WRITE(IU06,*) 'FILE_TRANSFER: COPY FILE ',
     & FILENAMEF(1:LEN_TRIM(FILENAMEF)), ' TO ',
     & FILENAMET(1:LEN_TRIM(FILENAMET))

      KLENG=((NSUP-NINF+1)*NANG*NFRE*NPRECR)/9+8+1
      ALLOCATE(KARRAY(KLENG)) 

!       1. OPEN INPUT AND OUTPUT FILES
!          ---------------------------
      MODE='r'
      CALL PBOPEN(KUNITF,FILENAMEF,MODE,KRET)
      IF(KRET.LT.0) THEN
        WRITE (IU06,*) '******************************************'
        WRITE (IU06,*) '*                                        *'
        WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN       *'
        WRITE (IU06,*) '*   IN FILE_TRANSFER                     *'
        IF(KRET.EQ.-1)
     &   WRITE (IU06,*) ' COULD NOT OPEN FILE ',FILENAMEF
        IF(KRET.EQ.-2) WRITE (IU06,*) ' INVALID FILENAME ',FILENAMEF
        IF(KRET.EQ.-3) WRITE (IU06,*) ' INVALID OPEN MODE SPECIFIED'
        WRITE (IU06,*) '*                                        *'
        WRITE (IU06,*) '******************************************'
        CALL ABORT1
      ENDIF

      MODE='w'
      CALL PBOPEN(KUNITO,FILENAMET,MODE,KRET)
      IF(KRET.LT.0) THEN
        WRITE (IU06,*) '******************************************'
        WRITE (IU06,*) '*                                        *'
        WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN       *'
        WRITE (IU06,*) '*   IN FILE_TRANSFER                     *'
        IF(KRET.EQ.-1)
     &   WRITE (IU06,*) 'COULD NOT OPEN FILE ',FILENAMET
        IF(KRET.EQ.-2) WRITE (IU06,*) 'INVALID FILENAME ',FILENAMET
        IF(KRET.EQ.-3) WRITE (IU06,*) 'INVALID OPEN MODE SPECIFIED'
        WRITE (IU06,*) '*                                        *'
        WRITE (IU06,*) '******************************************'
        CALL ABORT1
      ENDIF

!       2. COPY FILES
!          ----------
      DO WHILE (LMORE)

        CALL PBREAD2(KUNITF,KARRAY,KLENG,KRET)
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBREAD2      *'
          WRITE (IU06,*) '*   IN FILE_TRANSFER                     *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          CALL ABORT1
        ENDIF

        IF(KRET.LT.KLENG) THEN 
          KOUNT=KRET 
          LMORE=.FALSE.
        ELSE
          KOUNT=KLENG 
          LMORE=.TRUE.
        ENDIF

        CALL PBWRITE(KUNITO,KARRAY,KOUNT,KRET)
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '*****************************************'
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*    ERROR FOLLOWING CALL TO PBWRITE    *'
          WRITE (IU06,*) '*    IN FILE_TRANSFER                   *'
          WRITE (IU06,*) '*      FILE ',FILENAMET
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*****************************************'
          CALL ABORT1
        ENDIF

      ENDDO

!       3. CLOSE FILES
!          -----------
      CALL PBCLOSE(KUNITF,KRET) 
      IF(KRET.LT.0) THEN
        WRITE (IU06,*) '************************************'
        WRITE (IU06,*) '*                                  *'
        WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
        WRITE (IU06,*) '* IN FILE_TRANSFER                 *'
        WRITE (IU06,*) '* FILE ',FILENAMEF
        WRITE (IU06,*) '*                                  *'
        WRITE (IU06,*) '************************************'
        CALL ABORT1
      ENDIF

      CALL PBCLOSE(KUNITO,KRET) 
      IF(KRET.LT.0) THEN
        WRITE (IU06,*) '*******************************************'
        WRITE (IU06,*) '*                                         *'
        WRITE (IU06,*) '*     ERROR FOLLOWING CALL TO PBCLOSE     *'
        WRITE (IU06,*) '*     IN FILE_TRANSFER                    *'
        WRITE (IU06,*) '*     FILE ',FILENAMET
        WRITE (IU06,*) '*                                         *'
        WRITE (IU06,*) '*******************************************'
        CALL ABORT1
      ENDIF

      DEALLOCATE(KARRAY)

      RETURN
      END SUBROUTINE FILE_TRANSFER
