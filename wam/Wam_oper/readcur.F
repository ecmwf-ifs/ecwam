      SUBROUTINE READCUR (IU02, CDTCR)

! ----------------------------------------------------------------------

!**** *READCUR* - ROUTINE TO READ A CURRENT FIELD.

!     SUSANNE HASSELMANN        JUNE     1990.
!     H. GUNTHER    ECMWF/GKSS  DECEMBER 1990  MODIFIED FOR CYCLE_4.

!*    PURPOSE.
!     --------
!       INPUT OF A CURRENT FIELD.

!**   INTERFACE.
!     ----------
!       *CALL* *READCUR (IU02, CDTCR)*
!          *IU02*   - INPUT UNIT FOR CURRENTS.
!          *CDTCR*  - DATE/TIME OF THE DATA READ.
!     METHOD.
!     -------

!       FORMATED READ; CURRENT GRID IS LAT / LONG GRID.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.
!       *ADJUST*    - CORRECTS LONGITUDE INPUT.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCURG  , ONLY : DUCLO    ,DUCLA    ,UCLO     ,UCLH     ,
     &            UCLW     ,UCLE     ,UCUR     ,VCUR     ,KRCI     ,
     &            KCCI     ,IPCUR
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IU02
      CHARACTER(LEN=14), INTENT(OUT) :: CDTCR 


      INTEGER :: I, J
      CHARACTER(LEN=14), PARAMETER :: ZERO='              '

! ----------------------------------------------------------------------

!*    1. READ CURRENT FILE HEADER.
!        -------------------------


      READ (IU02, '(6F10.5,2I3)',END=3001,ERR=3001)
     &    UCLO, UCLH, UCLW, UCLE, DUCLO, DUCLA, KRCI, KCCI
      IF (ITEST.GT.0) THEN
        WRITE (IU06,*) '  '
        WRITE (IU06,*) ' SUB. READCUR: FILE HEADER READ'
        WRITE (IU06,*) ' GRID SPECIFICATION ARE:'
        WRITE (IU06,*) ' NUMBER OF COLUMNS IN GRID KCCI = ', KCCI
        WRITE (IU06,*) ' NUMBER OF ROWS    IN GRID KRCI = ', KRCI
        WRITE (IU06,*) ' SOUTHERN MOST LATITUDE    UCLO = ', UCLO
        WRITE (IU06,*) ' NORTHERN MOST LATITUDE    UCLH = ', UCLH
        WRITE (IU06,*) ' WESTERN MOST LONGITUDE    UCLW = ', UCLW
        WRITE (IU06,*) ' EASTERN MOST LONGITUDE    UCLE = ', UCLE
        WRITE (IU06,*) ' LATITUDE  STEP IS        DUCLA = ', DUCLA
        WRITE (IU06,*) ' LONGITUDE STEP IS        DUCLO = ', DUCLO
      ENDIF

      ALLOCATE(UCUR(KCCI,KRCI))
      ALLOCATE(VCUR(KCCI,KRCI))

!*    1.2  CHECK DIMENSIONS
!     ----------------------

      IF (KRCI.LT.2 .OR. KCCI.LT.2) THEN
        WRITE (IU06,*) ' *****************************************'
        WRITE (IU06,*) ' *                                       *'
        WRITE (IU06,*) ' *     FATAL ERROR IN SUB. READCUR       *'
        WRITE (IU06,*) ' *     ===========================       *'
        WRITE (IU06,*) ' * LESS THAN 2 LATITUDES OR LONGITUDES   *'
        WRITE (IU06,*) ' * IN CURRENT INPUT.                     *'
        WRITE (IU06,*) ' *                                       *'
        WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
        WRITE (IU06,*) ' *                                       *'
        WRITE (IU06,*) ' *****************************************'
        CALL ABORT1
      ENDIF

      CALL ADJUST (UCLW, UCLE)
      IPCUR=0
      IF (MOD(UCLE+DUCLO-UCLW+720.,360.).LT.0.0000001) IPCUR=1

! ----------------------------------------------------------------------

!*    2. READ CURRENT FIELD.
!        -------------------

      CDTCR = ZERO
      DO J=1,KRCI
        READ(IU02,'(8E9.3)',END=3002,ERR=3002) (UCUR(I,J),I=1,KCCI)
      ENDDO

      DO J=1,KRCI
        READ(IU02,'(8E9.3)',END=3003,ERR=3003) (VCUR(I,J),I=1,KCCI)
      ENDDO

      RETURN

! ----------------------------------------------------------------------

!*    3. ERROR HANDLING.
!        ---------------

 3001 CONTINUE
      WRITE (IU06,*) ' *****************************************'
      WRITE (IU06,*) ' *                                       *'
      WRITE (IU06,*) ' *      FATAL ERROR IN SUB. READCUR      *'
      WRITE (IU06,*) ' *      ===========================      *'
      WRITE (IU06,*) ' * READ ERROR OR EOF ON CURRENT FILE.    *'
      WRITE (IU06,*) ' * FILE HEADER EXPECTED                  *'
      WRITE (IU06,*) ' *                                       *'
      WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
      WRITE (IU06,*) ' *                                       *'
      WRITE (IU06,*) ' *****************************************'
      CALL ABORT1
 3002 CONTINUE
      WRITE (IU06,*) ' *****************************************'
      WRITE (IU06,*) ' *                                       *'
      WRITE (IU06,*) ' *      FATAL ERROR IN SUB. READCUR      *'
      WRITE (IU06,*) ' *      ===========================      *'
      WRITE (IU06,*) ' * READ ERROR OR EOF ON CURRENT FILE.    *'
      WRITE (IU06,*) ' * U - COMPONETS EXPECTED                *'
      WRITE (IU06,*) ' *                                       *'
      WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
      WRITE (IU06,*) ' *                                       *'
      WRITE (IU06,*) ' *****************************************'
      CALL ABORT1
 3003 CONTINUE
      WRITE (IU06,*) ' *****************************************'
      WRITE (IU06,*) ' *                                       *'
      WRITE (IU06,*) ' *      FATAL ERROR IN SUB. READCUR      *'
      WRITE (IU06,*) ' *      ===========================      *'
      WRITE (IU06,*) ' * READ ERROR OR EOF ON CURRENT FILE.    *'
      WRITE (IU06,*) ' * V - COMPONETS EXPECTED                *'
      WRITE (IU06,*) ' *                                       *'
      WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
      WRITE (IU06,*) ' *                                       *'
      WRITE (IU06,*) ' *****************************************'
      CALL ABORT1

      END SUBROUTINE READCUR
