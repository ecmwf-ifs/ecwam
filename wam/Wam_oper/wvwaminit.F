      SUBROUTINE WVWAMINIT (LLCOUPLED, IULOG, LLRNL,
     &                      NLON, NLAT, RSOUTW, RNORTW)

!****  *WVWAMINIT* - INITIALISES A FEW VARIABLE AND READ WAVE MODEL
!                    NAMELIST (if LLRNL is true) and
!                    CONFIGURATION PARAMETERS AND TABLES.


!     J. BIDLOT    ECMWF   FEBRUARY 2008: MOVE THIS PART OUT OF MPDECOMP 

!*    INTERFACE.
!     ----------
!     CALL *WVWAMINIT* (LLCOUPLED,IULOG,LLRNL,
!                       NLON, NLAT, RSOUTW, RNORTW)

!      *LLCOUPLED*  LOGICAL  TRUE IF ATMOSPHERIC COUPLED RUN.
!      *IULOG*      INTEGER  STANDARD OUTPUT UNIT USED BY IFS.
!      *LLRNL*      LOGICAL  IF TRUE READS WAM NAMELIST.
!      *NLON*       INTEGER  MAXIMUM NUMBER OF LONGITUDES USED BY WAM. 
!      *NLAT*       INTEGER  MAXIMUM NUMBER OF LATITUDES USED BY WAM. 
!      *RSOUTW*     REAL     MOST SOUTHERN LATITUDE IN WAM.
!      *RNORTW*     REAL     MOST NORTHERn LATITUDE IN WAM.

! ----------------------------------------------------------------------

      USE YOWMAP   , ONLY : AMOSOP   ,AMONOP
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPRECR   ,NPRECI   ,
     &            KTAG
      USE YOWPARAM , ONLY : NGX      ,NGY      ,KWAMVER 
      
      USE YOWTEST  , ONLY : IU06
      USE YOWUNIT  , ONLY : IREADG   ,NPROPAGS ,IU07     ,IU08     ,
     &            LWVWAMINIT
      USE YOWSTAT  , ONLY : IPROPAGS
      USE MPL_MODULE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWUNPOOL, ONLY : LLUNSTR

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IULOG
      INTEGER, INTENT(OUT) :: NLON, NLAT
      INTEGER :: IREAD, LFILE
      INTEGER :: I4(2)
      INTEGER :: I_GET_UNIT

      REAL, INTENT(OUT) :: RSOUTW, RNORTW
      REAL :: X4(2)
      REAL :: ZHOOK_HANDLE

      CHARACTER(LEN=1) :: C1 
      CHARACTER(LEN=80) :: FILENAME
      CHARACTER(LEN=80) :: LOGFILENAME

      LOGICAL, INTENT(IN) :: LLRNL
      LOGICAL, INTENT(IN) :: LLCOUPLED
      LOGICAL :: LLEXIST
      LOGICAL :: LFRST

      DATA LFRST /.TRUE./

! ----------------------------------------------------------------------

#ifdef ECMWF 
      IF (LHOOK) CALL DR_HOOK('WVWAMINIT',0,ZHOOK_HANDLE)
#endif

      LWVWAMINIT=.FALSE.

      IRANK = MPL_MYRANK()
      NPROC = MPL_NPROC()

!     STANDARD OUTPUT UNIT
!     --------------------

      IF(LLCOUPLED) THEN
        IU06=IULOG
      ELSE
        IU06=65
        LLEXIST=.TRUE.
        DO WHILE (LLEXIST)
          IU06=IU06+1
          LOGFILENAME='logfile.%p'
          CALL EXPAND_STRING(IRANK,NPROC,0,0,LOGFILENAME,1)
          INQUIRE(UNIT=IU06, OPENED=LLEXIST)
          IF(.NOT. LLEXIST) THEN
            OPEN(IU06,FILE=LOGFILENAME,STATUS='UNKNOWN')
            WRITE(IU06,*) ' STDOUT OF PE ', IRANK
            WRITE(IU06,*) ' ============ '
            WRITE(IU06,*) ' '
          ENDIF
        ENDDO
      ENDIF

!     GET CONTROLLING FLAG FROM INPUT NAMELIST 
!     ----------------------------------------
      IF(LFRST .AND. LLRNL) THEN
        CALL MPUSERIN
        LFRST=.FALSE.
      ENDIF
      IF(.NOT.LLRNL) THEN
        IPROPAGS=0
      ENDIF

      WRITE(IU06,*) ' WAM SOFTWARE VERSION: ', KWAMVER

!     DETERMINE BYTE STORAGE REPRESENTATION OF REAL NUMBERS
!     -----------------------------------------------------

      X4=1.
      NPRECR = KIND(X4)
      I4=1
      NPRECI = KIND(I4)

!     READ GRID INPUT FROM PREPROC 
!     ----------------------------

      IREAD=IREADG

!     CONSTANT FILES INPUT UNIT
!     -------------------------
      IF(IRANK.EQ.IREAD) THEN
        FILENAME='wam_grid_tables'
        LFILE=0
        LLEXIST=.FALSE.
        IF (FILENAME.NE. ' ') LFILE=LEN_TRIM(FILENAME)
        INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
        IF(.NOT. LLEXIST) THEN
          WRITE(IU06,*) '************************************'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '*  FATAL ERROR IN SUB. WVWAMINIT   *'
          WRITE(IU06,*) '*  =============================   *'
          WRITE(IU06,*) '*  WAVE MODEL INPUT FILE ',FILENAME(1:LFILE),
     &    ' IS MISSING !!!!'
          WRITE(*,*) '*  WAVE MODEL INPUT FILE ',FILENAME(1:LFILE),
     &    ' IS MISSING !!!!'
          WRITE(IU06,*) '*                                  *'
         WRITE(IU06,*) '************************************'
          CALL ABORT1
        ENDIF
        IU07 = I_GET_UNIT(IU06, FILENAME(1:LFILE) , 'r', 'u', 0)

        IF(IPROPAGS.LT.0 .OR. IPROPAGS.GT.NPROPAGS) THEN
          WRITE(IU06,*) '************************************'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '*  FATAL ERROR IN SUB. WVWAMINIT   *'
          WRITE(IU06,*) '*  WRONG VALUE FOR IPROPAGS:       *'
          WRITE(IU06,*) '*  IPROPAGS= ',IPROPAGS
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '************************************'
          CALL ABORT1
        ENDIF

        WRITE(C1,'(I1)') IPROPAGS 
        FILENAME='wam_subgrid_'//C1

        LFILE=0
        LLEXIST=.FALSE.
        IF (FILENAME.NE. ' ') LFILE=LEN_TRIM(FILENAME)
        INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
        IF(.NOT. LLEXIST) THEN
          WRITE(IU06,*) '************************************'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '*  FATAL ERROR IN SUB. WVWAMINIT   *'
          WRITE(IU06,*) '*  =============================   *'
          WRITE(IU06,*) '*  WAVE MODEL INPUT FILE ',FILENAME(1:LFILE),
     &    ' IS MISSING !!!!'
          WRITE(*,*) '*  WAVE MODEL INPUT FILE ',FILENAME(1:LFILE),
     &    ' IS MISSING !!!!'
          WRITE(IU06,*) '*                                  *'
          WRITE(IU06,*) '************************************'
          CALL ABORT1
        ENDIF
        IU08(IPROPAGS) = I_GET_UNIT(IU06, FILENAME(1:LFILE) ,'r','u',0)
      ENDIF

      KTAG=1

      CALL MPL_BARRIER(CDSTRING='WVWAMINIT:')

      CALL READPRE (IU07)
      CALL MPL_BARRIER(CDSTRING='WVWAMINIT:')

      WRITE(IU06,*) ' WVWAMINT: WAVE MODEL CONFIGURATION READ IN'
      CALL FLUSH (IU06)

!     RETURN NECESSARY PARAMETERS
      NLON=NGX
      NLAT=NGY
      RSOUTW=AMOSOP
      RNORTW=AMONOP

#ifdef ECMWF 
      IF (LHOOK) CALL DR_HOOK('WVWAMINIT',1,ZHOOK_HANDLE)
#endif

      RETURN
      END SUBROUTINE WVWAMINIT
