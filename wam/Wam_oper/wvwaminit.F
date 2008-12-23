      SUBROUTINE WVWAMINIT (LLCOUPLED,IULOG,NLON,NLAT,RSOUTW,RNORTW)

!****  *WVWAMINIT* - INITIALISES A FEW VARIABLE AND READ WAVE MODEL
!                    CONFIGURATION PARAMETERS AND TABLES.


!     J. BIDLOT    ECMWF   FEBRUARY 2008: MOVE THIS PART OUT OF MPDECOMP 

!*    INTERFACE.
!     ----------
!     CALL *WVWAMINIT* (LLCOUPLED,IULOG,NLON,NLAT)

!      *LLCOUPLED*  LOGICAL  TRUE IF ATMOSPHERIC COUPLED RUN.
!      *IULOG*      INTEGER  STANDARD OUTPUT UNIT USED BY IFS.
!      *NLON*       INTEGER  MAXIMUM NUMBER OF LONGITUDES USED BY WAM. 
!      *NLAT*       INTEGER  MAXIMUM NUMBER OF LATITUDES USED BY WAM. 

! ----------------------------------------------------------------------

      USE YOWMAP   , ONLY : AMOSOP   ,AMONOP
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPRECR   ,NPRECI   ,
     &            MBXSIZE  ,KTAG
      USE YOWPARAM , ONLY : NGX      ,NGY      ,KWAMVER 
      USE YOWTEST  , ONLY : IU06
      USE YOWUNIT  , ONLY : IU07     ,IU08

      USE MPL_MODULE
      USE MPL_DATA_MODULE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: IULOG, NLON, NLAT
      INTEGER :: IREAD, LFILE
      INTEGER :: I4(2)
      INTEGER :: I_GET_UNIT

      REAL :: RSOUTW, RNORTW
      REAL :: X4(2)
      REAL :: ZHOOK_HANDLE

      CHARACTER*9 :: CMBXSIZE
      CHARACTER*80 :: FILENAME
      CHARACTER*80 :: LOGFILENAME

      LOGICAL :: LLCOUPLED
      LOGICAL :: LLEXIST

      IF (LHOOK) CALL DR_HOOK('WVWAMINIT',0,ZHOOK_HANDLE)

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

      WRITE(IU06,*) ' WAM SOFTWARE VERSION: ', KWAMVER
      CALL FLUSH (IU06)

!     CONSTANT FILES INPUT UNIT
!     -------------------------
      FILENAME='wam_grid'
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
     &  ' IS MISSING !!!!'
        WRITE(*,*) '*  WAVE MODEL INPUT FILE ',FILENAME(1:LFILE),
     &  ' IS MISSING !!!!'
        WRITE(IU06,*) '*                                  *'
       WRITE(IU06,*) '************************************'
        CALL ABORT1
      ENDIF
      IU07 = I_GET_UNIT(IU06, FILENAME(1:LFILE) , 'r', 'u', 0)

      FILENAME='wam_subgrid'
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
     &  ' IS MISSING !!!!'
        WRITE(*,*) '*  WAVE MODEL INPUT FILE ',FILENAME(1:LFILE),
     &  ' IS MISSING !!!!'
        WRITE(IU06,*) '*                                  *'
        WRITE(IU06,*) '************************************'
        CALL ABORT1
      ENDIF
      IU08 = I_GET_UNIT(IU06, FILENAME(1:LFILE) , 'r', 'u', 0)


!     DETERMINE BYTE STORAGE REPRESENTATION OF REAL NUMBERS
!     -----------------------------------------------------

      X4=1.
      NPRECR = KIND(X4)
      I4=1
      NPRECI = KIND(I4)

!     GET MESSAGE PASSING BUFFER SIZE
      CALL GETENV('MBX_SIZE',CMBXSIZE)
      READ(CMBXSIZE,'(I9)') MBXSIZE

!     READ GRID INPUT FROM PREPROC 
!     ----------------------------

      IREAD=1

      KTAG=1

      IF(IRANK.EQ.IREAD) THEN
        CALL GSTATS(1771,0)
        CALL READPRE (IU07)
        CLOSE (UNIT=IU07)
        CALL GSTATS(1771,1)
      ENDIF

!     SEND INFORMATION FROM READPRE TO ALL PE's
      CALL GSTATS(694,0)
      CALL MPBCASTGRID(IU06,IREAD,KTAG)
      CALL GSTATS(694,1)

      WRITE(IU06,*) ' WVWAMINT: WAVE MODEL CONFIGURATION READ IN'
      CALL FLUSH (IU06)

!     RETURN NECESSARY PARAMETERS
      NLON=NGX
      NLAT=NGY
      RSOUTW=AMOSOP
      RNORTW=AMONOP

      IF (LHOOK) CALL DR_HOOK('WVWAMINIT',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE WVWAMINIT 
