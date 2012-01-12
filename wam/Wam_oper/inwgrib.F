      SUBROUTINE INWGRIB (FILNM, IREAD,
     &                    CDATE, IPARAM, KZLEV, FIELD)

! -----------------------------------------------------------------     

!***  *INWGRIB* - READS IN AND UNPACKS WAVE MODEL GRIB FIELDS
!                 IT WILL ALSO INTERPOLATE THE DATA TO THE MODEL
!                 GRID (IF NECESSARY).
!                 (Except for wave spectra. See *GETSPEC*)
!                 MESSAGE PASSING WILL BE USE TO DISTRIBUTED THE
!                 CODED DATA TO ALL PE'S BEFORE DECODING.

!      J. BIDLOT    ECMWF    APRIL 2010. 

!**   INTERFACE.                                                        
!     ----------                                                        

!      *CALL INWGRIB* (FILNM, IREAD,
!    &                 CDATE, IPARAM, KZLEV, FIELD)

!        *FILNM*  - DATA INPUT FILENAME. 
!        *IREAD*  - ACCESS TO FILE ONLY FOR PE=IREAD.
!        *CDATE*  - DATE/TIME OF THE DATA READ.         
!        *IPARAM* - PARAMETER NUMBER. 
!        *KZLEV*  - REFERENCE LEVEL IN FULL METER
!                   SHOULD BE 0 EXCEPT FOR 233, 245 AND 249 WHERE IT
!                   MIGHT BE DIFFERENT THAN ZERO. 
!        *FIELD*  - UNPACKED DATA.

!     EXTERNALS.                                                        
!     ----------                                                        

!     GRIBAPI 

! ----------------------------------------------------------------------

      USE YOWGRID  , ONLY : NLONRGG
      USE YOWGRIBHD, ONLY : PPEPS    ,PPREC
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NIBLO
      USE YOWMAP   , ONLY : IRGG     ,AMOSOP   ,
     &            AMONOP   ,AMOWEP   ,AMOEAP   ,XDELLA   ,ZDELLO
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPRECI 
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : NWAM_BLKS
      USE YOWTEST  , ONLY : IU06     ,ITEST

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE MPL_MODULE
      USE GRIB_API_INTERFACE

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: NBIT

      INTEGER :: IREAD, IPARAM, KZLEV, IFORP
      INTEGER :: LFILE, KFILE_HANDLE, KGRIB_HANDLE
      INTEGER :: IRET, KBYTES, ISIZE, IUNIT
      INTEGER :: KK, MM
      INTEGER :: IBUF(2) 
      INTEGER, ALLOCATABLE :: INGRIB(:)

      REAL :: ZHOOK_HANDLE
      REAL :: FIELD(NGX,NGY)                                          

      CHARACTER(LEN=14) :: CDATE
      CHARACTER(LEN=24) :: FILNM

      LOGICAL :: LLEXIST
      LOGICAL :: LLRESIZING=.FALSE.

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('INWGRIB',0,ZHOOK_HANDLE)

      NBIT=NIBLO

      LLEXIST=.FALSE.
      LFILE = LEN_TRIM(FILNM)
      INQUIRE(FILE=FILNM(1:LFILE),EXIST=LLEXIST)
      IF(.NOT. LLEXIST) THEN
        WRITE (IU06,*) '*************************************'
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
        WRITE (IU06,*) '*  IN INWGRIB:                      *'
        WRITE (IU06,*) '*  COULD NOT FIND FILE ',FILNM
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*************************************'
        CALL ABORT1
      ENDIF

!     READ DATA ON PE IREAD
      IF(IRANK.EQ.IREAD) THEN
1121    CALL IGRIB_OPEN_FILE(KFILE_HANDLE,FILNM(1:LFILE),'r')

1021    ISIZE=NBIT
        KBYTES=ISIZE*NPRECI
        IF(.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))
          CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE,INGRIB,KBYTES,IRET)
        IF(IRET.EQ.JPGRIB_BUFR_TOO_SHORT) THEN
          CALL KGRIBSIZE(IU06, IUNIT, KBYTES, NBIT, 'INWGRIB')
          DEALLOCATE(INGRIB)
          LLRESIZING=.TRUE.
          GOTO 1021
        ELSEIF(LLRESIZING .AND. IRET.NE.JPGRIB_END_OF_FILE) THEN
!         LOOP UNTIL YOU HAVE EXPLORE THE SIZE FOR THE WHOLE FILE.
          DEALLOCATE(INGRIB)
          GOTO 1021
        ELSEIF(LLRESIZING .AND. IRET.EQ.JPGRIB_END_OF_FILE) THEN
!         WE SHOULD HAVE THE MAXIMUM SIZE NECESSARY, START ALL OVER.
          DEALLOCATE(INGRIB)
          LLRESIZING=.FALSE.
          CALL IGRIB_CLOSE_FILE(KFILE_HANDLE)
          GOTO 1121
        ELSEIF(IRET.EQ.JPGRIB_END_OF_FILE) THEN
          WRITE(IU06,*) '**********************************'
          WRITE(IU06,*) '* INWGRIB: END OF FILE ENCOUNTED'
          WRITE(IU06,*) '**********************************'
          CALL ABORT1
          ELSEIF(IRET.NE.JPGRIB_SUCCESS) THEN
          WRITE(IU06,*) '**********************************'
          WRITE(IU06,*) '* INWGRIB: FILE HANDLING ERROR'
          WRITE(IU06,*) '**********************************'
          CALL ABORT1
        ENDIF
      ENDIF

      WRITE(IU06,*) ' SUB. INWGRIB - READ FROM ',FILNM

      CALL MPL_BARRIER(CDSTRING='INWGRIB: DATA READ IN')

!     SEND GRIB DATA TO THE OTHER PE'S
      IF(NPROC.GT.1) THEN
        CALL GSTATS(619,0)
        IF(IRANK.EQ.IREAD) THEN
          IBUF(1)=ISIZE
          IBUF(2)=KBYTES
        ENDIF
        CALL MPL_BROADCAST(IBUF(1:2),KROOT=IREAD,KTAG=1,
     &                     CDSTRING='INWGRIB IBUF:')
        IF(IRANK.NE.IREAD) THEN
          ISIZE=IBUF(1)
          KBYTES=IBUF(2)
          ALLOCATE(INGRIB(ISIZE))
        ENDIF

        CALL MPL_BROADCAST(INGRIB(1:ISIZE),KROOT=IREAD,KTAG=2,
     &                     CDSTRING='INWGRIB: INGRIB')
        CALL GSTATS(619,1)
      ENDIF

!     DECODE THE GRIB DATA
!     (and interpolate to model grid if necessary)

      KGRIB_HANDLE=-99
      CALL IGRIB_NEW_FROM_MESSAGE(KGRIB_HANDLE,INGRIB)

      CALL GRIB2WGRID (IU06, ITEST, NWAM_BLKS, 
     &                 KGRIB_HANDLE, INGRIB, ISIZE,
     &                 NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                 AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                 ZMISS, PPREC, PPEPS,
     &                 CDATE, IFORP, IPARAM, KZLEV, KK, MM, FIELD)

      CALL IGRIB_RELEASE(KGRIB_HANDLE)

      IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)

      IF (LHOOK) CALL DR_HOOK('INWGRIB',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE INWGRIB
