      SUBROUTINE READWIND (CDTWIR, FILNM, LLNOTOPENED, IREAD)

!***  *READWIND* - PROGRAM TO FORMAT WIND FIELDS FROM MARS 
!                  FROM INPUT FILE sfcwind.

!     PETER JANSSEN      ECMWF        FEBRUARY 1987 


!    MODIFIED BY:                                                 
!    ------------                                                       

!    LIANA ZAMBRESKY     GKSS/ECMWF   JULY 1988

!        1. PASSING OF ARRAYS AS FORMAL PARMETERS TO SAVE MEMORY
!        2. PROCESSING WINDS FOR ONLY ONE WIND TIME STEP

!    PEDRO VITERBO       ECMWF        OCTOBER 1988

!        DEALING WITH MARS PACKED DATA.

!    H. GUNTHER      ECMWF/GKSS       OCTOBER 1991

!        GRIB EDITION 1.

!     J. DOYLE    NRL/ECMWF     OCT. 1996
!                 DECODES GRID WIND FIELDS WHICH ARE PASSED AS ARRAYS
!                 WHEN MODEL IS COUPLED TO ATMOPHERIC MODEL

!    B. HANSEN    ECMWF 1997
!                 RENAME SUBROUTINE TO MATCH THE NAME OF THE FILE.
!                 USE INCLUDE FILES FOR ALL COMMON BLOCKS.

!    J. BIDLOT    ECMWF 1998 RENAMED READWND_MARS AS READWIND

!    S. ABDALLA   ECMWF 2001
!                 GETS ATMOSPHERIC PARAMETERS FOR GUSTINESS
!                 AND AIR DENSITY COMPUTATIONS
!
!    P.A.E.M. JANSSEN   ECMWF 2003
!                 USE NEUTRAL WINDS RATHER THAN REAL WINDS
!                 TRANSFER AIR DENSITY CALCULATION TO IFS

!    J. BIDLOT    ECMWF 2008  REMOVE PART FOR ATMOSPHERIC FIELDS
!                             TO *IFSTOWAM*

!    J. BIDLOR    ECMWF 2010 USE GRIB API

!     PURPOSE                                                       
!     -------                                                      

!     *READWIND*
!         READS TARGET FILES FROM THE MARS ARCHIVE AND INFERS FROM THIS 
!             A) THE DATE OF THE WIND FIELD
!             B) THE DEFINITION OF THE MARS GRID
!             C) THE U AND V COMPONENT OF THE WIND FIELD
!             D) THE TYPE OF WIND FIELD AND
!             E) THE UNPACKED WINDS

!     INTERFACE
!     ---------

!     *CALL* *READWIND (CDTWIR, FILNM, LLNOTOPENED, IREAD)*
!
!        *CDTWIR* - DATE/TIME OF THE DATA READ
!        *FILNM*    FILENAME OF INPUT FILE
!        *LLNOTOPENED*  TRUE IF THE INPUT FILE HAS TO BE OPENED
!        *IREAD*    PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!                   (IF NEEDED)


!     EXTERNALS
!     ---------

!     *ABORT1*         TERMINATES PROCESSING
!     *INCDATE*        INCREMENTS DATE
!     *KGRIBSIZE*
!     *FLUSH*
!     *GRB2WGRD*       TRANSFORM GRIB FIELDS TO WAM GRID FIELDS

! --------------------------------------------------------------------- 

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWGRID  , ONLY : IGL      ,NLONRGG
      USE YOWICE   , ONLY : LICERUN  ,IPARAMCI ,LICETH
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,IXLG     , 
     &            KXLT
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPRECI
      USE YOWPARAM , ONLY : NIBLO    ,CLDOMAIN ,
     &            SWAMPWIND,SWAMPWIND2,DTNEWWIND,LTURN90 ,LWDINTS  ,
     &            SWAMPCIFR
      USE YOWSTAT  , ONLY : CDATEA   ,IDELWI   ,LADEN    ,LGUST    ,
     &            NPROMA_WAM
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWNDG  , ONLY : DLAM     ,DPHI     ,RLATS    ,RLATN    ,
     &            RLONL    ,RLONR    ,KCOL     ,KROW     ,ICODE    ,
     &            IWPER    ,ICOORD
      USE YOWWIND  , ONLY : NXFF     ,NYFF     ,FIELDG   ,IUNITW   ,
     &            NBITW    ,CWDFILE  ,LLWSWAVE ,LLWDWAVE
      USE YOWWIND  , ONLY : FIELDG_coupl
      USE YOWPCONS , ONLY : RAD      ,ZMISS    ,ROAIR    ,WSTAR0
      USE YOWPD, ONLY : MNP => npa
      USE YOWUNPOOL ,ONLY : LLUNSTR
      USE MPL_MODULE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE GRIB_API_INTERFACE

! --------------------------------------------------------------------- 

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IREAD 
      CHARACTER(LEN=14), INTENT(INOUT) :: CDTWIR
      CHARACTER(LEN=24), INTENT(INOUT) :: FILNM
      LOGICAL, INTENT(INOUT) :: LLNOTOPENED

      INTEGER :: NFLD  
      INTEGER :: I, J, IVAR, JSN
      INTEGER :: ISIZE
      INTEGER :: IFORP, IPARAM, KZLEV, IDM 
      INTEGER :: IWTIME, IDTTURN
      INTEGER :: LNAME 
      INTEGER :: I_GET_UNIT
      INTEGER :: IRET
      INTEGER :: KGRIB_HANDLE
      INTEGER :: IDUM(2)
      INTEGER, ALLOCATABLE :: KGRIB(:)
      INTEGER :: NLONRGG_LOC(NYFF)
      INTEGER(KIND=JPKSIZE_T) :: KBYTES

      REAL :: ZDUM
      REAL :: UWIND, VWIND, WSPEED, WTHETA
      REAL :: ZHOOK_HANDLE
      REAL, ALLOCATABLE :: WORK(:,:)  

      CHARACTER(LEN=14) :: CDTTURN
      CHARACTER(LEN=14), SAVE :: CWDDATE

      LOGICAL :: LLABORT
      LOGICAL :: LLEXIST
      LOGICAL, ALLOCATABLE :: LLNOTREAD(:)
! --------------------------------------------------------------------  

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('READWIND',0,ZHOOK_HANDLE)
#endif

      IF (ITEST.GT.1) THEN
        WRITE(IU06,*) ' SUB. READWIND - STARTING ' 
        CALL FLUSH(IU06)
      ENDIF

      IF(LLUNSTR) THEN
        NLONRGG_LOC(:)=MNP
      ELSE
        NLONRGG_LOC(:)=NLONRGG(:)
      ENDIF

      IF (.NOT.(CLDOMAIN.EQ.'s' .OR. LWDINTS)) THEN
!       IF IT IS NOT A SWAMP CASE OR A CASE FOR WHICH A TIME SERIES
!       IS PRESCRIBED.
!       
!       WIND FIELDS ARE SUPPLIED FROM INPUT FILE sfcwindin
!       AS IS THE CASE IN THE STAND ALONE VERSION OR AT INITIAL INPUT
!       TIME IN COUPLED MODEL.

!       1.0 INITIALISE VALUES. 
!           ------------------ 

        ICOORD = 1

        FILNM='sfcwindin'
     
        IF(IRANK.EQ.IREAD) THEN
          LLEXIST=.FALSE.
          LNAME = LEN_TRIM(FILNM)
          INQUIRE(FILE=FILNM(1:LNAME),EXIST=LLEXIST)
          IF(.NOT. LLEXIST) THEN
            WRITE (IU06,*) '*************************************'
            WRITE (IU06,*) '*                                   *'
            WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
            WRITE (IU06,*) '*  IN READWIND:                     *'
            WRITE (IU06,*) '*  COULD NOT FIND FILE ',FILNM
            WRITE (IU06,*) '*                                   *'
            WRITE (IU06,*) '*************************************'
            WRITE (*,*) '*************************************'
            WRITE (*,*) '*                                   *'
            WRITE (*,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
            WRITE (*,*) '*  IN READWIND:                     *'
            WRITE (*,*) '*  COULD NOT FIND FILE ',FILNM
            WRITE (*,*) '*                                   *'
            WRITE (*,*) '*************************************'
            CALL ABORT1
          ENDIF

          IF (LLNOTOPENED) THEN
            CALL IGRIB_OPEN_FILE(IUNITW,FILNM(1:LNAME),'r') 
            LLNOTOPENED = .FALSE.
          ENDIF
        ENDIF  


! --------------------------------------------------------------------  

!*      2.0 READ MARS U and V WIND COMPONENTS  IN GRIB CODE FORMAT.
!           AND SEA ICE FACTION (OR SST) IF LICERUN=TRUE
!           WIND SPEED FROM PREVIOUS WAVE MODEL RUNS (IF LLWSWAVE)
!           WIND DIRECTION FROM PREVIOUS WAVE MODEL RUNS (IF LLWDWAVE)
!           SURFACE AIR DENSITY (IF LADEN AND NOT COUPLED TO IFS)
!           GUSTINESS (IF LGUST AND NOT COUPLED TO IFS)
!           -------------------------------------------------------

        IWPER = 1

        IF(LICERUN) THEN
          NFLD=3
        ELSE
          NFLD=2
          IPARAMCI=31
        ENDIF

        IF(LICETH) NFLD=NFLD+1
        IF(LLWSWAVE) NFLD=NFLD+1
        IF(LLWDWAVE) NFLD=NFLD+1

!       AIR DENSITY AND GUSTINESS ARE ONLY PROVIDED FROM A FILE IF STAND ALONE RUN
!!      for coupled run we re-initialise witht he value provided by IFS, even though
!!      strickly speaking, it should be provided as part of the restart but it is better than
!!      setting it to a constant !!
        IF(.NOT.LWCOU .AND. LADEN) NFLD=NFLD+1

!!      for coupled run we re-initialise witht he value provided by IFS, even though
!!      strickly speaking, it should be provided as part of the restart but it is better than
!!      setting it to a constant !!
        IF(.NOT.LWCOU .AND. LGUST) NFLD=NFLD+1

        ALLOCATE(WORK(NXFF,NYFF))
        ALLOCATE(LLNOTREAD(NFLD))
        LLNOTREAD=.TRUE.

        LLABORT=.FALSE.


!       LOOP OVER INPUT

        WND: DO IVAR=1,NFLD
          IF (ITEST.GT.1) THEN
            WRITE(IU06,*) ' SUB. READWIND - LOOP OVER IVAR ',IVAR,NFLD
            CALL FLUSH(IU06)
          ENDIF

2002      CONTINUE

          IF(IRANK.EQ.IREAD) THEN
1021        ISIZE=NBITW
            KBYTES=ISIZE*NPRECI
            IF(.NOT.ALLOCATED(KGRIB)) ALLOCATE(KGRIB(ISIZE))
            CALL IGRIB_READ_FROM_FILE(IUNITW,KGRIB,KBYTES,IRET)
            IF(IRET.EQ.JPGRIB_BUFFER_TOO_SMALL) THEN
!!!           *IGRIB_READ_FROM_FILE* does not read through the file if
!!!            the size is too small, so figure out the size and read again.
              CALL KGRIBSIZE(IU06, KBYTES, NBITW, 'READWIND')
              DEALLOCATE(KGRIB)
              GOTO 1021
            ELSEIF(IRET.EQ.JPGRIB_END_OF_FILE) THEN
              WRITE(IU06,*) '**********************************'
              WRITE(IU06,*) '* READWIND: END OF FILE ENCOUNTED'
              WRITE(IU06,*) '**********************************'
              CALL ABORT1
            ELSEIF(IRET.NE.JPGRIB_SUCCESS) THEN
              WRITE(IU06,*) '**********************************'
              WRITE(IU06,*) '* READWIND: FILE HANDLING ERROR'
              WRITE(IU06,*) '**********************************'
              CALL ABORT1
            ENDIF
          ENDIF

          IF (ITEST.GT.1) THEN
            WRITE(IU06,*) ' SUB. READWIND - READ FROM ',FILNM
          ENDIF

          CALL GSTATS(622,0)
          IF(LMESSPASS) CALL MPL_BARRIER(CDSTRING='READWIND: KGRIB ')

!         SEND GRIB DATA TO THE OTHER PE'S
!         --------------------------------
          IF(LMESSPASS .AND. NPROC.GT.1) THEN
            IF(IRANK.EQ.IREAD) THEN
              IDUM(1)=ISIZE
            ENDIF
            CALL MPL_BROADCAST(IDUM(1:1),KROOT=IREAD,KTAG=IVAR,
     &                         CDSTRING='READWIND IDUM:')
            IF(IRANK.NE.IREAD) THEN
              ISIZE=IDUM(1)
              ALLOCATE(KGRIB(ISIZE))
            ENDIF

            CALL MPL_BROADCAST(KGRIB(1:ISIZE),KROOT=IREAD,KTAG=IVAR,
     &                         CDSTRING='READWIND KGRIB:')

          ENDIF
          CALL GSTATS(622,1)

! ----------------------------------------------------------------------

!*        3.0 UNPACK MARS FIELDS.
!             -------------------
          KGRIB_HANDLE=-99
          CALL IGRIB_NEW_FROM_MESSAGE(KGRIB_HANDLE,KGRIB)
          ZDUM=0.
          CALL GRIB2WGRID (IU06, ITEST, NPROMA_WAM, 
     &                     KGRIB_HANDLE, KGRIB, ISIZE,
     &                     LLUNSTR,
     &                     NXFF, NYFF, NLONRGG_LOC,
     &                     IRGG, XDELLA, ZDELLO,
     &                     FIELDG%XLON, FIELDG%YLAT,
     &                     ZMISS, ZDUM, ZDUM,
     &                     CDTWIR, IFORP, IPARAM, KZLEV,IDM,IDM,WORK)

          CALL IGRIB_RELEASE(KGRIB_HANDLE)

          IF (IPARAM.EQ.165 .OR. IPARAM.EQ.33 .OR. IPARAM.EQ.131 .OR.
     &        IPARAM.EQ.180 ) THEN

            IF(LLNOTREAD(IVAR)) THEN
              DO J=1,NYFF
                JSN=NYFF-J+1
                DO I=1,NLONRGG_LOC(JSN)
                  FIELDG(I,J)%UWND=WORK(I,J)
                ENDDO
              ENDDO
              LLNOTREAD(IVAR)=.FALSE.
            ELSE
              LLABORT=.TRUE.
            ENDIF

            IF (IPARAM.EQ.180) THEN
              ICODE = 2
            ELSE
              ICODE = 3
            ENDIF
          ELSEIF (IPARAM.EQ.166 .OR. IPARAM.EQ.34 .OR. IPARAM.EQ.132
     &            .OR. IPARAM.EQ.181 ) THEN

            IF(LLNOTREAD(IVAR)) THEN 
              DO J=1,NYFF
                JSN=NYFF-J+1
                DO I=1,NLONRGG_LOC(JSN)
                  FIELDG(I,J)%VWND=WORK(I,J)
                ENDDO
              ENDDO
              LLNOTREAD(IVAR)=.FALSE.
            ELSE
              LLABORT=.TRUE.
            ENDIF

            IF (IPARAM.EQ.181) THEN
              ICODE = 2
            ELSE
              ICODE = 3
            ENDIF

          ELSEIF (IPARAM.EQ.31 .OR. IPARAM.EQ.139) THEN
             IPARAMCI=IPARAM

            IF(LLNOTREAD(IVAR)) THEN 
              DO J=1,NYFF
                JSN=NYFF-J+1
                DO I=1,NLONRGG_LOC(JSN)
                  FIELDG(I,J)%CIFR=WORK(I,J)
                ENDDO
              ENDDO
              LLNOTREAD(IVAR)=.FALSE.
            ELSEIF (.NOT.LICERUN .AND.
     &              (IPARAM.EQ.31.OR.IPARAM.EQ.139) ) THEN
!             SKIP SEA ICE MASK INFORMATION AS IT IS NOT NEEDED
              GOTO 2002
            ELSE
              LLABORT=.TRUE.
            ENDIF

          ELSEIF (IPARAM.EQ.92) THEN

            IF(LLNOTREAD(IVAR)) THEN 
              DO J=1,NYFF
                JSN=NYFF-J+1
                DO I=1,NLONRGG_LOC(JSN)
                  FIELDG(I,J)%CITH=WORK(I,J)
                ENDDO
              ENDDO
              LLNOTREAD(IVAR)=.FALSE.
            ELSEIF (.NOT.LICETH .AND.
     &              (IPARAM.EQ.92) ) THEN
!             SKIP SEA ICE THICKNESS INFORMATION AS IT IS NOT NEEDED
              GOTO 2002
            ELSE
              LLABORT=.TRUE.
            ENDIF

          ELSEIF (IPARAM.EQ.245) THEN

            IF (LLWSWAVE) THEN
              IF(LLNOTREAD(IVAR)) THEN 
                DO J=1,NYFF
                  JSN=NYFF-J+1
                  DO I=1,NLONRGG_LOC(JSN)
                    IF(WORK(I,J).EQ.ZMISS) THEN
                      FIELDG(I,J)%WSWAVE=0. 
                    ELSE
                      FIELDG(I,J)%WSWAVE=WORK(I,J)
                    ENDIF
                  ENDDO
                ENDDO
                LLNOTREAD(IVAR)=.FALSE.
              ELSE
                LLABORT=.TRUE.
              ENDIF
            ENDIF

          ELSEIF (IPARAM.EQ.249) THEN

            IF (LLWDWAVE) THEN
              IF(LLNOTREAD(IVAR)) THEN 
                DO J=1,NYFF
                  JSN=NYFF-J+1
                  DO I=1,NLONRGG_LOC(JSN)
                    IF(WORK(I,J).EQ.ZMISS) THEN
                      FIELDG(I,J)%WDWAVE=0. 
                    ELSE
!                     re-convert to WAM convention
                      FIELDG(I,J)%WDWAVE=RAD*(WORK(I,J)-180.)
                    ENDIF
                  ENDDO
                ENDDO
                LLNOTREAD(IVAR)=.FALSE.
              ELSE
                LLABORT=.TRUE.
              ENDIF
            ENDIF

!!!       215 is a pseudo parameter number for surface air density
!!!       adapt accordingly if such parameter were to be produced
          ELSEIF (IPARAM.EQ.215) THEN

            IF(LLNOTREAD(IVAR)) THEN 
              DO J=1,NYFF
                JSN=NYFF-J+1
                DO I=1,NLONRGG_LOC(JSN)
                  IF(WORK(I,J).EQ.ZMISS) THEN
                    FIELDG(I,J)%AIRD=ROAIR
                  ELSE
                    FIELDG(I,J)%AIRD=WORK(I,J)
                  ENDIF
                ENDDO
              ENDDO
              LLNOTREAD(IVAR)=.FALSE.
            ELSE
              LLABORT=.TRUE.
            ENDIF

!!!       216 is a pseudo parameter number for gustiness (Zi/L) 
!!!       adapt accordingly if such parameter were to be produced
          ELSEIF (IPARAM.EQ.216) THEN

            IF(LLNOTREAD(IVAR)) THEN 
              DO J=1,NYFF
                JSN=NYFF-J+1
                DO I=1,NLONRGG_LOC(JSN)
                  IF(WORK(I,J).EQ.ZMISS) THEN
                    FIELDG(I,J)%ZIDL = WSTAR0
                  ELSE
                    FIELDG(I,J)%ZIDL=WORK(I,J)
                  ENDIF
                ENDDO
              ENDDO
              LLNOTREAD(IVAR)=.FALSE.
            ELSE
              LLABORT=.TRUE.
            ENDIF

          ELSE
            WRITE(IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*) ' +                                        +'
            WRITE(IU06,*) ' +    WARNING ERROR IN SUB. READWIND      +'
            WRITE(IU06,*) ' +    ==============================      +'
            WRITE(IU06,*) ' + SUSPICIOUS WIND OR SEA ICE FIELD PARAM +'
            WRITE(IU06,*) ' + PARAM IS = ', IPARAM
            WRITE(IU06,*) ' +                                        +'
            WRITE(IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++'
            CALL ABORT1
          ENDIF

          IF(LLABORT) THEN
            WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*) ' +                                         +'
            WRITE(IU06,*) ' +    WARNING ERROR IN SUB. READWIND       +'
            WRITE(IU06,*) ' +    ==============================       +'
            WRITE(IU06,*) ' + WIND OR SEA ICE FIELD PARAM READ TWICE  +'
            WRITE(IU06,*) ' + PARAM IS = ', IPARAM
            WRITE(IU06,*) ' + SUSPECT                                 +'
            WRITE(IU06,*) ' + INCOMPLETE LIST OF INPUT PARAMETERS     +'
            WRITE(IU06,*) ' + THE LIST IS :                           +'
            WRITE(IU06,*) ' + 165 or 33 or 131 or 180 U-WIND COMPONENT+'
            WRITE(IU06,*) ' + 166 or 34 or 132 or 181 V-WIND COMPONENT+'
            IF(LICERUN) THEN
            WRITE(IU06,*) ' + 31 or 139 SEA ICE FRACTION OR SST       +'
            ENDIF
            IF(LICETH) THEN
            WRITE(IU06,*) ' + 92 SEA ICE THICKNESS                    +'
            ENDIF
            IF(LLWSWAVE) THEN
            WRITE(IU06,*) ' + 245  WIND SPEED FROM WAVE MODEL         +'
            ENDIF
            IF(LLWDWAVE) THEN
            WRITE(IU06,*) ' + 249  WIND DIRECTION FROM WAVE MODEL     +'
            ENDIF
            IF(LADEN) THEN
            WRITE(IU06,*) ' + 215  AIR DENSITY AT THE SURFACE         +'
            ENDIF
            IF(LGUST) THEN
            WRITE(IU06,*) ' + 216  GUSTINESS (Zi/L)                   +'
            ENDIF
            WRITE(IU06,*) ' +                                         +'
            WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++'
            CALL ABORT1
          ENDIF

!*      3.2.1 DETERMINE GRID PARAMETERS.                
!             -------------------------- 

          KCOL = NXFF
          KROW = NYFF
          RLATN = AMONOP
          RLATS = AMOSOP
          RLONL = AMOWEP
          RLONR = AMOEAP
          IF(KCOL.NE.1) THEN
            DLAM  = (RLONR-RLONL)/REAL(KCOL-1)
          ELSE
            DLAM  = 0.
          ENDIF
          IF(KROW.NE.1) THEN
            DPHI  = (RLATN-RLATS)/REAL(KROW-1)
          ELSE
            DPHI  = 0.
          ENDIF

          DEALLOCATE(KGRIB)
        ENDDO WND

        DEALLOCATE(WORK)
        DEALLOCATE(LLNOTREAD)

      ELSE
!       SWAMP CASE : CREATE YOUR OWN WIND FIELD
!       =======================================
!       for the swamp case create your own wind fields
!       or if the file windforcing_time_series is present then
!       the wind forcing will be extracted from the wind speed
!       and direction time series specified herewith.

        IF (LLNOTOPENED) THEN
           RLATS = AMOSOP
           RLATN = AMONOP 
           RLONL = AMOWEP 
           RLONR = AMOEAP 
           ICODE = 3
           CDTWIR = CDATEA
           CALL INCDATE(CDTWIR,-IDELWI)
           KCOL = NXFF
           KROW = NYFF
           IWPER = 0
           ICOORD = 1
           IF(KCOL.NE.1) THEN
             DLAM  = (RLONR-RLONL)/REAL(KCOL-1)
           ELSE
             DLAM  = 0.
           ENDIF
           IF(KROW.NE.1) THEN
             DPHI  = (RLATN-RLATS)/REAL(KROW-1)
           ELSE
             DPHI  = 0.
           ENDIF

           IF(LWDINTS) THEN
             IUNITW=I_GET_UNIT(IU06,CWDFILE, 'r', 'f', 0)
             OPEN(IUNITW,FILE=CWDFILE,FORM='FORMATTED')
             CWDDATE=CDATEA
           ELSE
             IUNITW=-1
           ENDIF
  
           LLNOTOPENED = .FALSE.
  
        ENDIF

        CALL INCDATE(CDTWIR,IDELWI)

        IF(LWDINTS) THEN

           READ(IUNITW,*,END=110,ERR=120) IWTIME,WSPEED,WTHETA

           CALL INCDATE(CWDDATE,IWTIME)
           WRITE(IU06,'(a28,a14,1x,f6.2,1x,f6.1)')
     &     '  WIND INPUT TIME SERIES AT ', CWDDATE,WSPEED,WTHETA

           WTHETA=RAD*WTHETA
           UWIND=-WSPEED*SIN(WTHETA)
           VWIND=-WSPEED*COS(WTHETA)
        ELSE
!         THE TRADITIONAL ONE GRID POINT FORCING

          CDTTURN = CDATEA
          IDTTURN=DTNEWWIND*3600
          CALL INCDATE(CDTTURN,IDTTURN)
          IF(SWAMPWIND2.LE.0. .OR. CDTWIR.LE.CDTTURN) THEN
            WRITE(IU06,*) ' READWIND - SWAMP CASE WITH WIND SPEED = ',
     &                    SWAMPWIND 
            UWIND=0.0
            VWIND=SWAMPWIND
          ELSE
            WRITE(IU06,*) ' READWIND - SWAMP CASE WITH WIND SPEED = ',
     &                    SWAMPWIND2 
            IF(LTURN90) THEN
              UWIND=SWAMPWIND2
              VWIND=0.0
            ELSE
              UWIND=0.0
              VWIND=SWAMPWIND2
            ENDIF
          ENDIF
        ENDIF

        IPARAMCI=31

#ifdef MODEL_COUPLING_ATM_WAV
        DO J=1,NYFF
          DO I=1,NXFF
            FIELDG(I,J)%UWND=FIELDG_coupl(I,J)%UWND
            FIELDG(I,J)%VWND=FIELDG_coupl(I,J)%VWND
            FIELDG(I,J)%CITH=FIELDG_coupl(I,J)%CITH
            FIELDG(I,J)%AIRD=FIELDG_coupl(I,J)%AIRD
            FIELDG(I,J)%ZIDL=FIELDG_coupl(I,J)%ZIDL
            FIELDG(I,J)%UCUR=FIELDG_coupl(I,J)%UCUR
            FIELDG(I,J)%VCUR=FIELDG_coupl(I,J)%VCUR
          ENDDO
        ENDDO
#else
        DO J=1,NYFF
          DO I=1,NXFF
            FIELDG(I,J)%UWND=UWIND
            FIELDG(I,J)%VWND=VWIND
          ENDDO
        ENDDO
#endif

!!!! impose the northern part of the swamp domain to be covered with
!!!! a sea ice cover=SWAMPCIFR
        DO J=1,NYFF/2
          DO I=1,NXFF
            FIELDG(I,J)%CIFR=SWAMPCIFR
          ENDDO
        ENDDO
        DO J=NYFF/2+1,NYFF
          DO I=1,NXFF
            FIELDG(I,J)%CIFR=0.
          ENDDO
        ENDDO

      ENDIF

#if !defined MODEL_COUPLING_ATM_WAV && defined NETCDF_OUTPUT_WAM
        DO J=1,NYFF
          DO I=1,NXFF
            FIELDG_coupl(I,J)%UWND=FIELDG(I,J)%UWND
            FIELDG_coupl(I,J)%VWND=FIELDG(I,J)%VWND
            FIELDG_coupl(I,J)%CITH=FIELDG(I,J)%CITH
            FIELDG_coupl(I,J)%AIRD=FIELDG(I,J)%AIRD
            FIELDG_coupl(I,J)%ZIDL=FIELDG(I,J)%ZIDL
            FIELDG_coupl(I,J)%UCUR=FIELDG(I,J)%UCUR
            FIELDG_coupl(I,J)%VCUR=FIELDG(I,J)%VCUR
          ENDDO
        ENDDO
#endif

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('READWIND',1,ZHOOK_HANDLE)
#endif
      RETURN

110   CONTINUE
      WRITE(IU06,*) '****************************************'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*  ERROR IN READWIND:                  *'
      WRITE(IU06,*) '*  END OF FILE REACHED FOR ',CWDFILE
      WRITE(IU06,*) '*  THE WIND TIME SERIES SHOULD AT LEAST*'
      WRITE(IU06,*) '*  BE AS LONG AS THE INTENDED RUN !!!! *'
      WRITE(IU06,*) '*  THE RUN HAS ABORTED                 *'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '****************************************'
      CALL ABORT1

120   CONTINUE
      WRITE(IU06,*) '****************************************'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*  ERROR IN READWIND:                  *'
      WRITE(IU06,*) '*  ERROR READING FILE ',CWDFILE
      WRITE(IU06,*) '*  THE RUN HAS ABORTED                 *'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '****************************************'
      CALL ABORT1

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('READWIND',1,ZHOOK_HANDLE)
#endif
      RETURN
      END SUBROUTINE READWIND
