      SUBROUTINE READWIND (CDTWIR, FILNM, LLNOTOPENED,IREAD)

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
!     *PBGRIB*
!     *PBOPEN*
!     *PBCLOSE*
!     *KGRIBSIZE*
!     *FLUSH*
!     *GRB2WGRD*       TRANSFORM GRIB FIELDS TO WAM GRID FIELDS

! --------------------------------------------------------------------- 

      USE YOWCOUT  , ONLY : LGRIBAPI_I
      USE YOWCOUP  , ONLY : LWCOU
      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK3    ,PSEK3
      USE YOWGRID  , ONLY : IGL      ,NLONRGG
      USE YOWICE   , ONLY : LICERUN  ,IPARAMCI
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,IXLG     , 
     &            KXLT
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPRECI
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NIBLO    ,CLDOMAIN ,
     &            SWAMPWIND,SWAMPWIND2,DTNEWWIND,LTURN90 ,LWDINTS
      USE YOWSTAT  , ONLY : CDATEA   ,IDELWI   ,LADEN    ,LGUST    ,
     &            NWAM_BLKS
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWNDG  , ONLY : DLAM     ,DPHI     ,RLATS    ,RLATN    ,
     &            RLONL    ,RLONR    ,KCOL     ,KROW     ,ICODE    ,
     &            IWPER    ,ICOORD
      USE YOWWIND  , ONLY : NC       ,NR       ,FIELDG   ,IUNITW   ,
     &            NBITW    ,CWDFILE  ,LLWSWAVE ,LLWDWAVE
      USE YOWPCONS , ONLY : RAD      ,ZMISS    ,ROAIR    ,WSTAR0

      USE MPL_MODULE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE GRIB_API_INTERFACE

! --------------------------------------------------------------------- 

      IMPLICIT NONE

      INTEGER :: IREAD 
      INTEGER :: NFLD  
      INTEGER :: I, J, IVAR, JSN
      INTEGER :: ISIZE, KRET, KPLENG, KLEN, KLENG 
      INTEGER :: IFORP, IPARAM, KZLEV, IDM 
      INTEGER :: IWTIME, IDTTURN
      INTEGER :: LNAME 
      INTEGER :: I_GET_UNIT
      INTEGER :: KBYTES, IRET, KDUM
      INTEGER :: KGRIB_HANDLE
      INTEGER :: IDUM(2)
      INTEGER, ALLOCATABLE :: KGRIB(:)

      REAL :: ZDUM
      REAL :: UWIND, VWIND, WSPEED, WTHETA
      REAL :: PMISS 
      REAL :: ZHOOK_HANDLE
      REAL, ALLOCATABLE :: WORK(:,:)  

      CHARACTER*24 :: FILNM
      CHARACTER*14 :: CDTWIR, CDTTURN
      CHARACTER*14, SAVE :: CWDDATE

      LOGICAL :: LLNOTOPENED
      LOGICAL :: LLABORT
      LOGICAL :: LLEXIST
      LOGICAL, ALLOCATABLE :: LLNOTREAD(:)

! --------------------------------------------------------------------  

      IF (LHOOK) CALL DR_HOOK('READWIND',0,ZHOOK_HANDLE)

      IF (ITEST.GT.1) THEN
        WRITE(IU06,*) ' SUB. READWIND - STARTING ' 
        CALL FLUSH(IU06)
      ENDIF
      NC=NGX
      NR=NGY

      IF(.NOT.(CLDOMAIN.EQ.'s' .OR. LWDINTS)) THEN
!       IF IT IS NOT A SWAMP CASE OR A CASE FOR WHICH A TIME SERIES
!       IS PRESCRIBED.

!       WIND FIELDS ARE SUPPLIED FROM INPUT FILE sfcwindin
!       AS IS THE CASE IN THE STAND ALONE VERSION OR AT INITIAL INPUT
!       TIME IN COUPLED MODEL.

!       1.0 INITIALISE VALUES. 
!           ------------------ 

        ICOORD = 1

        FILNM='sfcwindin'
     
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
          CALL ABORT1
        ENDIF

        IF(IRANK.EQ.IREAD) THEN
          IF (LLNOTOPENED) THEN
            IF(LGRIBAPI_I) THEN
              CALL IGRIB_OPEN_FILE(IUNITW,FILNM(1:LNAME),'r')
            ELSE
              CALL PBOPEN(IUNITW,FILNM,'r',KRET)
              IF(KRET.LT.0) THEN
                WRITE (IU06,*) '**************************************'
                WRITE (IU06,*) '*                                    *'
                WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO PBOPEN    *'
                WRITE (IU06,*) '*  IN READWIND                       *'
                IF(KRET.EQ.-1)
     &           WRITE (IU06,*) ' COULD NOT OPEN FILE ',FILNM
                IF(KRET.EQ.-2) WRITE (IU06,*)'INVALID FILENAME ',FILNM
                IF(KRET.EQ.-3) WRITE (IU06,*)'INVALID OPEN MODE '
                WRITE (IU06,*) '*                                    *'
                WRITE (IU06,*) '**************************************'
                CALL ABORT1
              ENDIF
            ENDIF
 
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
          DO J=1,NGY
            DO I=1,NGX
              FIELDG(I,J)%CIFR = 0.
            ENDDO
          ENDDO
          NFLD=2
          IPARAMCI=31
        ENDIF
        IF(LLWSWAVE) NFLD=NFLD+1
        IF(LLWDWAVE) NFLD=NFLD+1

!       AIR DENSITY AND GUSTINESS ARE ONLY PROVIDED FROM A FILE IF STAND ALONE RUN
        IF(.NOT.LWCOU .AND. LADEN) THEN
          NFLD=NFLD+1 
        ELSE
          DO J=1,NGY
            DO I=1,NGX
              FIELDG(I,J)%AIRD=ROAIR ! sets air density field to the default value.
            ENDDO
          ENDDO
        ENDIF

        IF(.NOT.LWCOU .AND. LGUST) THEN
          NFLD=NFLD+1 
        ELSE
          DO J=1,NGY
            DO I=1,NGX
              FIELDG(I,J)%ZIDL = WSTAR0
            ENDDO
          ENDDO
        ENDIF

        ALLOCATE(WORK(NGX,NGY))
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

          IF(LGRIBAPI_I) THEN
!           USE GRIB API
            IF(IRANK.EQ.IREAD) THEN
1021          ISIZE=NBITW
              KBYTES=ISIZE*NPRECI
              IF(.NOT.ALLOCATED(KGRIB)) ALLOCATE(KGRIB(ISIZE))
              CALL IGRIB_READ_FROM_FILE(IUNITW,KGRIB,KBYTES,IRET)
              IF(IRET.EQ.JPGRIB_BUFR_TOO_SHORT) THEN
!!!             *IGRIB_READ_FROM_FILE* does not read through the file if
!!!              the size is too small, so figure out the size and read again.
                KDUM=-1
                CALL KGRIBSIZE(IU06, KDUM, KBYTES, NBITW, 'READWIND')
                DEALLOCATE(KGRIB)
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

          ELSE
!           USE GRIBEX
            ISIZE=NBITW
            PMISS=PSEK3(2)

            IF(IRANK.EQ.IREAD) THEN
2001          KPLENG=ISIZE*KIND(ISIZE)
              IF(.NOT.ALLOCATED(KGRIB)) ALLOCATE(KGRIB(ISIZE))

              CALL GSTATS(1740,0)
              CALL PBGRIB(IUNITW,KGRIB,KPLENG,KLEN,KRET)
              CALL GSTATS(1740,1)

              IF (KRET.EQ.-1) THEN
                WRITE(IU06,*) ' ***************************************'
                WRITE(IU06,*) ' *    FATAL ERROR IN SUB. READWIND     *'
                WRITE(IU06,*) ' *    ============================     *'
                WRITE(IU06,*) ' * PBGRIB EOF ON FILE ',FILNM 
                WRITE(IU06,*) ' *       NO MORE WINDS                 *'
                WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS  *'
                WRITE(IU06,*) ' ***************************************'
                CALL ABORT1
              ELSE IF (KRET.EQ.-2) THEN
                WRITE(IU06,*) ' ***************************************'
                WRITE(IU06,*) ' *    FATAL ERROR IN SUB. READWIND     *'
                WRITE(IU06,*) ' *    ============================     *'
                WRITE(IU06,*) ' * PBGRIB ERROR ON FILE ',FILNM
                WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS  *'
                WRITE(IU06,*) ' ***************************************'
                CALL ABORT1
              ELSE IF (KRET.EQ.-3) THEN
                DEALLOCATE(KGRIB)
                CALL KGRIBSIZE(IU06, IUNITW, KLEN, ISIZE, 'READWIND')
                GOTO 2001
              ENDIF

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
              IDUM(2)=KLEN
            ENDIF
            CALL MPL_BROADCAST(IDUM(1:2),KROOT=IREAD,KTAG=IVAR,
     &                         CDSTRING='READWIND IDUM:')
            IF(IRANK.NE.IREAD) THEN
              ISIZE=IDUM(1)
              KLEN=IDUM(2)
              ALLOCATE(KGRIB(ISIZE))
            ENDIF

            CALL MPL_BROADCAST(KGRIB(1:ISIZE),KROOT=IREAD,KTAG=IVAR,
     &                         CDSTRING='READWIND KGRIB:')

          ENDIF
          CALL GSTATS(622,1)

! ----------------------------------------------------------------------

!*        3.0 UNPACK MARS FIELDS.
!             -------------------
          IF(LGRIBAPI_I) THEN
            KGRIB_HANDLE=-99
            CALL IGRIB_NEW_FROM_MESSAGE(KGRIB_HANDLE,KGRIB)
            ZDUM=0.
            CALL GRIB2WGRID (IU06, ITEST, NWAM_BLKS, 
     &                       KGRIB_HANDLE, KGRIB, ISIZE,
     &                       NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                       AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                       ZMISS, ZDUM, ZDUM,
     &                       CDTWIR, IFORP, IPARAM, KZLEV,IDM,IDM,WORK)

            CALL IGRIB_RELEASE(KGRIB_HANDLE)

          ELSE
            KLENG = KLEN
            ZDUM=0.
            CALL GRB2WGRD (IU06, ITEST, NKSEK0, NKSEK1, NKSEK2, NKSEK3,
     &                     NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &                     KGRIB, KLENG ,
     &                     NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                     AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                     ZDUM, ZDUM,
     &                     CDTWIR, IFORP, IPARAM, KZLEV, IDM, IDM, WORK)
          ENDIF

          IF (IPARAM.EQ.165 .OR. IPARAM.EQ.33 .OR. IPARAM.EQ.131 .OR.
     &        IPARAM.EQ.180 ) THEN

            IF(LLNOTREAD(IVAR)) THEN 
              DO J=1,NGY
                JSN=NGY-J+1
                DO I=1,NLONRGG(JSN)
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
              DO J=1,NGY
                JSN=NGY-J+1
                DO I=1,NLONRGG(JSN)
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
              DO J=1,NGY
                JSN=NGY-J+1
                DO I=1,NLONRGG(JSN)
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

          ELSEIF (IPARAM.EQ.245) THEN

            IF (LLWSWAVE) THEN
              IF(LLNOTREAD(IVAR)) THEN 
                DO J=1,NGY
                  JSN=NGY-J+1
                  DO I=1,NLONRGG(JSN)
                    IF(WORK(I,J).EQ.PMISS) THEN
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
                DO J=1,NGY
                  JSN=NGY-J+1
                  DO I=1,NLONRGG(JSN)
                    IF(WORK(I,J).EQ.PMISS) THEN
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
              DO J=1,NGY
                JSN=NGY-J+1
                DO I=1,NLONRGG(JSN)
                  IF(WORK(I,J).EQ.PMISS) THEN
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
              DO J=1,NGY
                JSN=NGY-J+1
                DO I=1,NLONRGG(JSN)
                  IF(WORK(I,J).EQ.PMISS) THEN
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

          KCOL = NGX
          KROW = NGY
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
           KCOL = NC
           KROW = NR
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

           READ(IUNITW,100,END=110,ERR=120) IWTIME,WSPEED,WTHETA
100        FORMAT(i8,1x,f6.2,1x,f6.1)

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

        DO J=1,NR
          DO I=1,NC
            FIELDG(I,J)%UWND=UWIND
            FIELDG(I,J)%VWND=VWIND
          ENDDO
        ENDDO

        DO J=1,NR
          DO I=1,NC
            FIELDG(I,J)%AIRD=ROAIR
            FIELDG(I,J)%ZIDL=WSTAR0
            FIELDG(I,J)%CIFR=0.
            FIELDG(I,J)%UCUR=0.
            FIELDG(I,J)%VCUR=0.
          ENDDO
       ENDDO

      ENDIF

      IF (LHOOK) CALL DR_HOOK('READWIND',1,ZHOOK_HANDLE)
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

      IF (LHOOK) CALL DR_HOOK('READWIND',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE READWIND
