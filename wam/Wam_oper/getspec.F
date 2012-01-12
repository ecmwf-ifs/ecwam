      SUBROUTINE GETSPEC(FL,SL,NSTART,NEND,IREAD,U10OLD,THWOLD)
! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      SEPTEMBER 1997 
!     J. BIDLOT    ECMWF      MARCH 2010: modified to use gribapi 

!*    PURPOSE.
!     --------
!     READ THE SPECTRA FROM DISK.

!**   INTERFACE.
!     ----------
!     *CALL* *GETSPEC(FL,SL,NSTART,NEND,IREAD,U10OLD,THWOLD)
!     *FL*        ARRAY CONTAINING THE SPECTRA CONTRIBUTION ON EACH PE
!     *SL*        SPARE ARRAY TO MANIPUATE SPECTRA
!     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *IREAD*     PROCESSOR WHICH WILL ACCESS THE FILE ON DISK 
!     *U10OLD*    WIND SPEED. (used with fetch law to fill empty 
!                 sea points)
!     *THWOLD*    WIND DIRECTION (RADIANS).


!     METHOD.
!     -------

!     IN CASE THE INPUT SPECTRA ARE IN GRIB THEN EITHER THEY ARE READ
!     FROM A FILE ON DISK USING THE PBGRIB SOFTWARE.
!     GRIB SPECTRA WILL BE DECODED. IN CASE OF BINARY
!     DATA, USE READFL TO READ IN THE SPECTRA DEPENDING ON THE USE OF
!     THE PBIO SOFTWARE OR NOT. PBIO WILL LIMIT THE SIZE OF THE ARRAY
!     NECESSARY TO READ THE INPUT SPECTRA. THE READING IS ONLY DONE ON
!     PE 1, THEREFORE THE RELEVANT INFORMATION IS SENT TO THE OTHER
!     PE'S USING MPDISTRIBFL

!     EXTERNALS.
!     ----------
!     GETENV
!     GRIBEX
!     GRSTNAME
!     MPDISTRIBFL
!     MPDISTRIBSCFLD
!     MPL_BARRIER
!     PBOPEN
!     PBGRIB
!     PBCLOSE 
!     READFL

!     REFERENCE.
!     ----------
!     NONE

! ----------------------------------------------------------------------
      USE YOWCOUT  , ONLY : KDEL     ,MDEL     ,LGRIBAPI_I
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3  , 
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK3    ,PSEK3   ,
     &            NTENCODE ,PPEPS    ,PPREC
      USE YOWGRID  , ONLY : DELPHI   ,NLONRGG
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IRGG     ,AMOSOP   ,
     &            AMONOP   ,AMOWEP   ,AMOEAP   ,XDELLA   ,ZDELLO

      USE YOWMESPAS, ONLY : LMESSPASS,LPBIOIN  ,LFDBIOIN ,LGRIBIN
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,
     &            KTAG     ,NPRECR   ,NPRECI
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,
     &            NBLO     ,NIBLO    ,CLDOMAIN
      USE YOWPCONS , ONLY : G        ,DEG      ,R        ,ZMISS    ,
     &            EPSMIN
      USE YOWSTAT  , ONLY : CDATEF   ,CDTPRO   ,NWAM_BLKS,IREFRA   ,
     &            LNSESTART
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH    ,LRESTARTED
      USE MPL_MODULE
      USE FDBSUBS_MOD
      USE GRIB_API_INTERFACE

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: NBIT

      INTEGER :: ISEND, IREAD
      INTEGER :: IJ, K, M, IC, ICC, JSN, IDUM, IX, IY, ID, MR, KR, IP
      INTEGER :: KLOOP, MLOOP, KINF, KSUP, MINF, MSUP
      INTEGER :: IERR, KRCOUNT, KFROM, KRTAG
      INTEGER :: IRA, IG
      INTEGER :: IUNIT
      INTEGER :: IBREAD, NBREAD, NBREAD_AGAIN 
      INTEGER :: IPARAM, KZLEV,KK,MM
      INTEGER :: IYYYY, JCONS, IFORP, KDEXN
      INTEGER :: ISTEP, ISTEP_LOCAL
      INTEGER :: IUGRS, KRET, IPLENG, ISIZE, KLEN, ILENG, KWORD
      INTEGER :: KBYTES, IRET
      INTEGER :: LFILE, KFILE_HANDLE, KGRIB_HANDLE
      INTEGER :: JKGLO,KIJS,KIJL,NPROMA
      INTEGER :: IPROC, ITAG, IREQ, IST, IEND, KSEND
      INTEGER :: ISENDREQ(NPROC)
      INTEGER, DIMENSION(NPROC) :: NSTART, NEND
      INTEGER, ALLOCATABLE :: INGRIB(:)
      INTEGER, ALLOCATABLE :: IPOINTER(:)
      INTEGER, ALLOCATABLE :: KSEC0(:),KSEC1(:),KSEC2(:),KSEC3(:)
      INTEGER, ALLOCATABLE :: KSEC4(:)

      REAL :: ZTHETA,ZFRE
      REAL :: RMONOP
      REAL, DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL, SL
      REAL, DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD
      REAL, ALLOCATABLE, DIMENSION(:) :: PSEC2,PSEC3,PSEC4,WORK
      REAL, ALLOCATABLE, DIMENSION(:) :: ZRECVBUF
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: RFL
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FIELD

      CHARACTER(LEN=14) :: CDATE 
      CHARACTER(LEN=120) :: FILENAME

      LOGICAL :: LLEPSMIN
      LOGICAL :: LFRSDECODE, LPBOPEN, LPBCLOSE, LLEXIST
      LOGICAL :: LLRESIZING=.FALSE.

! ----------------------------------------------------------------------

      LFRSDECODE=.TRUE.

      NBIT=NIBLO

      IF (ITEST > 3 ) THEN
        WRITE(IU06,*) ' SUB. GETSPEC:'
        WRITE(IU06,*) ' ABOUT TO READ SPECTRUM IN........'
        WRITE(IU06,*) ' CDTPRO    =', CDTPRO
        WRITE(IU06,*) ' LGRIBIN   =', LGRIBIN
        WRITE(IU06,*) ' LPBIOIN   =', LPBIOIN
        WRITE(IU06,*) ' LFDBIOIN  =', LFDBIOIN
        WRITE(IU06,*) ' LRESTARTED=', LRESTARTED
        WRITE(IU06,*) ' IRANK=', IRANK
        WRITE(IU06,*) ' IREAD=', IREAD
        CALL FLUSH(IU06)
      ENDIF
C
      LPBOPEN = .TRUE.
      LPBCLOSE = .TRUE.
      ISEND=IREAD

      IF(LNSESTART .AND. .NOT.LRESTARTED) THEN
!     BY-PASSED INPUT BY STARTING WITH SPECTRA AT NOISE LEVEL
!     =======================================================

        NPROMA=(NEND(IRANK)-NSTART(IRANK))/NWAM_BLKS+1
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ,K,M)
        DO JKGLO=NSTART(IRANK),NEND(IRANK),NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,NEND(IRANK))
          DO IJ=KIJS,KIJL
            DO K=1,NANG
              DO M=1,NFRE
                FL(IJ,K,M) = EPSMIN 
              ENDDO
           ENDDO
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO

         DO K=1,NANG
           DO M=1,NFRE
             FL(NINF-1,K,M) = 0. 
           ENDDO
         ENDDO

      ELSEIF(LGRIBIN.AND..NOT.LRESTARTED) THEN
!     INPUT SPECTRA ARE IN GRIB FORMAT
!     ================================

        IF (LFDBIOIN) THEN
          WRITE(IU06,*)'**************************************'
          WRITE(IU06,*)'*                                    *'
          WRITE(IU06,*)'*GETSPEC : THIS OPTION IS OBSOLETE !!*'
          WRITE(IU06,*)'*PROGRAM WILL ABORT                  *'
          WRITE(IU06,*)'*                                    *'
          WRITE(IU06,*)'**************************************'
          CALL ABORT1
        ENDIF

        FILENAME='specwavein'
        LFILE=LEN_TRIM(FILENAME)

        INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
        IF(.NOT.LLEXIST) THEN
          WRITE(IU06,*)'**************************************'
          WRITE(IU06,*)'*                                    *'
          WRITE(IU06,*)'*GETSPEC : GRIB SPECTRA NOT FOUND IN *'
          WRITE(IU06,*)  FILENAME
          WRITE(IU06,*)'*PROGRAM WILL ABORT                  *'
          WRITE(IU06,*)'*                                    *'
          WRITE(IU06,*)'**************************************'
          CALL ABORT1
        ENDIF

        NBREAD=0
        NBREAD_AGAIN=0
1121    CONTINUE
        IF(LGRIBAPI_I) THEN
!       USE GRIB API
!       ------------
!         CONNECT INPUT PE (IREAD) WITH INPUT FILE
          IF (IRANK.EQ.IREAD) THEN
            CALL IGRIB_OPEN_FILE(KFILE_HANDLE,FILENAME(1:LFILE),'r')
          ENDIF
          IF (.NOT.ALLOCATED(WORK)) ALLOCATE(WORK(NIBLO))
        ELSE
!       USE GRIBEX
!       ----------
          IF (.NOT.ALLOCATED(KSEC0)) ALLOCATE(KSEC0(NKSEK0))
          IF (.NOT.ALLOCATED(KSEC1)) ALLOCATE(KSEC1(NKSEK1))
          IF (.NOT.ALLOCATED(KSEC2)) ALLOCATE(KSEC2(NKSEK2))
          IF (.NOT.ALLOCATED(KSEC3)) ALLOCATE(KSEC3(NKSEK3))
          IF (.NOT.ALLOCATED(KSEC4)) ALLOCATE(KSEC4(NKSEK4))

          IF (.NOT.ALLOCATED(PSEC2)) ALLOCATE(PSEC2(NPSEK2))
          IF (.NOT.ALLOCATED(PSEC3)) ALLOCATE(PSEC3(NPSEK3))

          IF (.NOT.ALLOCATED(IPOINTER))
     &           ALLOCATE(IPOINTER(1:NEND(NPROC)))
  
          IF (.NOT.ALLOCATED(WORK)) ALLOCATE(WORK(NIBLO))

          IF (IRANK.EQ.IREAD) THEN
!           CONNECT INPUT PE (IREAD) WITH INPUT FILE
            CALL PBOPEN(IUGRS,FILENAME(1:LFILE),'r',KRET)
            IF(KRET.LT.0) THEN 
              WRITE (IU06,*) '***********************************'
              WRITE (IU06,*) '*                                 *'
              WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBOPEN  *'
              WRITE (IU06,*) '* IN GETSPEC                      *'
              IF(KRET.EQ.-1)
     &           WRITE (IU06,*) ' COULD NOT OPEN FILE ',FILENAME
              IF(KRET.EQ.-2)
     &         WRITE (IU06,*) ' INVALID FILENAME ',FILENAME
              IF(KRET.EQ.-3)
     &         WRITE (IU06,*)' INVALID OPEN MODE SPECIFIED'
              WRITE (IU06,*) '*                                 *'
              WRITE (IU06,*) '***********************************'
              CALL ABORT1
            ENDIF
          ENDIF

        ENDIF

!       GET GRIB DATA FROM (NFRE*NANG) FIELDS


!       READ NPROC-1 FIELDS AND SEND THEM SUCCESSIVELY TO ALL OTHER PE'S
!       FOR DECODING (IF IN MESSAGE PASSING MODE)

        ISIZE=NBIT
        ISTEP=MAX(NPROC-1,1)

        ALL_FILE : DO IC=1,NFRE*NANG,ISTEP 

          ISTEP_LOCAL=ISTEP
          IF(IC+ISTEP.GT.NFRE*NANG ) ISTEP_LOCAL=NFRE*NANG-IC+1

          ALL_DECODE_PE : DO IDUM=1,ISTEP_LOCAL
            IF(NPROC.EQ.1) THEN
              KSEND=1
            ELSEIF(IDUM.LT.IREAD) THEN
              KSEND=IDUM
            ELSE
              KSEND=IDUM+1
            ENDIF
            M=(((IC-1)+IDUM-1)/NANG)+1
            K=(IC-1)+IDUM-(M-1)*NANG

!           DATA ARE READ IN ON PE IREAD

            IF(LGRIBAPI_I) THEN
!           USE GRIB API
!             LOAD THE DATA
              IF (IRANK.EQ.IREAD) THEN
1021            ISIZE=NBIT
                KBYTES=ISIZE*NPRECI
                IF(.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))
                NBREAD=NBREAD+1
                CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE,INGRIB,
     &                                    KBYTES,IRET)
                IF(IRET.EQ.JPGRIB_BUFR_TOO_SHORT) THEN
                  IF(.NOT.LLRESIZING) NBREAD_AGAIN=NBREAD
                  CALL KGRIBSIZE(IU06, IUGRS, KBYTES, NBIT, 'GETSPEC')
                  DEALLOCATE(INGRIB)
                  LLRESIZING=.TRUE.
                  GOTO 1021
                ELSEIF(LLRESIZING .AND. IRET.NE.JPGRIB_END_OF_FILE) THEN
!                 LOOP UNTIL YOU HAVE EXPLORE THE SIZE FOR THE WHOLE FILE.
                  DEALLOCATE(INGRIB)
                  GOTO 1021
                ELSEIF(LLRESIZING .AND. IRET.EQ.JPGRIB_END_OF_FILE) THEN
!                 WE SHOULD HAVE THE MAXIMUM SIZE NECESSARY, START ALL OVER.
                  DEALLOCATE(INGRIB)
                  LLRESIZING=.FALSE.
                  CALL IGRIB_CLOSE_FILE(KFILE_HANDLE)
                  CALL IGRIB_OPEN_FILE(KFILE_HANDLE,
     &                                 FILENAME(1:LFILE),'r')
                  ISIZE=NBIT
                  IF(.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))
                  DO IBREAD=1,NBREAD_AGAIN
                    KBYTES=ISIZE*NPRECI
                    CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE,INGRIB,
     &                                        KBYTES,IRET)
                  ENDDO
                  NBREAD=IBREAD-1
                  NBREAD_AGAIN=0
                ELSEIF(IRET.EQ.JPGRIB_END_OF_FILE) THEN
                  WRITE(IU06,*) '**********************************'
                  WRITE(IU06,*) '* GETSPEC: END OF FILE ENCOUNTED'
                  WRITE(IU06,*) '**********************************'
                  CALL ABORT1
                ELSEIF(IRET.NE.JPGRIB_SUCCESS) THEN
                  WRITE(IU06,*) '**********************************'
                  WRITE(IU06,*) '* GETSPEC: FILE HANDLING ERROR'
                  WRITE(IU06,*) '**********************************'
                  CALL ABORT1
                ENDIF
              ENDIF

            ELSE
              IF (IRANK.EQ.IREAD) THEN
!               READ GRIB DATA FROM FILE (IF REQUESTED)

1000            IPLENG=ISIZE*NPRECI
                IF(.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))

                CALL GSTATS(1772,0)
                CALL PBGRIB(IUGRS,INGRIB,IPLENG,ILENG,KRET)
                CALL GSTATS(1772,1)
                IF (KRET.EQ.-1) THEN
                  WRITE (IU06,*) ' SUB. GETSPEC : '
                  WRITE (IU06,*) ' REACHED EOF ON UNIT ',IUGRS
                  WRITE (IU06,*) ' PROGRAM WILL ABORT '
                  CALL ABORT1
                ELSEIF (KRET .EQ. -2) THEN
                  WRITE (IU06,*) ' SUB. GETSPEC : '
                  WRITE (IU06,*) 
     &               ' FILE HANDLING ERROR ON UNIT ',IUGRS
                  WRITE (IU06,*) ' PROGRAM WILL ABORT '
                  CALL ABORT1
                ELSEIF (KRET .EQ. -3) THEN
                  DEALLOCATE(INGRIB)
                  CALL KGRIBSIZE(IU06, IUGRS, ILENG, ISIZE, 'GETSPEC')
                  GOTO 1000
                ENDIF

              ENDIF
            ENDIF


            IF(LMESSPASS) THEN
!             IN CASE OF MESSAGE PASSING THE DECODING WILL OCCUR ON KSEND
!             SEND DATA TO KEND: 
              CALL GSTATS(623,0)
              IF(IRANK.EQ.IREAD.AND.NPROC.NE.1) THEN
!               SEND GRIB DATA SIZE TO PE KSEND
                ITAG=(M-1)*NANG+K
                CALL MPL_SEND(ISIZE,KDEST=KSEND,KTAG=ITAG,KERROR=IERR,
     &            CDSTRING='GETSPEC 0:')
                IF(IERR.LT.0) CALL MPL_ABORT
     &                    ('MPL_SEND ERROR AT 1 in GETSPEC' )
              ENDIF

              IF (IRANK.EQ.KSEND.AND.NPROC.NE.1) THEN
!               RECEIVE GRIB DATA SIZE FROM IREAD
                ITAG=(M-1)*NANG+K
                CALL MPL_RECV(ISIZE,KSOURCE=IREAD,KTAG=ITAG,
     &            KOUNT=KRCOUNT,KRECVTAG=KRTAG,KERROR=IERR,
     &            CDSTRING='GETSPEC 0:')
                IF(IERR.LT.0) CALL MPL_ABORT
     &                           ('MPL_RECV ERROR AT 1 in GETSPEC' )
                IF(KRTAG.NE.ITAG) CALL MPL_ABORT
     &          ('MPL_RECV ERROR AT 1 in GETSPEC:  MISMATCHED TAGS' )

                IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)
                ALLOCATE(INGRIB(ISIZE))
              ENDIF

              IF(IRANK.EQ.IREAD.AND.NPROC.NE.1) THEN
!               SEND GRIB DATA TO PE KSEND
                ITAG=NFRE*NANG+(M-1)*NANG+K
                CALL MPL_SEND(INGRIB(1:ISIZE),KDEST=KSEND,KTAG=ITAG,
     &           KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(1),
     &           KERROR=IERR,CDSTRING='GETSPEC 1:')
                IF(IERR.LT.0) CALL MPL_ABORT
     &                      ('MPL_SEND ERROR AT 2 in GETSPEC' )
              ENDIF

              IF (IRANK.EQ.KSEND.AND.NPROC.NE.1) THEN
!               RECEIVED GRIB DATA FROM PE IREAD
                ITAG=NFRE*NANG+(M-1)*NANG+K
                CALL MPL_RECV(INGRIB(1:ISIZE),KSOURCE=IREAD,KTAG=ITAG,
     &            KOUNT=KRCOUNT,KRECVTAG=KRTAG,KERROR=IERR,
     &            CDSTRING='GETSPEC 1:')
                IF(IERR.LT.0) CALL MPL_ABORT
     &                       ('MPL_RECV ERROR AT 2 in GETSPEC ' )
                IF(KRCOUNT.NE.ISIZE) CALL MPL_ABORT
     &          ('MPL_RECV ERROR in 2 in GETSPEC:MISMATCHED MSG LENGTH')
                IF(KRTAG.NE.ITAG) CALL MPL_ABORT
     &          ('MPL_RECV ERROR in 2 in GETSPEC:MISMATCHED TAGS' )

              ENDIF

              IF(IRANK.EQ.IREAD.AND.NPROC.NE.1) THEN
                CALL MPL_WAIT(INGRIB(1:ISIZE),KREQUEST=ISENDREQ(1), 
     &           CDSTRING='GETSPEC: WAIT FOR SEND')
                IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)
              ENDIF
              CALL GSTATS(623,1)

            ENDIF


!           DECODE THE GRIB DATA ON PE KSEND 
            IF (IRANK.EQ.KSEND.OR.NPROC.EQ.1) THEN

              IF(LGRIBAPI_I) THEN
!             USE GRIB API

                KGRIB_HANDLE=-99
                CALL IGRIB_NEW_FROM_MESSAGE(KGRIB_HANDLE,INGRIB)

                IF(.NOT.ALLOCATED(FIELD)) ALLOCATE(FIELD(NGX,NGY))
                CALL GRIB2WGRID (IU06, ITEST, NWAM_BLKS, 
     &                         KGRIB_HANDLE, INGRIB, ISIZE,
     &                         NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                         AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                         ZMISS, PPREC, PPEPS,
     &                         CDATE, IFORP, IPARAM, KZLEV,KK,MM,FIELD)

                CALL IGRIB_RELEASE(KGRIB_HANDLE)

                IF (CDATE.NE.CDTPRO) THEN
                  WRITE(IU06,*)'**********************************'
                  WRITE(IU06,*)'*                                *'
                  WRITE(IU06,*)'* FATAL ERROR IN SUB GETSPEC     *'
                  WRITE(IU06,*)'* ===========================    *'
                  WRITE(IU06,*)'*                                *'
                  WRITE(IU06,*)'* REQUESTED DATE IS NOT EQUAL TO *'
                  WRITE(IU06,*)'* RETRIEVED DATE.                *'
                  WRITE(IU06,*)'* IN FILE: ',FILENAME
                  WRITE(IU06,*)'* CDATE = ',CDATE
                  WRITE(IU06,*)'* CDTPRO = ',CDTPRO
                  WRITE(IU06,*)'*                                *'
                  WRITE(IU06,*)'**********************************'
                  CALL ABORT1
                ENDIF
                IF(K.NE.KK) THEN
                  WRITE(IU06,*) '************************************'
                  WRITE(IU06,*) '* FATAL ERROR IN SUB. GETSPEC      *'
                  WRITE(IU06,*) '* REQUESTED AND DECODED DIRECTIONAL*'
                  WRITE(IU06,*) '* INDEX ARE DIFFERENT :            *'
                  WRITE(IU06,*) '* REQUESTED : ',K 
                  WRITE(IU06,*) '* DECODED   : ',KK
                  WRITE(IU06,*) '*                                  *'
                  WRITE(IU06,*) '************************************'
                  CALL ABORT1
                ENDIF
                IF(M.NE.MM) THEN
                  WRITE(IU06,*) '************************************'
                  WRITE(IU06,*) '* FATAL ERROR IN SUB. GETSPEC      *'
                  WRITE(IU06,*) '* REQUESTED AND DECODED FREQUENCY  *'
                  WRITE(IU06,*) '* INDEX ARE DIFFERENT :            *'
                  WRITE(IU06,*) '* REQUESTED : ',M 
                  WRITE(IU06,*) '* DECODED   : ',MM
                  WRITE(IU06,*) '*                                  *'
                  WRITE(IU06,*) '************************************'
                  CALL ABORT1
                ENDIF

                NPROMA=NEND(NPROC)/NWAM_BLKS+1
!$OMP           PARALLEL DO SCHEDULE(STATIC)
!$OMP+          PRIVATE(JKGLO,KIJS,KIJL,IJ,IX,IY)
                DO JKGLO=1,NEND(NPROC),NPROMA
                  KIJS=JKGLO
                  KIJL=MIN(KIJS+NPROMA-1,NEND(NPROC))
                  DO IJ=KIJS,KIJL
                    IX = IXLG(IJ,1)
                    IY = NGY- KXLT(IJ,1) +1
                    IF (FIELD(IX,IY) .NE. ZMISS) THEN
                      WORK(IJ) = FIELD(IX,IY)
                    ELSE
                      WORK(IJ) = EPSMIN
                    ENDIF
                  ENDDO
                ENDDO
!$OMP           END PARALLEL DO
                DEALLOCATE(FIELD)

              ELSE

                KSEC3=KSEK3
                PSEC3=PSEK3
                ALLOCATE(PSEC4(NTENCODE))

                CALL GSTATS(1703,0)
                CALL GRIBEX (KSEC0,KSEC1,KSEC2,PSEC2,KSEC3,PSEC3,KSEC4,
     &                     PSEC4, NTENCODE,INGRIB,ISIZE,KWORD,'D',KRET)
                CALL GSTATS(1703,1)

                CALL GSTATS(1886,0)
                  IF(KRET.GT.0) THEN 
                  WRITE(IU06,*) '************************************'
                  WRITE(IU06,*) '* FATAL ERROR IN SUB. GETSPEC      *'
                  WRITE(IU06,*) '* GRIBEX ERROR - KRET = ',KRET
                  WRITE(IU06,*) '* WHILE DECODING INGRB             *'
                  WRITE(IU06,*) '*                                  *'
                  WRITE(IU06,*) '*    PROGRAM ABORTS                *'
                  WRITE(IU06,*) '*                                  *'
                  WRITE(IU06,*) '************************************'
                  CALL ABORT1
                ENDIF

                IYYYY=(KSEC1(21)-1)*100+KSEC1(10)
                WRITE (CDATE(1:12),'(I4.4,4I2.2)') IYYYY,KSEC1(11),
     &                                 KSEC1(12), KSEC1(13),KSEC1(14)
                CDATE(13:14)='00'

                IF(KSEC1(39).EQ.9.OR.KSEC1(39).EQ.6.OR.
     &             KSEC1(39).EQ.10.OR.KSEC1(39).EQ.11)
     &          THEN
                  IF (KSEC1(15).EQ.0) THEN
                    JCONS = 60
                  ELSEIF (KSEC1(15).EQ.1) THEN
                    JCONS = 3600
                  ELSEIF (KSEC1(15).EQ.2) THEN
                    JCONS = 86400
                  ELSEIF (KSEC1(15).EQ.10) THEN
                    JCONS = 10800
                  ELSEIF (KSEC1(15).EQ.11) THEN
                    JCONS = 21600
                  ELSEIF (KSEC1(15).EQ.12) THEN
                    JCONS = 43200
                  ELSE
                    WRITE(IU06,*)'**********************************'
                    WRITE(IU06,*)'* FATAL ERROR IN SUB GETSPEC     *'
                    WRITE(IU06,*)'* SUSPICIOUS TIME UNIT           *' 
                    WRITE(IU06,*)'**********************************'
                    CALL ABORT1
                  ENDIF
!                 DETERMINE FORECAST PERIOD.
                  IFORP = JCONS*KSEC1(16)

                  CALL INCDATE (CDATE,IFORP)
                ENDIF
  
                IF (CDATE.NE.CDTPRO) THEN
                  WRITE(IU06,*)'**********************************'
                  WRITE(IU06,*)'*                                *'
                  WRITE(IU06,*)'* FATAL ERROR IN SUB GETSPEC     *'
                  WRITE(IU06,*)'* ===========================    *'
                  WRITE(IU06,*)'*                                *'
                  WRITE(IU06,*)'* REQUESTED DATE IS NOT EQUAL TO *'
                  WRITE(IU06,*)'* RETRIEVED DATE.                *'
                  WRITE(IU06,*)'* IN FILE: ',FILENAME
                  WRITE(IU06,*)'* CDATE = ',CDATE
                  WRITE(IU06,*)'* CDTPRO = ',CDTPRO
                  WRITE(IU06,*)'*                                *'
                  WRITE(IU06,*)'**********************************'
                  CALL ABORT1
                ENDIF

                IF(K.NE.KSEC1(44)) THEN
                  WRITE(IU06,*) '************************************'
                  WRITE(IU06,*) '* FATAL ERROR IN SUB. GETSPEC      *'
                  WRITE(IU06,*) '* REQUESTED AND DECODED DIRECTIONAL*'
                  WRITE(IU06,*) '* INDEX ARE DIFFERENT :            *'
                  WRITE(IU06,*) '* REQUESTED : ',K 
                  WRITE(IU06,*) '* DECODED   : ',KSEC1(44) 
                  WRITE(IU06,*) '*                                  *'
                  WRITE(IU06,*) '************************************'
                  CALL ABORT1
                ENDIF

                ZTHETA = REAL(KSEC1(49+K))/KSEC1(48) 
                IF(ABS(ZTHETA-(TH(K)*DEG)).GT.1.) THEN
                  WRITE(IU06,*) '************************************'
                  WRITE(IU06,*) '* FATAL ERROR IN SUB. GETSPEC      *'
                  WRITE(IU06,*) '* REQUESTED AND DECODED DIRECTIONS *'
                  WRITE(IU06,*) '* ARE DIFFERENT :                  *'
                  WRITE(IU06,*) '* ZTHETA = ',ZTHETA
                  WRITE(IU06,*) '* TH(K)*DEG = ',TH(K)*DEG
                  WRITE(IU06,*) '*                                  *'
                  WRITE(IU06,*) '*    PROGRAM ABORTS                *'
                  WRITE(IU06,*) '*                                  *'
                  WRITE(IU06,*) '************************************'
                  CALL ABORT1
                ENDIF

                IF(M.NE.KSEC1(45)) THEN
                  WRITE(IU06,*) '************************************'
                  WRITE(IU06,*) '* FATAL ERROR IN SUB. GETSPEC      *'
                  WRITE(IU06,*) '* REQUESTED AND DECODED FREQUENCY  *'
                  WRITE(IU06,*) '* INDEX ARE DIFFERENT :            *'
                  WRITE(IU06,*) '* REQUESTED : ',M 
                  WRITE(IU06,*) '* DECODED   : ',KSEC1(45) 
                  WRITE(IU06,*) '*                                  *'
                  WRITE(IU06,*) '************************************'
                  CALL ABORT1
                ENDIF


                ZFRE = REAL(KSEC1(49+NANG+M))/KSEC1(49) 
                IF(ABS(ZFRE-FR(M)).GT.0.01*FR(M)) THEN
                  WRITE(IU06,*) '************************************'
                  WRITE(IU06,*) '* FATAL ERROR IN SUB. GETSPEC      *'
                  WRITE(IU06,*) '* REQUESTED AND DECODED FREQUENCIES*'
                  WRITE(IU06,*) '* ARE DIFFERENT :                  *'
                  WRITE(IU06,*) '* ZFRE = ',ZFRE
                  WRITE(IU06,*) '* FR(M) = ',FR(M)
                  WRITE(IU06,*) '*                                  *'
                  WRITE(IU06,*) '*    PROGRAM ABORTS                *'
                  WRITE(IU06,*) '*                                  *'
                  WRITE(IU06,*) '************************************'
                  CALL ABORT1
                ENDIF

                IF(LFRSDECODE) THEN
                RMONOP=REAL(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &                0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)

                  IF(IRGG.EQ.1 .OR. CLDOMAIN == 'm') THEN
                    KDEXN = 0
                  ELSEIF(CLDOMAIN == 's') THEN
                    KDEXN = 0
                  ELSE
                    KDEXN = (NINT((RMONOP - AMONOP ) / XDELLA))*NGX
                  ENDIF

                  IF (IRGG.EQ.1) THEN
                    DO IJ = 1,NEND(NPROC) 
                      IX = IXLG(IJ,1)
                      IY = NGY- KXLT(IJ,1) +1
                      IPOINTER(IJ)=0
                      DO ICC=1,IY-1
                        JSN=NGY-ICC+1
                        IPOINTER(IJ) = NLONRGG(JSN) + IPOINTER(IJ)
                      ENDDO
                      IPOINTER(IJ) = IPOINTER(IJ) + IX 
                    ENDDO
                  ELSE
                    DO IJ = 1,NEND(NPROC)
                      IX = IXLG(IJ,1)
                      IY = NGY- KXLT(IJ,1) +1
                      IPOINTER(IJ) = (IY-1)*NGX + IX + KDEXN
                    ENDDO
                  ENDIF
                  LFRSDECODE=.FALSE.
                ENDIF

                IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)

                DO IJ = 1,NEND(NPROC) 
                  WORK(IJ) = PSEC4(IPOINTER(IJ)) 
                ENDDO
                IF(ALLOCATED(PSEC4)) DEALLOCATE(PSEC4)

                NPROMA=NEND(NPROC)/NWAM_BLKS+1
!$OMP           PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ)
                DO JKGLO=1,NEND(NPROC),NPROMA
                  KIJS=JKGLO
                  KIJL=MIN(KIJS+NPROMA-1,NEND(NPROC))
                  DO IJ=KIJS,KIJL
                    IF (WORK(IJ) .NE. ZMISS) THEN
                      WORK(IJ) = 10.**(WORK(IJ)-ABS(PPREC))- PPEPS
                    ELSE
                      WORK(IJ) = EPSMIN
                    ENDIF
                  ENDDO
                ENDDO
!$OMP           END PARALLEL DO


                CALL GSTATS(1886,1)

              ENDIF  ! decode gribex or gribapi

            ENDIF ! end decode on KSEND

          ENDDO ALL_DECODE_PE

!         SEND THE RESPECTIVE CONTRIBUTIONS OF WORK TO EACH PE.
          CALL GSTATS(623,0)
          IREQ=0
          DO IDUM=1,ISTEP_LOCAL
            IF(NPROC.EQ.1) THEN
              KSEND=1
            ELSEIF(IDUM.LT.IREAD) THEN
              KSEND=IDUM
            ELSE
              KSEND=IDUM+1
            ENDIF
            M=(((IC-1)+IDUM-1)/NANG)+1
            K=(IC-1)+IDUM-(M-1)*NANG

            ITAG=2*NFRE*NANG+(M-1)*NANG+K

            IF(IRANK.EQ.KSEND) THEN
!             SEND TO ALL OTHER TASKS
              DO IP=1,NPROC-1
                IPROC=MOD(IRANK+IP-1,NPROC)+1
                IREQ=IREQ+1
                IST=NSTART(IPROC)
                IEND=NEND(IPROC)
                CALL MPL_SEND(WORK(IST:IEND),
     &                        KDEST=IPROC,KTAG=ITAG,
     &                        KMP_TYPE=JP_NON_BLOCKING_STANDARD,
     &                        KREQUEST=ISENDREQ(IREQ),
     &                        CDSTRING='GETSPEC: SENDING WORK' )
              ENDDO

            ENDIF
          ENDDO

!         RECEIVE THE RESPECTIVE CONTRIBUTIONS OF WORK TO EACH PE.
          IPROC=IRANK
          IST=NSTART(IPROC)
          IEND=NEND(IPROC)
          ALLOCATE(ZRECVBUF(IST:IEND))
          DO IDUM=1,ISTEP_LOCAL
            IF(NPROC.EQ.1) THEN
              KSEND=1
            ELSEIF(IDUM.LT.IREAD) THEN
              KSEND=IDUM
            ELSE
              KSEND=IDUM+1
            ENDIF

            M=(((IC-1)+IDUM-1)/NANG)+1
            K=(IC-1)+IDUM-(M-1)*NANG

            IF(IRANK.EQ.KSEND) THEN
!             SAVE LOCAL CONTRIBUTION
              DO IJ = NSTART(IRANK),NEND(IRANK) 
                FL(IJ,K,M)=WORK(IJ)
              ENDDO
            ELSE
!             RECEIVE INFORMATION FROM KSEND (that sets M and K)
              CALL MPL_RECV(ZRECVBUF(IST:IEND),
     &                      KFROM=KFROM,KRECVTAG=KRTAG,
     &                      KMP_TYPE=JP_BLOCKING_STANDARD,
     &                      CDSTRING='GETSPEC: RECEIVING WORK' )

              IF(KFROM.LT.IREAD) THEN
                ID=KFROM
              ELSE
                ID=KFROM-1
              ENDIF
              MR=(((IC-1)+ID-1)/NANG)+1
              KR=(IC-1)+ID-(MR-1)*NANG
              ITAG=2*NFRE*NANG+(MR-1)*NANG+KR
              IF(KRTAG.NE.ITAG) THEN
                WRITE(0,*)'MPL_RECV ERROR in GETSPEC: MISMATCHED TAGS'
                WRITE(0,*)'IRANK = ',IRANK
                WRITE(0,*)'KFROM = ',KFROM
                WRITE(0,*)'KRTAG = ',KRTAG
                WRITE(0,*)'ITAG  = ',ITAG
                WRITE(0,*)' ABORTING !!!!'
                CALL ABORT1
              ENDIF

              DO IJ = IST,IEND 
                FL(IJ,KR,MR)=ZRECVBUF(IJ) 
              ENDDO

            ENDIF

            FL(NINF-1,K,M)=0.
          ENDDO
          DEALLOCATE(ZRECVBUF)

!         ENSURE ALL SENDS ARE FINISHED.
          IF(IREQ.GT.0) THEN
            CALL MPL_WAIT(WORK(:),KREQUEST=ISENDREQ(1:IREQ),
     &                    CDSTRING='GETSPEC: WAIT SENDING WORK')
          ENDIF
          CALL GSTATS(623,1)


!         MAKE SURE THAT ALL RECEIVE ARE FINISHED ON ALL PE'S
!         BEFORE PROCESSING ANOTHER BATCH.
          CALL MPL_BARRIER(CDSTRING='GETSPEC:')

        ENDDO ALL_FILE

        IF(LGRIBAPI_I) THEN
          IF (IRANK.EQ.IREAD) THEN
            CALL IGRIB_CLOSE_FILE(KFILE_HANDLE)
          ELSE
          ENDIF

          IF(ALLOCATED(WORK)) DEALLOCATE(WORK)
        ELSE

          IF (IRANK.EQ.IREAD) THEN
            CALL PBCLOSE(IUGRS,KRET)
            IF(KRET.LT.0) THEN 
              WRITE (IU06,*) '**********************************'
              WRITE (IU06,*) '*                                *'
              WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE '
              WRITE (IU06,*) '* IN GETSPEC                     *'
              WRITE (IU06,*) '* FILE ',FILENAME
              WRITE (IU06,*) '*                                *'
              WRITE (IU06,*) '**********************************'
              CALL ABORT1
            ENDIF
          ENDIF

          IF(ALLOCATED(KSEC0)) DEALLOCATE(KSEC0)
          IF(ALLOCATED(KSEC1)) DEALLOCATE(KSEC1)
          IF(ALLOCATED(KSEC2)) DEALLOCATE(KSEC2)
          IF(ALLOCATED(KSEC3)) DEALLOCATE(KSEC3)
          IF(ALLOCATED(KSEC4)) DEALLOCATE(KSEC4)
          IF(ALLOCATED(PSEC2)) DEALLOCATE(PSEC2)
          IF(ALLOCATED(PSEC3)) DEALLOCATE(PSEC3)
          IF(ALLOCATED(IPOINTER)) DEALLOCATE(IPOINTER)
          IF(ALLOCATED(WORK)) DEALLOCATE(WORK)

        ENDIF

!       THE INPUT SPECTRA ARE IN THE ABSOLUTE FRAME OF REFERENCE
!       THEY NEED TO BE CONVERTED TO THE RELATIVE FRAME
!       (WITH RESPECT TO AMBIENT CURRENTS)
        IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
          IG=1
          IRA=-1

          NPROMA=(NEND(IRANK)-NSTART(IRANK)+1)/NWAM_BLKS+1

          CALL GSTATS(1236,0)
!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1) 
!$OMP+    PRIVATE(JKGLO,KIJS,KIJL,K,M,IJ)
          DO JKGLO=NSTART(IRANK),NEND(IRANK),NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,NEND(IRANK))

            CALL INTPOL (FL, SL, KIJS, KIJL, IG, IRA)

            DO M=1,NFRE
              DO K=1,NANG
                DO IJ=KIJS,KIJL
                  FL(IJ,K,M) = SL(IJ,K,M)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1236,1)

        ENDIF


      ELSE
 
!     BINARY INPUT:
!     =============
         IF (ITEST.GE.2) THEN
           WRITE(IU06,*)''
           WRITE(IU06,*)' IN GETSPEC :'
           IF(MDEL.EQ.1) THEN
             WRITE(IU06,*)'  ',MDEL,
     &                  ' FREQUENCY IS DISTRIBUTED AT ONE TIME'
           ELSE
             WRITE(IU06,*)'  ',MDEL,
     &                    ' FREQUENCIES ARE DISTRIBUTED EACH TIME'
           ENDIF
         ENDIF

         IF (IRANK.EQ.ISEND)
     &       CALL GRSTNAME(CDTPRO,CDATEF,'BLS',ICPLEN,CPATH,FILENAME)

         IUNIT=0

         DO MLOOP=1,NFRE,MDEL
           MINF=MLOOP
           MSUP=MIN(MLOOP+MDEL-1,NFRE)
           DO KLOOP=1,NANG,KDEL
             KINF=KLOOP
             KSUP=MIN(KLOOP+KDEL-1,NANG)

             ALLOCATE(RFL(0:NIBLO,KINF:KSUP,MINF:MSUP))
!            READ RESTART SPECTRA FROM PE ISEND (IREAD) 
             IF (IRANK.EQ.ISEND) THEN
               LPBOPEN = .FALSE.
               LPBCLOSE = .FALSE.
               IF(MINF.EQ.1 .AND. KINF.EQ.1) LPBOPEN = .TRUE.
               IF(MSUP.EQ.NFRE .AND. KSUP.EQ.NANG) LPBCLOSE = .TRUE.

               CALL READFL(RFL, KINF, KSUP, MINF, MSUP,
     &                     FILENAME, IUNIT, LPBOPEN, LPBCLOSE)
             ENDIF

             CALL MPDISTRIBFL(ISEND,KTAG,NSTART,NEND,KINF,KSUP,
     &                        MINF,MSUP,RFL)
             KTAG=KTAG+1

             IF (ITEST.GE.2)
     &        WRITE(IU06,*)
     &        'SUB. GETSPEC: RESTART SPECTRUM COLLECTED, PART:',MLOOP

!            KEEP CORRESPONDING CONTRIBUTION TO FL
             DO M=MINF,MSUP
               DO K=KINF,KSUP
                 DO IJ=NSTART(IRANK),NEND(IRANK)
                   FL(IJ,K,M) = RFL(IJ,K,M)
                 ENDDO
               ENDDO
             ENDDO

             DO M=MINF,MSUP
               DO K=KINF,KSUP
                 FL(NINF-1,K,M) = 0. 
               ENDDO
             ENDDO

             DEALLOCATE(RFL)
           ENDDO
         ENDDO

      ENDIF

      WRITE(IU06,*) ' SPECTRUM FILE READ IN...............',
     & ' CDTPRO  = ', CDTPRO
      WRITE(IU06,*) ' '
      CALL FLUSH (IU06)

      RETURN
      END SUBROUTINE GETSPEC
