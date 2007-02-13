      SUBROUTINE GETSPEC(FL,NSTART,NEND,IREAD,U10OLD,THWOLD)
! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      SEPTEMBER 1997 

!*    PURPOSE.
!     --------
!     READ THE SPECTRA FROM DISK.

!**   INTERFACE.
!     ----------
!     *CALL* *GETSPEC(FL,NSTART,NEND,IREAD,U10OLD,THWOLD)
!     *FL*        ARRAY CONTAINING THE SPECTRA CONTRIBUTION ON EACH PE
!     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *IREAD*     PROCESSOR WHICH WILL ACCESS THE FILE ON DISK 
!     *U10OLD*    WIND SPEED. (used with fetch law to fill empty 
!                 sea points)
!     *THWOLD*    WIND DIRECTION (RADIANS).


!     METHOD.
!     -------

!     IN CASE THE INPUT SPECTRA ARE IN GRIB THEN EITHER THEY ARE READ
!     FROM A FILE ON DISK USING THE PBGRIB SOFTWARE OR THEY ARE READ
!     DIRECTLY FROM FDB. GRIB SPECTRA WILL BE DECODED. IN CASE OF BINARY
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

      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3  , 
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK3    ,PSEK3   ,
     &            NTENCODE ,PPEPS    ,PPREC
      USE YOWGRID  , ONLY : DELPHI   ,NLONRGG
      USE YOWJONS  , ONLY : AJONS    ,BJONS    ,DJONS    ,EJONS
      USE YOWICE   , ONLY : FLMINFR  ,ICEMASK
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IRGG     ,AMOSOP   ,
     &            AMONOP   ,XDELLA
      USE YOWMESPAS, ONLY : LMESSPASS,LPBIOIN  ,LFDBIOIN ,LGRIBIN
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,
     &            KTAG     ,NPRECR   ,NPRECI
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,
     &            NBLO     ,NIBLO    ,CLDOMAIN
      USE YOWPCONS , ONLY : G        ,DEG      ,R        ,ZMISS    ,
     &            EPSMIN
      USE YOWSTAT  , ONLY : CDATEF   ,CDTPRO   ,CFDB2DSP ,YCLASS   ,
     &            YEXPVER
      USE YOWTABL  , ONLY : JUMAX    ,DELU
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH    ,LRESTARTED
      USE MPL_MODULE
      USE FDBSUBS

! ----------------------------------------------------------------------

      PARAMETER (NBIT = 100000)

      CHARACTER*  9 CMBXSIZE
      CHARACTER* 12 CDATE 
      CHARACTER* 120 FILENAME

!* fdb variables.

      CHARACTER*  4 CLEXPVER
      CHARACTER*  2 CLTYPE, CLCLASS, CLDIR, CLFRE

      LOGICAL LFRSTIME, LPBOPEN, LPBCLOSE, LLEXIST

      INTEGER ifdb
      SAVE ifdb

! ----------------------------------------------------------------------

      INTEGER,DIMENSION(NPROC) :: NSTART, NEND
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL
      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD

! ----------------------------------------------------------------------
!     define local allocatable arrays

      INTEGER,ALLOCATABLE :: INGRIB(:), IPOINTER(:)
      INTEGER,ALLOCATABLE :: KSEC0(:),KSEC1(:),KSEC2(:),KSEC3(:)
      INTEGER,ALLOCATABLE :: KSEC4(:)
      REAL, ALLOCATABLE,DIMENSION(:) :: PSEC2,PSEC3,PSEC4,WORK,WORKRCV

      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: RFL
      REAL,ALLOCATABLE :: EM(:), ET(:,:), STK(:)
      REAL :: ONE,ZTHETA,ZFRE

      INTEGER :: ISENDREQ

      LOGICAL :: LLEPSMIN

      DATA MDEL, LFRSTIME, ONE / 1, .TRUE., 1.0/

! ----------------------------------------------------------------------

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
      ONE=1.0
      ISEND=IREAD

      IF(LGRIBIN.AND..NOT.LRESTARTED) THEN

!       INPUT SPECTRA ARE IN GRIB FORMAT

        ALLOCATE(KSEC0(NKSEK0),KSEC1(NKSEK1))
        ALLOCATE(KSEC2(NKSEK2),KSEC3(NKSEK3),KSEC4(NKSEK4))
        ALLOCATE(PSEC2(NPSEK2),PSEC3(NPSEK3))
        ALLOCATE(IPOINTER(1:NEND(NPROC)))
        ALLOCATE(WORK(NIBLO))

        IF (IRANK.EQ.IREAD) THEN
          IF (.NOT.LFDBIOIN) THEN

!           CONNECT INPUT PE (IREAD) WITH INPUT FILE

            CALL GSTATS(1888,0)
            FILENAME='specwavein'
            LFILE=LEN_TRIM(FILENAME)
            INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
            CALL GSTATS(1888,1)
            IF(LLEXIST) THEN

              CALL PBOPEN(IUGRS,FILENAME(1:LFILE),'r',KRET)
              IF(KRET.LT.0) THEN 
                WRITE (IU06,*) '***********************************'
                WRITE (IU06,*) '*                                 *'
                WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBOPEN  *'
                WRITE (IU06,*) '* IN GETSPEC                      *'
                IF(KRET.EQ.-1)
     &             WRITE (IU06,*) ' COULD NOT OPEN FILE ',FILENAME
                IF(KRET.EQ.-2)
     &           WRITE (IU06,*) ' INVALID FILENAME ',FILENAME
                IF(KRET.EQ.-3)
     &           WRITE (IU06,*)' INVALID OPEN MODE SPECIFIED'
                WRITE (IU06,*) '*                                 *'
                WRITE (IU06,*) '***********************************'
                CALL ABORT1
              ENDIF
            ELSE
              WRITE(IU06,*)'**************************************'
              WRITE(IU06,*)'*                                    *'
              WRITE(IU06,*)'*GETSPEC : GRIB SPECTRA NOT FOUND IN *'
              WRITE(IU06,*)  FILENAME
              WRITE(IU06,*)'*PROGRAM WILL ABORT                  *'
              WRITE(IU06,*)'*                                    *'
              WRITE(IU06,*)'**************************************'
              CALL ABORT1
            ENDIF

          ELSE

!           INITIALIZE AND OPEN FDB.

            istat = iopenfdbsubs ('fdb', ifdb, 'r')
            IF (ITEST.GE.2)WRITE(iu06,'("  FDB OPENED status=",I3,
     &         " Reference=", I2)') istat, ifdb

            ICFDB2DSP = LEN_TRIM(CFDB2DSP)
            IF (ICFDB2DSP .GT. 0 ) THEN
              ISTAT = ISET_FDBSUBS_ROOT(IFDB, CFDB2DSP(1:ICFDB2DSP) )
              IF ( ISTAT .NE. 0 ) THEN
                WRITE(IU06,'("Error\ /Getspec/ Root dir specified:
     &                       ",A)') CFDB2DSP(1:ICFDB2DSP)
                WRITE(IU06,'("Error\ /Getspec/
     &                      iopenfdb return status:",I3)') istat
                kerr = 1
                CALL ABORT1
              ELSEIF (ITEST.GE.2) THEN
                WRITE(IU06,'("\ /Getspec/ iset_fdb_root status:",
     &                  i3,
     &          " Root: ", a25)') istat, CFDB2DSP(1:ICFDB2DSP)
              ENDIF
            ENDIF

!           SET FIXED PARAMETERS FOR FDB.

            istat = isetvalfdbsubs (ifdb, 'levty',  's'      )
            istat = isetvalfdbsubs (ifdb, 'param',  '251'    )
            istat = isetvalfdbsubs (ifdb, 'repres', 'll'     )
            istat = isetvalfdbsubs (ifdb, 'stream', 'wave'   )
            istat = isetvalfdbsubs (ifdb, 'step',   '000000' )
!---<       istat = isetvalfdbsubs (ifdb, 'fstep',  '000000' )

            kfcst = 0
            IF ( kfcst .EQ. 0 ) THEN
              cltype='an'
            ELSE
              cltype='fc'
            ENDIF

            clclass = yclass
            clexpver = yexpver

            CALL u2l1cr( cltype )
            CALL u2l1cr( clclass )
            CALL u2l1cr( clexpver )

            istat = isetvalfdbsubs (ifdb, 'expver', clexpver    )
            istat = isetvalfdbsubs (ifdb, 'date',   cdtpro(1:8))
            istat = isetvalfdbsubs (ifdb, 'type',   cltype     )
            istat = isetvalfdbsubs (ifdb, 'time',   cdtpro(9:10))
            istat = isetvalfdbsubs (ifdb, 'domain', cldomain )
!---<       istat = isetrankfdbsubs (ifdb, ipe)

            IF (ITEST.GE.2) THEN
              WRITE(IU06,'("\ /getspec/ PARAM: ", a8)') '251'
              WRITE(IU06,'("\ /getspec/ STEP:  ", a8)') '000000'
              WRITE(IU06,'("\ /getspec/ CLASS: ", a8)') clclass
              WRITE(IU06,'("\ /getspec/ EXPVER:", a8)') clexpver
              WRITE(IU06,'("\ /getspec/ DOMAIN:", a8)') cldomain
              WRITE(IU06,'("\ /getspec/ TY:   ", a8)') cltype
              WRITE(IU06,'("\ /getspec/ DATE: ", a8)') cdtpro(1:8)
              WRITE(IU06,'("\ /getspec/ TIME: ", a8)') cdtpro(9:10)
            ENDIF
          ENDIF
        ENDIF


!       GET GRIB DATA FROM (NFRE*NANG) FIELDS


!       READ NPROC-1 FIELDS AND SEND THEM SUCCESSIVELY TO ALL OTHER PE'S
!       FOR DECODING (IF IN MESSAGE PASSING MODE)

        ISIZE=NBIT
        ISTEP=MAX(NPROC-1,1)

        ALL_FILE : DO IC=1,NFRE*NANG,ISTEP 

          IST=ISTEP
          IF(IC+ISTEP.GT.NFRE*NANG ) IST=NFRE*NANG-IC+1

          ALL_DECODE_PE : DO IDUM=1,IST
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

            IF (LFDBIOIN.AND.IRANK.EQ.IREAD) THEN

!             READ GRIB DATA FROM FDB (IF REQUESTED)
!             SET VARIABLE PARAMETERS FOR FDB.

              WRITE( cldir, '(i2.2)' ) K 
              istat = isetvalfdbsubs( ifdb, 'direction', cldir )
              WRITE( clfre, '(i2.2)' ) M
              istat = isetvalfdbsubs( ifdb, 'frequency', clfre )
              istat = isetvalfdbsubs (ifdb, 'class',  clclass       )
              IF (ITEST.GE.2) THEN
                WRITE(IU06,'("\ /getspec/ DIR:   ", a8)') cldir
                WRITE(IU06,'("\ /getspec/ FREQ:  ", a8)') clfre
              ENDIF

!* try to estimate the size of the next product:

              icount = 1129120
              IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)
              ALLOCATE (ingrib(icount))
              IF (ITEST.GE.2) WRITE(IU06,'(" GETSPEC: allocate", I10,
     &          " words for ingrib", " nstart=", i8, " nend=", i8,
     &          " spectrum=" , i4) ')icount,nstart(irank),nend(irank),
     &          nang*nfre
              irstat = ireadfdbsubs (ifdb, ingrib, icount      )
              IF (irstat.NE.0 .OR. icount.LE.0) THEN
                WRITE (IU06,'("***********************************")')
                WRITE (IU06,'("*                                 *")')
                WRITE (IU06,'("*ERROR FOLLOWING CALL TO ireadfdb *")')
                WRITE (IU06,'("*IN GETSPEC                       *")')
                WRITE (IU06,'("*STAT CODE irstat =", I10)') irstat
                WRITE (IU06,'("*LENGTH OF MESSAGE=", I10)') icount
                WRITE (IU06,'("*                                 *")')
                WRITE (IU06,'("***********************************")')
                CALL ABORT1
              ELSEIF (ITEST.GE.2) THEN
                WRITE(IU06,'(" GETSPEC: readfdb read ",i10,
     &          " words (status:", i3, ")" )')icount, irstat
              ENDIF

              IF(M.EQ.NFRE.AND.K.EQ.NANG) THEN 
                CALL GSTATS(1787,0)
                istat = iclosefdbsubs(ifdb)
                CALL GSTATS(1787,1)
                IF (ISTAT .NE. 0) CALL ABORT1
              ENDIF

              isize = icount
              IF (ITEST.GE.2) CALL FLUSH (IU06)


            ELSEIF (.NOT.LFDBIOIN.AND.IRANK.EQ.IREAD) THEN

!             READ GRIB DATA FROM FILE (IF REQUESTED)

1000          IPLENG=ISIZE*NPRECI
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
                WRITE (IU06,*) ' FILE HANDLING ERROR ON UNIT ',IUGRS
                WRITE (IU06,*) ' PROGRAM WILL ABORT '
                CALL ABORT1
              ELSEIF (KRET .EQ. -3) THEN
                DEALLOCATE(INGRIB)
                CALL KGRIBSIZE(IU06, IUGRS, ILENG, ISIZE, 'GETSPEC')
                GOTO 1000
              ENDIF

            ENDIF


            IF(LMESSPASS) THEN

!             IN CASE OF MESSAGE PASSING THE DECODING OCCURS ON KSEND

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

!               RECEIVE GRIB DATA SIZE FROM PE IREAD

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
     &           KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ,
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
                CALL MPL_WAIT(INGRIB(1:ISIZE),KREQUEST=ISENDREQ, 
     &           CDSTRING='GETSPEC: WAIT FOR SEND')

                IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)

              ENDIF

              CALL GSTATS(623,1)
            ENDIF

            IF (IRANK.EQ.KSEND.OR.NPROC.EQ.1) THEN

!             DECODE THE GRIB DATA ON PE KSEND 

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
              WRITE (CDATE,'(I4.4,4I2.2)') IYYYY,KSEC1(11),KSEC1(12),
     &                                     KSEC1(13),KSEC1(14)

              IF(KSEC1(39).EQ.9.OR.KSEC1(39).EQ.6.OR.
     &           KSEC1(39).EQ.10.OR.KSEC1(39).EQ.11)
     &        THEN
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
!               DETERMINE FORECAST PERIOD.
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

              IF(LFRSTIME) THEN
                RMONOP=REAL(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &               0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)

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
                LFRSTIME=.FALSE.
              ENDIF

              IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)

              DO IJ = 1,NEND(NPROC) 
                WORK(IJ) = PSEC4(IPOINTER(IJ)) 
              ENDDO
              IF(ALLOCATED(PSEC4)) DEALLOCATE(PSEC4)

              DO IJ = 1,NEND(NPROC) 
                IF (WORK(IJ) .NE. ZMISS) THEN
                  WORK(IJ) = 10.**(WORK(IJ)-ABS(PPREC))- PPEPS
                ELSE
                  WORK(IJ) = EPSMIN
                ENDIF
              ENDDO
              CALL GSTATS(1886,1)

            ENDIF

          ENDDO ALL_DECODE_PE

!         DISTRIBUTE THE RESPECTIVE CONTRIBUTION OF WORK TO EACH PE.

          ALLOCATE(WORKRCV(NIBLO))

          DO IDUM=1,IST
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
              CALL MPDISTRIBSCFLD(KSEND,ITAG,NSTART,NEND,WORK)
              DO IJ = NSTART(IRANK),NEND(IRANK) 
                FL(IJ,K,M) = WORK(IJ)
              ENDDO
            ELSE
              CALL MPDISTRIBSCFLD(KSEND,ITAG,NSTART,NEND,WORKRCV)
              DO IJ = NSTART(IRANK),NEND(IRANK) 
                FL(IJ,K,M) = WORKRCV(IJ)
              ENDDO
            ENDIF

            FL(NINF-1,K,M) = 0.
          ENDDO

          IF(ALLOCATED(WORKRCV)) DEALLOCATE(WORKRCV)

        ENDDO ALL_FILE

        IF (.NOT.LFDBIOIN.AND.IRANK.EQ.IREAD) THEN
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


      ELSE

        IF(LMESSPASS) THEN

          IF(LFRSTIME) THEN

!         DETERMINE MAXIMUM NUMBER OF FREQUENCIES WHICH CAN BE USED
!         TO SEND SPECTRA FROM ONE PE  TO THE OTHERS
!         WITHOUT EXCEEDING THE MAXIMUM MAILBOX SIZE 

            CALL GETENV('MBX_SIZE',CMBXSIZE)
            READ(CMBXSIZE,'(I9)') MBXSIZE
            MINSIZE=NINT(REAL(NPROC-1)*NIBLO*NANG*NPRECR/NPROC)
            IF(MINSIZE.GT.MBXSIZE) THEN
              WRITE (IU06,*) '***************************************'
              WRITE (IU06,*) '*                                     *'
              WRITE (IU06,*) '* MESSAGE MAIL BOX SIZE IS TOO SMALL  *'
              WRITE (IU06,*) '* TO DISTRIBUTE THE RESTART SPECTRUM  *'
              WRITE (IU06,*) '* FILE                                *'
              WRITE (IU06,*) '* THE PRESCRIBED SIZE SIZE IS         *'
              WRITE (IU06,*) '* ', MBXSIZE,' BYTES                  *'
              WRITE (IU06,*) '* PLEASE INCREASE SIZE TO AT LEAST    *'
              WRITE (IU06,*) '* ', MINSIZE, 'BYTES                  *'
              WRITE (IU06,*) '* SEE ENVIRONMENT VARIABLE            *'
              WRITE (IU06,*) '* MBX_SIZE                            *' 
              WRITE (IU06,*) '*                                     *'
              WRITE (IU06,*) '***************************************'
              CALL ABORT1
            ENDIF

            MDEL=1
!CC          MDEL = NFRE 
!CC          DO WHILE (MINSIZE*MDEL.GE.MBXSIZE.AND.MDEL.GT.1)
!CC             MDEL=MDEL-1
!CC          ENDDO

            IF (ITEST.GE.2) THEN
              WRITE(IU06,*)''
              WRITE(IU06,*)' IN GETSPEC :'
              WRITE(IU06,*)'  MESSAGE MAILBOX SIZE OF',MBXSIZE,' BYTES'
              IF(MDEL.EQ.1) THEN
                WRITE(IU06,*)'  ',MDEL,
     &                     ' FREQUENCY IS DISTRIBUTED AT ONE TIME'
              ELSE
                WRITE(IU06,*)'  ',MDEL,
     &                       ' FREQUENCIES ARE DISTRIBUTED EACH TIME'
              ENDIF
              WRITE(IU06,*)'  WHICH SHOULD RESULT IN A MAXIMUM USAGE OF'
              WRITE(IU06,*)'  ',MINSIZE*MDEL,
     &                     ' BYTES FOR EACH MPDISTRIBFL CALL'
              WRITE(IU06,*)''
            ENDIF
          ENDIF


          IF(LPBIOIN) THEN

!           INPUT SPECTRA ARE IN A BINARY FILE AND WILL BE READ
!           USING THE PBIO SOFTWARE

            IF (IRANK.EQ.ISEND)
     &          CALL GRSTNAME(CDTPRO,CDATEF,'BLS',ICPLEN,CPATH,FILENAME)

            DO MLOOP=1,NFRE,MDEL

              MINF=MLOOP
              MSUP=MIN(MLOOP+MDEL-1,NFRE)
              ALLOCATE(RFL(0:NIBLO,NANG,MINF:MSUP))
!              READ RESTART SPECTRA FROM PE ISEND (IREAD) 
              IF (IRANK.EQ.ISEND) THEN

                LPBOPEN = .FALSE.
                LPBCLOSE = .FALSE.
                IF(MINF.EQ.1) LPBOPEN = .TRUE.
                IF(MSUP.EQ.NFRE) LPBCLOSE = .TRUE.

                CALL READFL(RFL,MINF,MSUP,FILENAME,IUNIT,
     &           LPBOPEN,LPBCLOSE)
              ENDIF
!              PUT A BARRIER TO AVOID FLOODING THE CROSSBAR 
              CALL MPL_BARRIER(CDSTRING='GETSPEC:')

              CALL MPDISTRIBFL(ISEND,KTAG,NSTART,NEND,MINF,MSUP,
     &         MINF,MSUP,RFL)
              KTAG=KTAG+1
              IF (ITEST.GE.2)
     &         WRITE(IU06,*)
     &         'SUB. GETSPEC: RESTART SPECTRUM COLLECTED, PART:',MLOOP

!              KEEP CORRESPONDING CONTRIBUTION TO FL
              DO M=MINF,MSUP
                DO K=1,NANG
                  DO IJ=NSTART(IRANK),NEND(IRANK)
                    FL(IJ,K,M) = RFL(IJ,K,M)
                  ENDDO
                ENDDO
              ENDDO

              DO M=MINF,MSUP
                DO K=1,NANG
                  FL(NINF-1,K,M) = 0. 
                ENDDO
              ENDDO

              IF(ALLOCATED(RFL)) DEALLOCATE(RFL)

            ENDDO

          ELSE

!           BINARY SPECTRA ARE READ USING FORTRAN READ

            ALLOCATE(RFL(0:NIBLO,NANG,NFRE))

!           READ RESTART SPECTRA FROM PE ISEND (IREAD)
            IF (IRANK.EQ.ISEND) THEN
              CALL GRSTNAME(CDTPRO,CDATEF,'BLS',ICPLEN,CPATH,FILENAME)
              IF (ITEST > 3 ) THEN
                 WRITE(IU06,*) ' BLS arranged to ', FILENAME
                 CALL FLUSH (IU06)
              ENDIF
              CALL READFL(RFL,1,NFRE,FILENAME,IUNIT,LPBOPEN,LPBCLOSE)
            ENDIF

            DO MLOOP=1,NFRE,MDEL
              MINF=MLOOP
              MSUP=MIN(MLOOP+MDEL-1,NFRE)
!             PUT A BARRIER TO AVOID FLOODING THE CROSSBAR 
              CALL MPL_BARRIER(CDSTRING='GETSPEC 2:')

              CALL MPDISTRIBFL(ISEND,KTAG,NSTART,NEND,MINF,MSUP,
     &                         1,NFRE,RFL)
              KTAG=KTAG+1
            ENDDO

!           KEEP CORRESPONDING CONTRIBUTION TO FL
            DO K=1,NANG
              DO M=1,NFRE
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  FL(IJ,K,M) = RFL(IJ,K,M)
                ENDDO
              ENDDO
            ENDDO

            DO K=1,NANG
              DO M=1,NFRE
                FL(NINF-1,K,M) = 0. 
              ENDDO
            ENDDO

            IF (ITEST.GE.2)
     &      WRITE(IU06,*)
     &      'SUB. GETSPEC: RESTART SPECTRUM COLLECTED'

            IF(ALLOCATED(RFL)) DEALLOCATE(RFL)

          ENDIF

        ELSE

!         NON MESSAGE PASSING READ OF SPECTRA

          CALL GRSTNAME(CDTPRO,CDATEF,'BLS',ICPLEN,CPATH,FILENAME)
          LPBOPEN = .TRUE.
          LPBCLOSE = .TRUE.

          CALL READFL(FL,1,NFRE,FILENAME,IUNIT,LPBOPEN,LPBCLOSE)
        ENDIF

      ENDIF

      IF(.NOT.LRESTARTED .AND. CLDOMAIN.NE.'s') THEN
!       TRY TO CREATE A JONSWAP SPECTRUM AT POINTS WHICH WERE ICE IN THE
!       PREVIOUS RUN AND WHICH ARE NOW OPEN SEA.
!       THEY ARE CHARCTERISED BY HAVING NO WAVE ENERGY BESIDE NOISE
!       AND HAVE AN ICEMASK = 1.
        DO K=1,NANG
          DO M=1,NFRE
            DO IJ=NSTART(IRANK),NEND(IRANK)
               IF(FL(IJ,K,M) .LE. EPSMIN) FL(IJ,K,M) = 0. 
            ENDDO
          ENDDO
        ENDDO

        LLEPSMIN=.FALSE.

        ALLOCATE(EM(NSTART(IRANK):NEND(IRANK)))
        ALLOCATE(ET(1,NFRE))
        ALLOCATE(STK(NANG))

        CALL SEMEAN(FL,NSTART(IRANK),NEND(IRANK),EM,LLEPSMIN)

        GAMMA=3.000000
        SA=7.000000E-02
        SB=9.000000E-02
        FETCH=MIN(5*DELPHI,250000.)
        GX = G * FETCH
        FPMAX=FR(NFRE-1)

        DO IJ=NSTART(IRANK),NEND(IRANK)
           IF(EM(IJ).LE.EPSMIN .AND. ICEMASK(IJ,1).EQ.1 ) THEN
!            find peak period and energy level from fetch law
!            the same formulation as in subroutine peak is used.
             IF (U10OLD(IJ,1) .GT. 0.1E-08) THEN
               U10 = U10OLD(IJ,1)
               GXU = GX/(U10*U10)
               UG = G/U10
               FPK = AJONS * GXU ** DJONS
               FPK = MAX(0.13, FPK)
               FPK = MIN(FPK, FPMAX/UG)
               ALPHJO = BJONS * FPK** EJONS
               ALPHJO = MAX(ALPHJO, 0.0081)
               FPK = FPK*UG
             ELSE
               FPK = 0.
               ALPHJO = 0.
             ENDIF

             CALL JONSWAP(ALPHJO, GAMMA, SA, SB, FPK, 1, 1, ET)

             THESK=THWOLD(IJ,1)

             CALL SPR (NANG, NANG, THESK, TH, STK)

             DO K=1,NANG
               DO M=1,NFRE
                 FL(IJ,K,M) = ET(1,M)*STK(K)
               ENDDO
             ENDDO

           ENDIF
        ENDDO

        DO K=1,NANG
          DO M=1,NFRE
            DO IJ=NSTART(IRANK),NEND(IRANK)
               XJ=U10OLD(IJ,1)/DELU
               JU=MIN(JUMAX, MAX(NINT(XJ),1))
               FLLOWEST=FLMINFR(JU,M)*MAX(0.,COS(TH(K)-THWOLD(IJ,1)))**2
               IF(FL(IJ,K,M) .LE. FLLOWEST) FL(IJ,K,M) = FLLOWEST
            ENDDO
          ENDDO
        ENDDO

        DEALLOCATE(EM)
        DEALLOCATE(ET)
        DEALLOCATE(STK)

      ENDIF


      WRITE(IU06,*) ' SPECTRUM FILE READ IN...............',
     & ' CDTPRO  = ', CDTPRO
      WRITE(IU06,*) ' '
      CALL FLUSH (IU06)

      RETURN
      END SUBROUTINE GETSPEC
