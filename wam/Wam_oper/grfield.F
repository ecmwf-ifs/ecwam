      SUBROUTINE GRFIELD(WHMOD,SIGMOD,DISTMAX)

!--------------------------------------------------------------------

!**** *GRFIELD* - READS GRIDDED MEASUREMENTS AND DISTRIBUTES THEM.
!                 ALSO COMPUTES THE ALTIMETER DATA CORRECTIONS.

!     B.HANSEN       ECMWF       JANUARY 1998
!     P.JANSSEN      ECMWF       FEBRUARY 1999 ADD SEA STATE DEPENDENT
!                                              CORRECTION TO ALTIMETER
!                                              WAVE HEIGHT.

!     PURPOSE.
!     --------

!       PREPARE THE DATA FOR OPTIMAL INTERPOLATION.

!**   INTERFACE.
!     ----------

!       *CALL* *GRFIELD(WHMOD,SIGMOD,DISTMAX)*

!          *WHMOD*      REAL     MODEL WAVE HEIGHT (FIRST GUESS).
!          *SIGMOD*     REAL     MODEL ERROR. 
!          *DISTMAX*    REAL     MAXIMUM SPREADING DISTANCE IN RADIAN !!

!     METHOD.
!     -------

!        READ PRE PREPARED SATELLITE DATA FROM FILE AND DISTRIBUTE TO
!        OTHER PES.

!      EXTERNALS.
!      ----------

!         MPBCASTGRDFLD  - BROADCAST 2-D DATA FIELD FROM 1 PE TO ALL
!                          OTHERS.
!         SKEWNESS       - COMPUTES PARAMETERS OF THE NEARLY-GAUSSIAN
!                          DISTRIBUTION OF OCEAN WAVES

!      MODIFICATONS.
!      -------------
!         NONE.

!-----------------------------------------------------------------------

      USE YOWALTAS , ONLY : NUMALT   ,NIJALT   ,NALTDT   ,NALTAVLB ,
     &            IJALT    ,INTLMAX  ,KMINLMAX ,KMAXLMAX ,NOBSPE   ,
     &           ALTSDTHRSH,ALTBGTHRSH,HSALTCUT,ALTDATA  ,IBUFRSAT ,
     &            XKAPPA2  ,HSCOEFCOR,HSCONSCOR ,LALTCOR ,LALTLRGR ,
     &            LALTGRDOUT
      USE YOWCOUT  , ONLY : IPFGTBL
      USE YOWGRIBHD, ONLY : KSEK1
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL      ,SINPH    ,
     &            COSPH
      USE YOWICE   , ONLY : ICEMASK  ,FLMIN
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,AMOWEP   ,AMONOP   ,
     &            XDELLA   ,ZDELLO
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,KTAG
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,
     &            NIBLO    ,CLDOMAIN ,LL1D
      USE YOWPCONS , ONLY : RAD      ,ZMISS
      USE YOWSTAT  , ONLY : CDTPRO   ,YEXPVER
      USE YOWSPEC  , ONLY : NSTART   ,NEND     ,FL1      ,IJ2NEWIJ
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWINTP  , ONLY : WHGTTG   ,WHGTAG   ,CWHTAG   , RANGAG
      USE MPL_MODULE

! ----------------------------------------------------------------------

      LOGICAL LLOPENED, LLTEST
      INTEGER IREAD, IDUM
      INTEGER,DIMENSION(NPROC) :: ISENDLEN
      INTEGER,ALLOCATABLE :: ISENDREQ(:)
      CHARACTER*2 :: CDUM
      CHARACTER* 3 CLFID
      CHARACTER*4 CSTRM
      CHARACTER*12 CTIME
      REAL :: DISTMAX, COSDISTMAX
      REAL :: DELLON, COSLON, DIS2, XLONGKJ
      REAL, ALLOCATABLE :: XLONG(:)
      REAL,ALLOCATABLE :: XKAPPA1(:), DELH_ALT(:)
      REAL,ALLOCATABLE :: FALT(:,:,:)
      INTEGER,ALLOCATABLE,DIMENSION(:) :: ICOMBUF
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: IJALTPE, ITEMP
      INTEGER :: NOBSPSAT(NUMALT)
      REAL,ALLOCATABLE, DIMENSION(:) :: ZCOMBUF, TEMP
      REAL,ALLOCATABLE :: ZCOMBUFS(:,:)
      REAL,ALLOCATABLE :: ZCOMBUFR(:)
      REAL,ALLOCATABLE,DIMENSION(:) :: WHGTAGB, CWHTAGB, RANGAGB, AVGWGT
      REAL,ALLOCATABLE,DIMENSION(:,:) :: ALTDATAPE, ZTEMP
      INTEGER :: NFRPEALT
      REAL :: WHMOD(NIBLO)
      LOGICAL :: LSQRT, LDUM
      LOGICAL :: LLCOMMON(NPROC)
      LOGICAL,ALLOCATABLE, DIMENSION(:) :: LLWITHIN 
      LOGICAL :: LLALT(NUMALT)


!*     VARIABLE     TYPE      PURPOSE.
!      --------     ----      --------

!     *IREAD*       INTGER    PROCESSOR WHICH WILL READ THE INPUT FILES.
!     *CLFID        CHARACTER THREE CHARACTER FILE IDENTIFIER.

      DATA CLFID, IREAD /'RFL', 1/

! ----------------------------------------------------------------------

!*    1. OPEN FILE AND READ.
!        -------------------

      IDUM=0
      NOBS=0
      NOBSPSAT=0

      IF((LMESSPASS .AND. IRANK .EQ. IREAD) .OR. .NOT.LMESSPASS) THEN
        IFIELD = 0
        IFAIL = 1

!*       FIND SPARE UNIT TO ATTACH INPUT FILE.

        UNIT: DO IUME=33,99,1
          INQUIRE ( UNIT=IUME, OPENED=LLOPENED)
          IF ( .NOT. LLOPENED) THEN
            IF(ITEST.GT.0) THEN
              WRITE(IU06,'("  GRFIELD : TAKES UNIT ", I2,
     &         " AS ALTIMETER DATA INPUT UNIT")') IUME
            ENDIF
            EXIT UNIT
          ENDIF
        ENDDO UNIT

!*       IF FOUND: SEEK INPUT FILE, READ HEADER AND FIELDS.
!*                 IFIELD MAY BE MISSING INDICATING THERE ARE NO FIELDS.

        IF (IUME.LT.100) THEN
          LLTEST=.FALSE.
          IF(ITEST.GT.0) LLTEST=.TRUE.
          CALL CONFILE (IU06, IUME, CDTPRO, CLFID, IFAIL, LLTEST)
          IF(IFAIL.EQ.0) THEN
            READ (IUME,ERR=9000,IOSTAT=IOS) CTIME, NOBS
            IF(NOBS.GT.0) THEN
!             the size of the 2nd dimension of IJALT is increased by 1
!             to keep track of the input order (see below).
              ALLOCATE(IJALT(NOBS,NIJALT+1))
              ALLOCATE(ALTDATA(NOBS,NALTDT))
              READ (IUME,ERR=9000,IOSTAT=IOS) 
     &                             ((IJALT(IOBS,IIJALT),IOBS = 1,NOBS),
     &                                                IIJALT = 1,NIJALT)
              READ (IUME,ERR=9000,IOSTAT=IOS) ALTDATA

!             REORGANISE THE IJ INDEX WHEN 2-D DECOMPOSITION
              IF(.NOT.LL1D .AND. NPROC.GT.1) THEN
                DO IOBS = 1, NOBS
                  IJALT(IOBS,1) = IJ2NEWIJ(IJALT(IOBS,1))
                ENDDO
              ENDIF

              DO IOBS = 1, NOBS
                IJALT(IOBS,NIJALT+1) = IOBS
              ENDDO

              DO ISAT=1,NUMALT
                DO IOBS = 1, NOBS
!             take only points that belong to the satellite
                  IF(IJALT(IOBS,2).EQ.IBUFRSAT(ISAT))
     &              NOBSPSAT(ISAT)=NOBSPSAT(ISAT)+1
                ENDDO
              ENDDO
!             check whether sum(NOBSPSAT)=NOBS
              NOBSDUM=0
              DO ISAT=1,NUMALT
                NOBSDUM=NOBSDUM+NOBSPSAT(ISAT)
              ENDDO
              IF(NOBSDUM.NE.NOBS) THEN
                WRITE(IU06,*) ' '
                WRITE(IU06,'("IN GRFIELD:  W A R N I N G  !!!!!!")')
                WRITE(IU06,*) 'THE TOTAL NUMBER OF ALT DATA DOES NOT '
                WRITE(IU06,*) 'CORRESPOND TO THE SUM OF THE DIFFERENT '
                WRITE(IU06,*) 'CONTRIBUTION '
                WRITE(IU06,*) 'NOBS = ',NOBS
                DO ISAT=1,NUMALT
                  WRITE(IU06,*) 'NOBSPSAT = ',NOBSPSAT(ISAT),ISAT
                ENDDO
                WRITE(IU06,*) ' '
                WRITE(IU06,*) 'NO DATA ASSIMILATION WILL TAKE PLACE !!!'
                WRITE(IU06,*) ' '
                WRITE(IU06,'("IN GRFIELD:  W A R N I N G  !!!!!!")')
                WRITE(IU06,*) ' '
                DO ISAT=1,NUMALT
                  NOBSPSAT(ISAT)=0
                ENDDO
                CALL FLUSH(IU06)
              ENDIF
            ENDIF
            CLOSE (IUME)
          ELSE
            WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!!!")')
            WRITE(IU06,'("NO FILES FOUND"                            )')
            WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!!!")')
            CALL FLUSH(IU06)
            NOBS=0
          ENDIF
        ELSE
!*       IF NOT FOUND: ISSUE WARNING MESSAGE.
          WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!!!!!")')
          WRITE(IU06,'("NO SPARE UNIT TO READ SATELLITE DATA FOUND"  )')
          WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!!!!!")')
          CALL FLUSH(IU06)
          NOBS=0
        ENDIF
        GOTO 9001
9000    WRITE(IU06,'("READING ERROR IN GRFIELD: ")')
        WRITE(IU06,'(25X,"I/O STATUS IS",I8)') IOS
        WRITE(IU06,'(25X," ABORT !!!!!!!!")')
        CALL ABORT1  
9001    CONTINUE
      ENDIF

!     SEND NOBS TO THE OTHER PE'S TO LET THEM KNOW WHETHER THERE ARE
!     DATA OR NOT.
      IF(LMESSPASS) THEN
        CALL MPBCASTINTFLD(IREAD,KTAG,NUMALT,1,NOBSPSAT)
        KTAG=KTAG+1
      ENDIF
      NOBS=0
      DO ISAT=1,NUMALT
        NOBS=NOBS+NOBSPSAT(ISAT)
      ENDDO

!     IF NO DATA, NOTHING TO BE DONE
      IF (NOBS.LE.0) THEN
        IF(ALLOCATED(WHGTAG)) WHGTAG = ZMISS
        IF(ALLOCATED(CWHTAG)) CWHTAG = ZMISS
        IF(ALLOCATED(RANGAG)) RANGAG = ZMISS
        NOBSPE = 0
        RETURN
      ENDIF

!     add a small number to DISTMAX as a safety precaution
!     i.e each PE will get slightly more data that they should
!     strictly have.
      COSDISTMAX=COS(DISTMAX+0.1)

! ----------------------------------------------------------------------

!*    2. SEND FIELDS TO ALL OTHER PES.
!        -----------------------------

      IF(LMESSPASS .AND. NPROC.GT.1) THEN
!       GET THE NUMBER OF DATA PER PE.
        IF(IRANK .EQ. IREAD) THEN
          NOBSPE = 0
          DO IR = 1, NPROC
            DO IOBS = 1, NOBS
              IF( IJALT(IOBS,1) .GE. NSTART(IR) .AND.
     &            IJALT(IOBS,1) .LE. NEND(IR) ) THEN
                NOBSPE(IR) = NOBSPE(IR) + 1
              ENDIF
            ENDDO
          ENDDO
        ENDIF
        CALL MPBCASTINTFLD(IREAD,KTAG,NPROC,1,NOBSPE)
        KTAG=KTAG+1


!       DETERMINE THE NUMBER PF PE'S EACH PE WILL NEED TO GET DATA
!       FROM
        IR=IRANK
        NFRPEALT = 0 
        NOBSPEMAX = 0
        DO JR = 1, IR-1
          IF(INTLMAX(IR,JR).EQ.1 .AND. NOBSPE(JR).GT.0) THEN
            NFRPEALT = NFRPEALT + 1
            NOBSPEMAX=MAX(NOBSPEMAX,NOBSPE(JR))
          ENDIF
        ENDDO
        DO JR = IR+1,NPROC
          IF(INTLMAX(IR,JR).EQ.1 .AND. NOBSPE(JR).GT.0) THEN
            NFRPEALT = NFRPEALT + 1
            NOBSPEMAX=MAX(NOBSPEMAX,NOBSPE(JR))
          ENDIF
        ENDDO


!       DISTRIBUTE THE DATA FOR EACH PE.
        IF(IRANK .EQ. IREAD) THEN
!         EXTRACT THE DATA THAT ARE ON EACH PE.
          DO IR = 1, NPROC
            IF(NOBSPE(IR) .GT. 0 .AND. IR.NE.IREAD ) THEN

              IZCOMLEN = (NIJALT+1+NALTDT) * NOBSPE(IR)
              ALLOCATE(ZCOMBUF(IZCOMLEN))
              ICOUNT = 0
              DO IIJALT = 1, NIJALT+1
                DO IOBS = 1, NOBS
                  IF( IJALT(IOBS,1) .GE. NSTART(IR) .AND.
     &                IJALT(IOBS,1) .LE. NEND(IR) ) THEN
                    ICOUNT = ICOUNT + 1
                    ZCOMBUF(ICOUNT) = FLOAT(IJALT(IOBS,IIJALT))
                  ENDIF
                ENDDO
              ENDDO
              DO IALTDT = 1, NALTDT
                DO IOBS = 1, NOBS
                  IF( IJALT(IOBS,1) .GE. NSTART(IR) .AND.
     &                IJALT(IOBS,1) .LE. NEND(IR) ) THEN
                    ICOUNT = ICOUNT + 1
                    ZCOMBUF(ICOUNT) = ALTDATA(IOBS,IALTDT)
                  ENDIF
                ENDDO
              ENDDO
              CALL GSTATS(618,0)
              CALL MPL_SEND(ZCOMBUF(1:ICOUNT),KDEST=IR,KTAG=KTAG,
     &                      CDSTRING='GRFIELD 1:')
              CALL GSTATS(618,1)
              DEALLOCATE(ZCOMBUF)

            ENDIF
          ENDDO

!         only keep the contribution that are on PE IREAD
          IR=IREAD
          IF(NOBSPE(IR).GT.0) THEN
            IZCOMLEN = (NIJALT+1+NALTDT) * NOBSPE(IR)
            ALLOCATE(ZCOMBUF(IZCOMLEN))
            ICOUNT = 0
            DO IIJALT = 1, NIJALT+1
              DO IOBS = 1, NOBS
                IF( IJALT(IOBS,1) .GE. NSTART(IR) .AND.
     &              IJALT(IOBS,1) .LE. NEND(IR) ) THEN
                  ICOUNT = ICOUNT + 1
                  ZCOMBUF(ICOUNT) = FLOAT(IJALT(IOBS,IIJALT))
                ENDIF
              ENDDO
            ENDDO
            DO IALTDT = 1, NALTDT
              DO IOBS = 1, NOBS
                IF( IJALT(IOBS,1) .GE. NSTART(IR) .AND.
     &            IJALT(IOBS,1) .LE. NEND(IR) ) THEN
                  ICOUNT = ICOUNT + 1
                  ZCOMBUF(ICOUNT) = ALTDATA(IOBS,IALTDT)
                ENDIF
              ENDDO
            ENDDO
          ENDIF

          DEALLOCATE(IJALT)
          DEALLOCATE(ALTDATA)

        ELSE
          CALL GSTATS(618,0)
          IF(NOBSPE(IRANK) .GT. 0 ) THEN
            IZCOMLEN = (NIJALT+1+NALTDT) * NOBSPE(IRANK)
            ALLOCATE(ZCOMBUF(IZCOMLEN))
            CALL MPL_RECV(ZCOMBUF(1:IZCOMLEN),KSOURCE=IREAD,KTAG=KTAG,
     &         KOUNT=KRCOUNT,KRECVTAG=KRTAG,CDSTRING='GRFIELD 1:')
            IF(KRCOUNT.NE.IZCOMLEN) CALL MPL_ABORT
     &        ('MPL_RECV ERROR in GRFIELD 1 :MISMATCHED MSG LENGTH')
            IF(KRTAG.NE.KTAG) CALL MPL_ABORT
     &        ('MPL_RECV ERROR in GRFIELD 1 : MISMATCHED TAGS' )
          ENDIF
          CALL GSTATS(618,1)
        ENDIF
        KTAG = KTAG + 1

!       KEEP THE DATA ON EACH PE.
        IF(NOBSPE(IRANK).GT.0) THEN
          ALLOCATE(IJALT(NOBSPE(IRANK),NIJALT+1))
          ALLOCATE(ALTDATA(NOBSPE(IRANK),NALTDT))
          ICOUNT = 0
          DO IIJALT = 1, NIJALT+1
            DO IOBS = 1, NOBSPE(IRANK)
              ICOUNT = ICOUNT + 1
              IJALT(IOBS,IIJALT) = NINT(ZCOMBUF(ICOUNT)) 
            ENDDO
          ENDDO
          DO IALTDT = 1, NALTDT
            DO IOBS = 1, NOBSPE(IRANK)
              ICOUNT = ICOUNT + 1
              ALTDATA(IOBS,IALTDT) = ZCOMBUF(ICOUNT) 
            ENDDO
          ENDDO
          DEALLOCATE(ZCOMBUF)
        ENDIF

      ELSE
        NOBSPE = NOBS
      ENDIF

!     FIND OUT WHICH SATELLITE IS PRESENT
      NALTAVLB=0
      DO ISAT=1,NUMALT
        IF(NOBSPSAT(ISAT).GT.0) THEN
          LLALT(ISAT)=.TRUE.
          NALTAVLB=NALTAVLB+1
        ELSE
          LLALT(ISAT)=.FALSE.
        ENDIF
        WRITE(IU06,'(2X,I2,A,I3,A,I7,A,L7)')
     &        ISAT,'. SATELLITE:',IBUFRSAT(ISAT),
     &     '  ... NOBSPSAT=',NOBSPSAT(ISAT),
     &     '  ... LLALT=',LLALT(ISAT)
      ENDDO

! ----------------------------------------------------------------------

!*    3. APPLY CORRECTION TO WAVE HEIGHT DATA ACCORDING TO 
!        SUPPLIED LINEAR REGRESSION COEFFICIENTS
!        --------------------------------------------------
      DO ISAT=1,NUMALT
        IF(LALTLRGR(ISAT) .AND. LLALT(ISAT)) THEN
          DO IOBS = 1, NOBSPE(IRANK)
            IF (ALTDATA(IOBS,1).GT.0. .AND. 
     &          IJALT(IOBS,2).EQ.IBUFRSAT(ISAT) ) THEN
              ALTDATA(IOBS,1) = 
     &        HSCOEFCOR(ISAT)*ALTDATA(IOBS,1)+HSCONSCOR(ISAT)
            ENDIF
          ENDDO
          IF(ITEST.GT.0) THEN
            WRITE(IU06,*)' GRFIELD : ALTIMETER ',IBUFRSAT(ISAT) 
            WRITE(IU06,*)'           ALTIMETER WAVE HEIGHTS CORRECTED '
            WRITE(IU06,*)'           Hsnew= ',HSCOEFCOR(ISAT),' Hs + ',
     &                                        HSCONSCOR(ISAT)
            WRITE(IU06,*)' '
          ENDIF
        ENDIF
      ENDDO

! ----------------------------------------------------------------------

!     4. APPLY SCREENING OF DATA BASED ON MODEL FIELD
!        --------------------------------------------

      DO ISAT=1,NUMALT

        WRITE(IU06,*)' '
        WRITE(IU06,*)'   ALTIMETER ',IBUFRSAT(ISAT)

        NREJICE=0
        NREJSD=0
        NREJBG=0
        NSAT=0

        IF(ALTSDTHRSH(ISAT).LT.0.) THEN
          IF(XDELLA.LT.0.3) THEN
            ALTSDTHRSH(ISAT) = 0.75
          ELSE
            ALTSDTHRSH(ISAT) = 0.6
          ENDIF
          WRITE(IU06,'(A,F6.2,A)')
     &        '    SUSPICIOUS DATA THRESHOLD NOT GIVEN. IT IS SET TO ',
     &         ALTSDTHRSH(ISAT)*100.,'%.'
        ENDIF

        IF(LLALT(ISAT)) THEN

          DO IOBS = 1, NOBSPE(IRANK)
!           take only points that belong to the satellite
            IF(IJALT(IOBS,2).EQ.IBUFRSAT(ISAT)) THEN

!             EXCLUDE POINTS WHICH ARE ON MODEL SEA ICE

              IJ=IJALT(IOBS,1)
              IF (ICEMASK(IJ,1).EQ.0) THEN
                ALTDATA(IOBS,1)=-3.
                NREJICE=NREJICE+1
              ENDIF

!              LOOK FOR SUSPICIOUS DATA
              IF (ALTDATA(IOBS,1).GT.0. .AND. WHMOD(IJ).GT.0. ) THEN
                THRSHLD=ABS(ALTDATA(IOBS,1)-WHMOD(IJ))/WHMOD(IJ)
                IF(THRSHLD.GE.ALTSDTHRSH(ISAT)) THEN
                  NREJSD=NREJSD+1
                ENDIF
              ENDIF

!              BACKGROUND CHECK (remove data that deviate too much
!                                from the first guess)

              IF (ALTDATA(IOBS,1).GT.0. .AND. WHMOD(IJ).GT.0. ) THEN
                THRSHLD=(ALTDATA(IOBS,2)/SIGMOD)*
     &                  (ALTDATA(IOBS,1)-WHMOD(IJ))**2/
     &                  (ALTDATA(IOBS,1)*WHMOD(IJ))
                HSCUT=MIN(HSALTCUT(ISAT),ALTDATA(IOBS,2))
                IF(THRSHLD.GE.ALTBGTHRSH(ISAT)) THEN
                  ALTDATA(IOBS,1)=-3.
                  NREJBG=NREJBG+1
                ELSE IF(ALTDATA(IOBS,1).LE.HSCUT .OR.
     &                  WHMOD(IJ).LE.HSCUT ) THEN
                  ALTDATA(IOBS,1)=-3.
                  NREJBG=NREJBG+1
                ENDIF
              ENDIF
            ENDIF
          ENDDO
          NSAT=1
        ENDIF

        NTOTREJICE=NREJICE
        NTOTREJSD=NREJSD
        NTOTREJBG=NREJBG
        NTOTSAT=NSAT

        IF(NPROC.GT.1) THEN
!         SEND NREJICE, NREJSD AND  NREJBG TO ALL OTHER PE's
          MPLENGTH=4
          ALLOCATE(ICOMBUF(MPLENGTH))

          DO IP=1,NPROC
            IF(IRANK.EQ.IP) THEN
              ICOMBUF(1)=NREJICE
              ICOMBUF(2)=NREJSD
              ICOMBUF(3)=NREJBG
              IF(LLALT(ISAT)) THEN
                NSAT=1
              ELSE
                NSAT=0
              ENDIF
              ICOMBUF(4)=NSAT
            ENDIF

            CALL GSTATS(618,0)
            CALL MPL_BROADCAST(ICOMBUF(1:MPLENGTH),KROOT=IP,
     &        KTAG=KTAG,CDSTRING='GRFIELD 2:')
            CALL GSTATS(618,1)

            IF(IRANK.NE.IP) THEN
              NTOTREJICE=NTOTREJICE+ICOMBUF(1)
              NTOTREJSD=NTOTREJSD+ICOMBUF(2)
              NTOTREJBG=NTOTREJBG+ICOMBUF(3)
              NTOTSAT=NTOTSAT+ICOMBUF(4)
            ENDIF

            KTAG = KTAG + 1
          ENDDO

          DEALLOCATE(ICOMBUF)
        ENDIF

        WRITE(IU06,*)'   NUMBER OF GRIDDED WAVE HEIGHTS= ',
     &                   NOBSPSAT(ISAT)
        IF(NOBSPSAT(ISAT).GT.0) THEN
          WRITE(IU06,*)'   NUMBER OF SUBAREA WITH DATA= ',NTOTSAT 
          WRITE(IU06,*)'   NUMBER OVER MODEL SEA ICE= ',NTOTREJICE 
          WRITE(IU06,*)'   NUMBER OF SUSPICIOUS DATA= ',NTOTREJSD
          WRITE(IU06,*)'   NUMBER REJECTED BY THE BACKGROUND CHECK= ',
     &                   NTOTREJBG
        ENDIF

        IF(NTOTSAT.GT.0) THEN
 
          IF(XDELLA.LT.0.3) THEN
            REJECT=0.75
          ELSE
            REJECT=0.5
          ENDIF

          REJRATIO=FLOAT(NTOTREJSD)/FLOAT(NOBS)

          IF(REJRATIO.GT.REJECT) THEN
            WRITE(IU06,*)'*********************************************'
            WRITE(IU06,*)' '
            WRITE(IU06,*)'   ALTIMETER ',IBUFRSAT(ISAT)
            WRITE(IU06,*)'  THE NUMBER OF SUSPICIOUS ALT DATA IS OVER ',
     &                   REJECT*100,' %'
            WRITE(IU06,*)' '
            WRITE(IU06,*)'  THEY ARE ALL REJECTED AS A PRECAUTION !' 
            WRITE(IU06,*)' '
            WRITE(IU06,*)'*********************************************'
            DO IOBS = 1, NOBSPE(IRANK)
              IF(IJALT(IOBS,2).EQ.IBUFRSAT(ISAT)) ALTDATA(IOBS,1)=-3.
            ENDDO
            IF(IRANK.EQ.1) THEN
              IUMAIL=I_GET_UNIT(IU06,'MM','S','F',0)
              WRITE(IUMAIL,*)'          !!!!!! WARNING !!!!!! '
              WRITE(IUMAIL,*)' NUMBER OF SUSPICIOUS DATA IS OVER',
     &                       REJECT*100,' % ',REJRATIO*100,' %'
              WRITE(IUMAIL,*)' IN EXPERIMENT ',YEXPVER 
              CALL WSTREAM_STRG(KSEK1(40),CSTRM,IDUM,IDUM,CDUM,IDUM,
     &                          LDUM)
              WRITE(IUMAIL,*)' STREAM        ',CSTRM
              WRITE(IUMAIL,*)' DOMAIN        ',CLDOMAIN
              WRITE(IUMAIL,*)' AT TIME       ',CDTPRO
              WRITE(IUMAIL,*)'          !!!!!! WARNING !!!!!! '
              CLOSE(IUMAIL)
              CALL SYSTEM (' mail wab < MM ' )
            ENDIF
          ENDIF

        ENDIF

      ENDDO ! END LOOP ON SATELLITE

! ----------------------------------------------------------------------

!*    5. APPLY CORRECTION TO WAVE HEIGHT DATA ACCORDING TO SROKOSZ.
!        ---------------------------------------------------------
!        SEE ECMWF TECH MEMO 269 or JANSSEN 2000: ECMWF wave modelling
!        and satellite altimeter wave data, Satellites, Oceanography and
!        Society, ed. D Halpern, Elsevier: Amsterdam. pp 35-56.

!     5.1 COLLOCATE WAMSPECTRA WITH ALT-HS OBSERVATIONS.
!         ---------------------------------------------

      IF(ALLOCATED(WHGTAG)) WHGTAG = ZMISS
      IF(ALLOCATED(CWHTAG)) CWHTAG = ZMISS
      IF(ALLOCATED(RANGAG)) RANGAG = ZMISS
      ALLOCATE(WHGTAGB(NSTART(IRANK):NEND(IRANK)))
      ALLOCATE(CWHTAGB(NSTART(IRANK):NEND(IRANK)))
      ALLOCATE(RANGAGB(NSTART(IRANK):NEND(IRANK)))
      ALLOCATE(AVGWGT(NSTART(IRANK):NEND(IRANK)))
      WHGTAGB=0.
      CWHTAGB=0.
      RANGAGB=0.
      AVGWGT=0.
 
      IF(NOBSPE(IRANK).GT.0) THEN
        NERS=NOBSPE(IRANK)
        ALLOCATE(XKAPPA1(NERS))
        ALLOCATE(DELH_ALT(NERS))
        ALLOCATE(FALT(NERS,NANG,NFRE))
         
        DO M=1,NFRE
          DO K=1,NANG
            DO IOBS = 1, NOBSPE(IRANK)
              IF (ALTDATA(IOBS,1).GT.0.) THEN
                IJ=IJALT(IOBS,1)
                FALT(IOBS,K,M)=FL1(IJ,K,M)
              ELSE
                FALT(IOBS,K,M)=FLMIN
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      
        CALL SKEWNESS(IU06,FALT,NERS,XKAPPA1,DELH_ALT)

!       5.3 APPLY SROKOSZ CORRECTION.
!           -------------------------
     
        DO ISAT=1,NUMALT
          IF(LLALT(ISAT)) THEN

            IF(LALTCOR(ISAT).AND.ITEST.GT.0) THEN
              WRITE(IU06,*) ' GRFIELD : ALTIMETER ',IBUFRSAT(ISAT)
              WRITE(IU06,*) '           WAVE HEIGHTS CORRECTED '
              WRITE(IU06,*) ' '
            ENDIF

            DO IOBS = 1, NOBSPE(IRANK)
              IF (ALTDATA(IOBS,1).GT.0. .AND.
     &            IJALT(IOBS,2).EQ.IBUFRSAT(ISAT)) THEN
                RCOR = DELH_ALT(IOBS)
                HS = ALTDATA(IOBS,1) 
                IF (LALTCOR(ISAT)) THEN
                  XN = HS**2/16. + XKAPPA2(ISAT)
                  HS_NEW = 16.*(XKAPPA1(IOBS)*XN - XKAPPA2(ISAT))
                  IF ( HS_NEW >  0.0 ) THEN
                    HS_NEW = SQRT(HS_NEW)
                    ALTDATA(IOBS,1) = HS_NEW
                  ELSE
                    HS_NEW = 0.0
                  ENDIF
                ELSE
                  HS_NEW = HS
                ENDIF

!               SAVE ALTIMETER DATA TO GRIB
                IF (LALTGRDOUT(ISAT)) THEN
                  IJ=IJALT(IOBS,1)
                  IF (IJ.GE.NSTART(IRANK) .AND. IJ.LE.NEND(IRANK)) THEN
!                   WHGTAGB CONTAINS THE ORIGINAL DATA
!                   A WEIGHTED AVERAGE IS COMPUTED IF MORE THAN ONE SATELLITE
                    IF(HS_NEW.GT.0.) THEN
                      WGT=1./ALTDATA(IOBS,2)**2
                      WHGTAGB(IJ) = WHGTAGB(IJ) + WGT*ALTDATA(IOBS,3)
                      CWHTAGB(IJ) = CWHTAGB(IJ) + WGT*HS_NEW
                      RANGAGB(IJ) = RANGAGB(IJ) + WGT*HS_NEW*RCOR 
                      AVGWGT(IJ) = AVGWGT(IJ) + WGT
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO

        DO IJ = NSTART(IRANK), NEND(IRANK)
          IF(WHGTAGB(IJ).GT.0.) THEN
            WHGTAGB(IJ) = WHGTAGB(IJ)/AVGWGT(IJ)
            CWHTAGB(IJ) = CWHTAGB(IJ)/AVGWGT(IJ)
            RANGAGB(IJ) = RANGAGB(IJ)/AVGWGT(IJ)
            RANGAGB(IJ) = RANGAGB(IJ)/CWHTAGB(IJ)
          ELSE
            WHGTAGB(IJ) = ZMISS 
            CWHTAGB(IJ) = ZMISS 
            RANGAGB(IJ) = ZMISS 
            RANGAGB(IJ) = ZMISS 
          ENDIF
        ENDDO

        IF( ALLOCATED(XKAPPA1) ) DEALLOCATE(XKAPPA1)
        IF( ALLOCATED(DELH_ALT) ) DEALLOCATE(DELH_ALT)
        IF( ALLOCATED(FALT) ) DEALLOCATE(FALT)
      ENDIF ! ON NOBSPE


!     COLLECT THE OUTPUT OF THE ALTIMETER CORRECTIONS ONTO THE OUTPUT PE
      ALLOCATE(TEMP(NIBLO))
      DO IJ = NSTART(IRANK), NEND(IRANK)
        TEMP(IJ) = WHGTAGB(IJ)
      ENDDO
      CALL MPGATHERSCFLD(IPFGTBL(22),KTAG,NSTART,NEND,TEMP)
      KTAG=KTAG+1
      IF(ALLOCATED(WHGTAG).AND.IPFGTBL(22).EQ.IRANK) THEN
        LSQRT=.FALSE.
        CALL MAKEGRID (TEMP , WHGTAG, 1, NEND(NPROC), 1, LSQRT, ZMISS)
      ENDIF

      DO IJ = NSTART(IRANK), NEND(IRANK)
        TEMP(IJ) = CWHTAGB(IJ)
      ENDDO
      CALL MPGATHERSCFLD(IPFGTBL(23),KTAG,NSTART,NEND,TEMP)
      KTAG=KTAG+1
      IF(ALLOCATED(CWHTAG).AND.IPFGTBL(23).EQ.IRANK) THEN
        LSQRT=.FALSE.
        CALL MAKEGRID (TEMP , CWHTAG, 1, NEND(NPROC), 1, LSQRT, ZMISS)
      ENDIF

      DO IJ = NSTART(IRANK), NEND(IRANK)
        TEMP(IJ) = RANGAGB(IJ)
      ENDDO
      CALL MPGATHERSCFLD(IPFGTBL(24),KTAG,NSTART,NEND,TEMP)
      KTAG=KTAG+1
      IF(ALLOCATED(RANGAG).AND.IPFGTBL(24).EQ.IRANK) THEN
        LSQRT=.FALSE.
        CALL MAKEGRID (TEMP , RANGAG, 1, NEND(NPROC), 1, LSQRT, ZMISS)
      ENDIF
      DEALLOCATE(TEMP)

      DEALLOCATE(WHGTAGB)
      DEALLOCATE(CWHTAGB)
      DEALLOCATE(RANGAGB)
      DEALLOCATE(AVGWGT)

      
!     GET ALL RELEVANT DATA ON EACH PE.
!     BY DISTRIBUTING THEM ONTO THE DIFFERENT PE'S

      IF(LMESSPASS .AND. NPROC.GT. 1) THEN
        IZCOMLENMAX=(NIJALT+1+NALTDT) * NOBSPE(IRANK)
        NDISPE=0
        DO IR = 1, NPROC
          IF(INTLMAX(IR,IRANK).EQ.1 .AND. IR.NE.IRANK) THEN
            IF(NOBSPE(IRANK) .GT. 0 ) THEN
              LLCOMMON(IR)=.TRUE.
              NDISPE=NDISPE+1
            ELSE
              LLCOMMON(IR)=.FALSE.
            ENDIF
          ELSE
            LLCOMMON(IR)=.FALSE.
          ENDIF
        ENDDO

        IF(IZCOMLENMAX.GT.0) THEN
!         PACKING THE CONTRIBUTION WHICH MIGHT BE COMMON
          ALLOCATE(ZCOMBUFS(IZCOMLENMAX,NDISPE))
          ALLOCATE(ISENDREQ(NDISPE))
          IDISPE=0
          DO IR = 1, NPROC
            IF(LLCOMMON(IR)) THEN
              IDISPE=IDISPE+1
              ICOUNT = 0
              DO IIJALT = 1, NIJALT+1
                DO IOBS = 1, NOBSPE(IRANK)
                  ICOUNT = ICOUNT + 1
                  ZCOMBUFS(ICOUNT,IDISPE) = FLOAT(IJALT(IOBS,IIJALT))
                ENDDO
              ENDDO
              DO IALTDT = 1, NALTDT
                DO IOBS = 1, NOBSPE(IRANK)
                  ICOUNT = ICOUNT + 1
                  ZCOMBUFS(ICOUNT,IDISPE) = ALTDATA(IOBS,IALTDT)
                ENDDO
              ENDDO
              ISENDLEN(IR)=ICOUNT
            ELSE
              ISENDLEN(IR)=0
            ENDIF
          ENDDO

!         SENDING THE CONTRIBUTION WHICH MIGHT BE COMMON
!         SEND NON BLOCKING THE BUFFERS

          IDISPE=0
          CALL GSTATS(618,0)
          DO IR = 1, NPROC
            IF(LLCOMMON(IR)) THEN
              IDISPE=IDISPE+1
              ICOUNT = ISENDLEN(IR) 
              CALL MPL_SEND(ZCOMBUFS(1:ICOUNT,IDISPE),
     &                      KDEST=IR,KTAG=KTAG,
     &                      KMP_TYPE=JP_NON_BLOCKING_STANDARD,
     &                      KREQUEST=ISENDREQ(IDISPE),
     &                      CDSTRING='GRFIELD 3:')
            ENDIF
          ENDDO
          CALL GSTATS(618,1)
        ENDIF

        ALLOCATE(XLONG(NSTART(IRANK):NEND(IRANK)))
        DO IJ=NSTART(IRANK),NEND(IRANK)
          XLONG(IJ) = REAL(IXLG(IJ,1)-1)*ZDELLO(KXLT(IJ,1))
        ENDDO

!       RECEIVING THE CONTRIBUTIONS WHICH MIGHT BE COMMON

        DO IPR = 1, NFRPEALT

          IZCOMLEN = (NIJALT+1+NALTDT) * NOBSPEMAX
          ALLOCATE(ZCOMBUFR(IZCOMLEN))

          CALL GSTATS(618,0)
          CALL MPL_RECV(ZCOMBUFR(1:IZCOMLEN),KFROM=IR,KTAG=KTAG,
     &           KOUNT=KRCOUNT,KRECVTAG=KRTAG,CDSTRING='GRFIELD 4:')
          CALL GSTATS(618,1)
          IF(KRTAG.NE.KTAG) CALL MPL_ABORT
     &    ('MPL_RECV ERROR in GRFIELD 3 : MISMATCHED TAGS' )


          ALLOCATE(IJALTPE(NOBSPE(IR),NIJALT+1))
          ICOUNT = 0
          DO IIJALT = 1, NIJALT+1
            DO IOBS = 1, NOBSPE(IR)
              ICOUNT = ICOUNT + 1
              IJALTPE(IOBS,IIJALT) = NINT(ZCOMBUFR(ICOUNT)) 
            ENDDO
          ENDDO

          ALLOCATE(ALTDATAPE(NOBSPE(IR),NALTDT))
          DO IALTDT = 1, NALTDT
            DO IOBS = 1, NOBSPE(IR)
              ICOUNT = ICOUNT + 1
              ALTDATAPE(IOBS,IALTDT) = ZCOMBUFR(ICOUNT) 
            ENDDO
          ENDDO

          DEALLOCATE(ZCOMBUFR)

!         how many observations are actually relevant to PE IRANK
          NOBSIN=0
          ALLOCATE(LLWITHIN(NOBSPE(IR)))
          DO IOBS = 1, NOBSPE(IR)
            IJJ=IJALTPE(IOBS,1)
            KJ = KXLT(IJJ,1)
            XLONGKJ = REAL(IXLG(IJJ,1)-1)*ZDELLO(KJ)
            LLWITHIN(IOBS) = .FALSE.
            IF(KJ .GE. KMINLMAX(IRANK) .AND.
     &         KJ .LE. KMAXLMAX(IRANK) ) THEN

               DO IJ=NSTART(IRANK),NEND(IRANK)
                   KI = KXLT(IJ,1)
                   DELLON = XLONGKJ-XLONG(IJ)
                   COSLON = COS(DELLON*RAD)
                   DIS2 = COSLON*COSPH(KJ)*COSPH(KI)+
     &                    SINPH(KJ)*SINPH(KI)
                   IF(DIS2.GE.COSDISTMAX) THEN
                     NOBSIN=NOBSIN+1
                     LLWITHIN(IOBS) = .TRUE.
                     EXIT
                   ENDIF
               ENDDO

            ENDIF
          ENDDO

          IF(NOBSIN.GT.0) THEN
            NOBSPEOLD = NOBSPE(IRANK)
            NOBSPE(IRANK) = NOBSPE(IRANK) + NOBSIN

            IF(NOBSPEOLD.GT.0) THEN
              ALLOCATE(ITEMP(NOBSPEOLD,NIJALT+1))
              ITEMP = IJALT
              DEALLOCATE(IJALT)
              ALLOCATE(IJALT(NOBSPE(IRANK),NIJALT+1))
              DO IIJALT = 1, NIJALT+1
                DO IOBS = 1, NOBSPEOLD
                  IJALT(IOBS,IIJALT) = ITEMP(IOBS,IIJALT)
                ENDDO
              ENDDO
              DEALLOCATE(ITEMP)

              ALLOCATE(ZTEMP(NOBSPEOLD,NALTDT))
              ZTEMP = ALTDATA
              DEALLOCATE(ALTDATA)
              ALLOCATE(ALTDATA(NOBSPE(IRANK),NALTDT))
              DO IALTDT = 1, NALTDT 
                DO IOBS = 1, NOBSPEOLD
                  ALTDATA(IOBS,IALTDT) = ZTEMP(IOBS,IALTDT) 
                ENDDO
              ENDDO
              DEALLOCATE(ZTEMP)
            ELSE
              ALLOCATE(IJALT(NOBSPE(IRANK),NIJALT+1))
              ALLOCATE(ALTDATA(NOBSPE(IRANK),NALTDT))
            ENDIF

            DO IIJALT = 1, NIJALT+1
              ICOUNT = 0
              DO IOBS = 1, NOBSPE(IR)
                IF(LLWITHIN(IOBS)) THEN
                  ICOUNT = ICOUNT + 1
                  IJALT(NOBSPEOLD+ICOUNT,IIJALT) =
     &                  IJALTPE(IOBS,IIJALT)
                ENDIF
              ENDDO
            ENDDO

            DO IALTDT = 1, NALTDT 
              ICOUNT = 0
              DO IOBS = 1, NOBSPE(IR)
                IF(LLWITHIN(IOBS)) THEN
                  ICOUNT = ICOUNT + 1
                  ALTDATA(NOBSPEOLD+ICOUNT,IALTDT) =
     &                    ALTDATAPE(IOBS,IALTDT)
                ENDIF
              ENDDO
            ENDDO
 
          ENDIF
          DEALLOCATE(LLWITHIN)
          DEALLOCATE(IJALTPE)
          DEALLOCATE(ALTDATAPE)

        ENDDO

        DEALLOCATE(XLONG)

! Wait any outstanding sends to complete

        CALL GSTATS(618,0)
        IF(IZCOMLENMAX.GT.0) THEN
          CALL MPL_WAIT(ZCOMBUFS(1:IZCOMLENMAX,1),KREQUEST=ISENDREQ,
     &                  CDSTRING='GRFIELD')
        ENDIF
        CALL GSTATS(618,1)
        IF(ALLOCATED(ISENDREQ)) DEALLOCATE(ISENDREQ)

        KTAG = KTAG + 1

        IF(ALLOCATED(ZCOMBUFS)) DEALLOCATE(ZCOMBUFS)

      ENDIF

!     REORDER THE DATA (to insure that the results are the same
!                       on any number of PE's)
      IF(NOBSPE(IRANK).GT.0) THEN
        ALLOCATE(ITEMP(NOBS,1))
        ITEMP=0
        DO IOBS = 1, NOBSPE(IRANK)
            ITEMP(IJALT(IOBS,NIJALT+1),1)=IJALT(IOBS,NIJALT+1)
        ENDDO
        ICOUNT=0
        DO IOBS = 1, NOBS
          IF(ITEMP(IOBS,1).GT.0) THEN
            ICOUNT=ICOUNT+1
            ITEMP(IOBS,1)=ICOUNT
          ENDIF
        ENDDO
        DO IOBS = 1, NOBSPE(IRANK)
          IJALT(IOBS,NIJALT+1)=ITEMP(IJALT(IOBS,NIJALT+1),1)
        ENDDO
        DEALLOCATE(ITEMP)

        ALLOCATE(ITEMP(NOBSPE(IRANK),NIJALT))
        DO IIJALT = 1, NIJALT
            DO IOBS = 1, NOBSPE(IRANK)
              ITEMP(IJALT(IOBS,NIJALT+1),IIJALT)=IJALT(IOBS,IIJALT)
            ENDDO
        ENDDO
        ALLOCATE(ZTEMP(NOBSPE(IRANK),NALTDT))
        DO IALTDT = 1, NALTDT 
          DO IOBS = 1, NOBSPE(IRANK)
            ZTEMP(IJALT(IOBS,NIJALT+1),IALTDT) = ALTDATA(IOBS,IALTDT)
          ENDDO
        ENDDO

        DEALLOCATE(IJALT)
        ALLOCATE(IJALT(NOBSPE(IRANK),NIJALT))
        DO IIJALT = 1, NIJALT
          DO IOBS = 1, NOBSPE(IRANK)
            IJALT(IOBS,IIJALT) = ITEMP(IOBS,IIJALT)
          ENDDO
        ENDDO
        DEALLOCATE(ITEMP)
        DO IALTDT = 1, NALTDT 
          DO IOBS = 1, NOBSPE(IRANK)
            ALTDATA(IOBS,IALTDT) = ZTEMP(IOBS,IALTDT)
          ENDDO
        ENDDO
        DEALLOCATE(ZTEMP)

      ENDIF

      RETURN
      END SUBROUTINE GRFIELD
