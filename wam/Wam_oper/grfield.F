      SUBROUTINE GRFIELD

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

!       *CALL* *GRFIELD*

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

      USE YOWALTAS , ONLY : NIJALT   ,NALTDT   ,IJALT    ,INTLMAX  ,
     &            KMINLMAX ,KMAXLMAX ,NOBSPE   ,ALTDATA
      USE YOWCOUT  , ONLY : IPFGTBL
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWICE   , ONLY : ICEMASK  ,FLMIN
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,AMOWEP   ,AMONOP   ,
     &            XDELLA   ,ZDELLO
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,KTAG
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,
     &            NIBLO    ,XKAPPA2  ,HSCOEFCOR,HSCONSCOR
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : CDTPRO   ,LALTCOR  ,LALTLRGR
      USE YOWSPEC  , ONLY : NSTART   ,NEND     ,FL1
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWINTP  , ONLY : WHGTTG   ,WHGTAG   ,CWHTAG   , RANGAG
      USE MPL_MODULE

! ----------------------------------------------------------------------

      LOGICAL LLOPENED, LLTEST
      INTEGER IREAD
      CHARACTER* 3 CLFID
      CHARACTER*12 CTIME
      REAL,ALLOCATABLE :: XKAPPA1(:), DELH_ALT(:)
      REAL,ALLOCATABLE :: FALT(:,:,:)
      INTEGER,ALLOCATABLE, DIMENSION(:) :: ICOMBUF, IKALT
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: IJALTPE, ITEMP
      REAL,ALLOCATABLE, DIMENSION(:) :: ZCOMBUF, TEMP
      REAL,ALLOCATABLE,DIMENSION(:) :: WHGTAGB, CWHTAGB, RANGAGB
      REAL,ALLOCATABLE,DIMENSION(:,:) :: ALTDATAPE, ZTEMP
      LOGICAL :: LSQRT


!*     VARIABLE     TYPE      PURPOSE.
!      --------     ----      --------

!     *IREAD*       INTGER    PROCESSOR WHICH WILL READ THE INPUT FILES.
!     *CLFID        CHARACTER THREE CHARACTER FILE IDENTIFIER.

      DATA CLFID, IREAD /'RFL', 1/

! ----------------------------------------------------------------------

!*    1. OPEN FILE AND READ.
!        -------------------

      NOBS=0

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
            READ (IUME,ERR=9000) CTIME, NOBS
            IF(NOBS.GT.0) THEN
!             the size of the second dimension of IJALT is increased by 1
!             to keep track of the input order (see below).
              ALLOCATE(IJALT(NOBS,NIJALT+1))
              ALLOCATE(ALTDATA(NOBS,NALTDT))
              READ (IUME,ERR=9000) ((IJALT(IOBS,IIJALT),IOBS = 1,NOBS),
     &                                                IIJALT = 1,NIJALT)
              READ (IUME,ERR=9000) ALTDATA
              DO IOBS = 1, NOBS
                IJALT(IOBS,NIJALT+1) = IOBS
              ENDDO
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
9000    WRITE(IU06,'("READING ERROR IN GRFIELD:  ABORT !!!!!!!!")')
        CALL ABORT1  
9001    CONTINUE
      ENDIF

!     SEND NOBS TO THE OTHER PE'S TO LET THEM KNOW WETHER THERE ARE
!     DATA OR NOT.
      IF(LMESSPASS) THEN
        CALL MPBCASTINTFLD(IREAD,KTAG,1,1,NOBS)
        KTAG=KTAG+1
      ENDIF

!     IF NO DATA, NOTHING TO BE DONE
      IF (NOBS.LE.0) THEN
        IF(ALLOCATED(WHGTAG)) WHGTAG = ZMISS
        IF(ALLOCATED(CWHTAG)) CWHTAG = ZMISS
        IF(ALLOCATED(RANGAG)) RANGAG = ZMISS
        NOBSPE = 0
        RETURN
      ENDIF

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

!       DISTRIBUTE THE DATA FOR EACH PE.
        IF(IRANK .EQ. IREAD) THEN
!         EXTRACT THE DATA THAT ARE ON EACH PE.
          DO IR = 1, NPROC
            IF(NOBSPE(IR) .GT. 0 .AND. IR.NE.IREAD ) THEN

              IICOMLEN = (NIJALT+1) * NOBSPE(IR)
              ALLOCATE(ICOMBUF(IICOMLEN))
              ICOUNT = 0
              DO IIJALT = 1, NIJALT+1
                DO IOBS = 1, NOBS
                  IF( IJALT(IOBS,1) .GE. NSTART(IR) .AND.
     &                IJALT(IOBS,1) .LE. NEND(IR) ) THEN
                    ICOUNT = ICOUNT + 1
                    ICOMBUF(ICOUNT) = IJALT(IOBS,IIJALT)
                  ENDIF
                ENDDO
              ENDDO
              CALL MPL_SEND(ICOMBUF(1:ICOUNT),KDEST=IR,KTAG=KTAG,
     &                      CDSTRING='GRFIELD 1:')
              DEALLOCATE(ICOMBUF)

              IZCOMLEN = NALTDT * NOBSPE(IR)
              ALLOCATE(ZCOMBUF(IZCOMLEN))
              ICOUNT = 0
              DO IALTDT = 1, NALTDT
                DO IOBS = 1, NOBS
                  IF( IJALT(IOBS,1) .GE. NSTART(IR) .AND.
     &                IJALT(IOBS,1) .LE. NEND(IR) ) THEN
                    ICOUNT = ICOUNT + 1
                    ZCOMBUF(ICOUNT) = ALTDATA(IOBS,IALTDT)
                  ENDIF
                ENDDO
              ENDDO
              CALL MPL_SEND(ZCOMBUF(1:ICOUNT),KDEST=IR,KTAG=KTAG+1,
     &                      CDSTRING='GRFIELD 2:')
              DEALLOCATE(ZCOMBUF)

            ENDIF
          ENDDO

!         only keep the contribution that are on PE IREAD
          IR=IREAD
          IF(NOBSPE(IR).GT.0) THEN
            IICOMLEN = (NIJALT+1) * NOBSPE(IR)
            ALLOCATE(ICOMBUF(IICOMLEN))
            ICOUNT = 0
            DO IIJALT = 1, NIJALT+1
              DO IOBS = 1, NOBS
                IF( IJALT(IOBS,1) .GE. NSTART(IR) .AND.
     &              IJALT(IOBS,1) .LE. NEND(IR) ) THEN
                  ICOUNT = ICOUNT + 1
                  ICOMBUF(ICOUNT) = IJALT(IOBS,IIJALT)
                ENDIF
              ENDDO
            ENDDO

            IZCOMLEN = NALTDT * NOBSPE(IR)
            ALLOCATE(ZCOMBUF(IZCOMLEN))
            ICOUNT = 0
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
          IF(NOBSPE(IRANK) .GT. 0 ) THEN
            IICOMLEN = (NIJALT+1) * NOBSPE(IRANK)
            ALLOCATE(ICOMBUF(IICOMLEN))
            CALL MPL_RECV(ICOMBUF(1:IICOMLEN),KSOURCE=IREAD,KTAG=KTAG,
     &         KOUNT=KRCOUNT,KRECVTAG=KRTAG,CDSTRING='GRFIELD 1:')
            IF(KRCOUNT.NE.IICOMLEN) CALL MPL_ABORT
     &        ('MPL_RECV ERROR in GRFIELD 1 :MISMATCHED MSG LENGTH')
            IF(KRTAG.NE.KTAG) CALL MPL_ABORT
     &        ('MPL_RECV ERROR in GRFIELD 1 : MISMATCHED TAGS' )


            IZCOMLEN = NALTDT * NOBSPE(IRANK)
            ALLOCATE(ZCOMBUF(IZCOMLEN))
            CALL MPL_RECV(ZCOMBUF(1:IZCOMLEN),KSOURCE=IREAD,KTAG=KTAG+1,
     &         KOUNT=KRCOUNT,KRECVTAG=KRTAG,CDSTRING='GRFIELD 2:')
            IF(KRCOUNT.NE.IZCOMLEN) CALL MPL_ABORT
     &        ('MPL_RECV ERROR in GRFIELD 2 :MISMATCHED MSG LENGTH')
            IF(KRTAG.NE.KTAG+1) CALL MPL_ABORT
     &        ('MPL_RECV ERROR in GRFIELD 2 : MISMATCHED TAGS' )

          ENDIF
        ENDIF
        KTAG = KTAG + 2

!       KEEP THE DATA ON EACH PE.
        IF(NOBSPE(IRANK).GT.0) THEN
          ALLOCATE(IJALT(NOBSPE(IRANK),NIJALT+1))
          ICOUNT = 0
          DO IIJALT = 1, NIJALT+1
            DO IOBS = 1, NOBSPE(IRANK)
              ICOUNT = ICOUNT + 1
              IJALT(IOBS,IIJALT) = ICOMBUF(ICOUNT) 
            ENDDO
          ENDDO
          DEALLOCATE(ICOMBUF)

          ALLOCATE(ALTDATA(NOBSPE(IRANK),NALTDT))
          ICOUNT = 0
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


! ----------------------------------------------------------------------

!     3. APPLY SCREENING OF DATA BASED ON MODEL FIELD
!        --------------------------------------------

!        EXCLUDE POINTS WHICH ARE ON MODEL ICE

      DO IOBS = 1, NOBSPE(IRANK)
        IJ=IJALT(IOBS,1)
        IF (ICEMASK(IJ,1).EQ.0) THEN
          ALTDATA(IOBS,1)=-3.
        ENDIF
      ENDDO

! ----------------------------------------------------------------------

!*    4. APPLY CORRECTION TO WAVE HEIGHT DATA ACCORDING TO 
!        SUPPLIED LINEAR REGRESSION COEFFICIENTS
!        --------------------------------------------------
!??????? multi-satellite issues
      IF(LALTLRGR) THEN
        DO IOBS = 1, NOBSPE(IRANK)
          IF (ALTDATA(IOBS,1).GT.0.)
     &        ALTDATA(IOBS,1) = HSCOEFCOR*ALTDATA(IOBS,1)+HSCONSCOR
        ENDDO
        IF(ITEST.GT.0) THEN
          WRITE(IU06,*)' GRFIELD : ALTIMETER WAVE HEIGHTS CORRECTED '
          WRITE(IU06,*)'           Hsnew= ',HSCOEFCOR,' Hs + ',HSCONSCOR
          WRITE(IU06,*)' '
        ENDIF
      ENDIF


!*    5. APPLY CORRECTION TO WAVE HEIGHT DATA ACCORDING TO SROKOSZ.
!        ---------------------------------------------------------

!     5.1 COLLOCATE WAMSPECTRA WITH ALT-HS OBSERVATIONS.
!         ---------------------------------------------

      IF(ALLOCATED(WHGTAG)) WHGTAG = ZMISS
      IF(ALLOCATED(CWHTAG)) CWHTAG = ZMISS
      IF(ALLOCATED(RANGAG)) RANGAG = ZMISS
      ALLOCATE(WHGTAGB(NSTART(IRANK):NEND(IRANK)))
      ALLOCATE(CWHTAGB(NSTART(IRANK):NEND(IRANK)))
      ALLOCATE(RANGAGB(NSTART(IRANK):NEND(IRANK)))
      WHGTAGB=ZMISS
      CWHTAGB=ZMISS
      RANGAGB=ZMISS
 
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
     
        IF(LALTCOR.AND.ITEST.GT.0) THEN
          WRITE(IU06,*) ' GRFIELD : ALTIMETER WAVE HEIGHTS CORRECTED '
          WRITE(IU06,*) ' '
        ENDIF

        DO IOBS = 1, NOBSPE(IRANK)
          IF (ALTDATA(IOBS,1).GT.0.) THEN
            RCOR = DELH_ALT(IOBS)

            HS = ALTDATA(IOBS,1) 
            XN = HS**2/16. + XKAPPA2
            HS_NEW = 16.*(XKAPPA1(IOBS)*XN - XKAPPA2)
            IF ( HS_NEW >  0.0 ) THEN
              HS_NEW = SQRT(HS_NEW)
            ELSE
              HS_NEW = 0.0
            ENDIF

!??????? multi-satellite issues
! we should only do it for one satellite !!!

            IJ=IJALT(IOBS,1)
            IF(HS.GT.0..AND.
     &         IJ.GE.NSTART(IRANK) .AND. IJ.LE.NEND(IRANK)) THEN
              WHGTAGB(IJ) = HS 
              CWHTAGB(IJ) = HS_NEW
              RANGAGB(IJ) = RCOR 
            ENDIF
!???????
            IF(LALTCOR .AND. HS_NEW >  0.0 ) ALTDATA(IOBS,1) = HS_NEW
          ENDIF
        ENDDO

      ENDIF

      IF( ALLOCATED(XKAPPA1) ) DEALLOCATE(XKAPPA1)
      IF( ALLOCATED(DELH_ALT) ) DEALLOCATE(DELH_ALT)
      IF( ALLOCATED(FALT) ) DEALLOCATE(FALT)

!     COLLECT THE OUTPUT OF THE ALTIMETER CORRECTIONS ONTO THE OUTPUT PE
      ALLOCATE(TEMP(NIBLO))
      DO IJ = NSTART(IRANK), NEND(IRANK)
        TEMP(IJ) = WHGTAGB(IJ)
      ENDDO
      CALL MPGATHERSCFLD(IPFGTBL(22),KTAG,NSTART,NEND,TEMP)
      KTAG=KTAG+1
      IF(ALLOCATED(WHGTAG).AND.IPFGTBL(22).EQ.IRANK) THEN
        LSQRT=.FALSE.
        CALL MAKEGRID (TEMP , WHGTAG, 1, NEND(NPROC), 1, LSQRT)
      ENDIF

      DO IJ = NSTART(IRANK), NEND(IRANK)
        TEMP(IJ) = CWHTAGB(IJ)
      ENDDO
      CALL MPGATHERSCFLD(IPFGTBL(23),KTAG,NSTART,NEND,TEMP)
      KTAG=KTAG+1
      IF(ALLOCATED(CWHTAG).AND.IPFGTBL(23).EQ.IRANK) THEN
        LSQRT=.FALSE.
        CALL MAKEGRID (TEMP , CWHTAG, 1, NEND(NPROC), 1, LSQRT)
      ENDIF

      DO IJ = NSTART(IRANK), NEND(IRANK)
        TEMP(IJ) = RANGAGB(IJ)
      ENDDO
      CALL MPGATHERSCFLD(IPFGTBL(24),KTAG,NSTART,NEND,TEMP)
      KTAG=KTAG+1
      IF(ALLOCATED(RANGAG).AND.IPFGTBL(24).EQ.IRANK) THEN
        LSQRT=.FALSE.
        CALL MAKEGRID (TEMP , RANGAG, 1, NEND(NPROC), 1, LSQRT)
      ENDIF
      DEALLOCATE(TEMP)

      DEALLOCATE(WHGTAGB)
      DEALLOCATE(CWHTAGB)
      DEALLOCATE(RANGAGB)

      
!     GET ALL RELEVANT DATA ON EACH PE.
!     BY DISTRIBUTING THEM ONTO THE DIFFERENT PE'S

      IF(LMESSPASS .AND. NPROC.GT. 1) THEN
!       SENDING THE CONTRIBUTION WHICH MIGHT BE COMMON
        DO IR = 1, NPROC
          IF(INTLMAX(IR,IRANK).EQ.1 .AND. IR.NE.IRANK) THEN
            IF(NOBSPE(IRANK) .GT. 0 ) THEN
              IICOMLEN = (NIJALT+1) * NOBSPE(IRANK)
              ALLOCATE(ICOMBUF(IICOMLEN))
              ICOUNT = 0
              DO IIJALT = 1, NIJALT+1
                DO IOBS = 1, NOBSPE(IRANK)
                  ICOUNT = ICOUNT + 1
                  ICOMBUF(ICOUNT) = IJALT(IOBS,IIJALT)
                ENDDO
              ENDDO
              CALL MPL_SEND(ICOMBUF(1:ICOUNT),KDEST=IR,KTAG=KTAG,
     &                      CDSTRING='GRFIELD 3:')

              DEALLOCATE(ICOMBUF)

              IZCOMLEN = NALTDT * NOBSPE(IRANK)
              ALLOCATE(ZCOMBUF(IZCOMLEN))
              ICOUNT = 0
              DO IALTDT = 1, NALTDT
                DO IOBS = 1, NOBSPE(IRANK)
                  ICOUNT = ICOUNT + 1
                  ZCOMBUF(ICOUNT) = ALTDATA(IOBS,IALTDT)
                ENDDO
              ENDDO
              CALL MPL_SEND(ZCOMBUF(1:ICOUNT),KDEST=IR,KTAG=KTAG+1,
     &                      CDSTRING='GRFIELD 4:')
              DEALLOCATE(ZCOMBUF)
            ENDIF
          ENDIF
        ENDDO

!       RECEIVING THE CONTRIBUTIONS WHICH MIGHT BE COMMON
        DO IR = 1, NPROC
          IF(INTLMAX(IR,IRANK).EQ.1 .AND. IR.NE.IRANK) THEN
            IF(NOBSPE(IR) .GT. 0 ) THEN

              IICOMLEN = (NIJALT+1) * NOBSPE(IR)
              ALLOCATE(ICOMBUF(IICOMLEN))

              CALL MPL_RECV(ICOMBUF(1:IICOMLEN),KSOURCE=IR,
     &           KTAG=KTAG, KOUNT=KRCOUNT,
     &           KRECVTAG=KRTAG,CDSTRING='GRFIELD 3:')
              IF(KRCOUNT.NE.IICOMLEN) CALL MPL_ABORT
     &        ('MPL_RECV ERROR in GRFIELD 3 :MISMATCHED MSG LENGTH')
              IF(KRTAG.NE.KTAG) CALL MPL_ABORT
     &        ('MPL_RECV ERROR in GRFIELD 3 : MISMATCHED TAGS' )

              ALLOCATE(IJALTPE(NOBSPE(IR),NIJALT+1))
              ICOUNT = 0
              DO IIJALT = 1, NIJALT+1
                DO IOBS = 1, NOBSPE(IR)
                  ICOUNT = ICOUNT + 1
                  IJALTPE(IOBS,IIJALT) = ICOMBUF(ICOUNT) 
                ENDDO
              ENDDO
              DEALLOCATE(ICOMBUF)

              IZCOMLEN = NALTDT * NOBSPE(IR)
              ALLOCATE(ZCOMBUF(IZCOMLEN))
              CALL MPL_RECV(ZCOMBUF(1:IZCOMLEN),KSOURCE=IR,KTAG=KTAG+1,
     &           KOUNT=KRCOUNT,KRECVTAG=KRTAG,CDSTRING='GRFIELD 4:')
              IF(KRCOUNT.NE.IZCOMLEN) CALL MPL_ABORT
     &        ('MPL_RECV ERROR in GRFIELD 4 :MISMATCHED MSG LENGTH')
              IF(KRTAG.NE.KTAG+1) CALL MPL_ABORT
     &        ('MPL_RECV ERROR in GRFIELD 4 : MISMATCHED TAGS' )

              ALLOCATE(ALTDATAPE(NOBSPE(IR),NALTDT))
              ICOUNT = 0
              DO IALTDT = 1, NALTDT
                DO IOBS = 1, NOBSPE(IR)
                  ICOUNT = ICOUNT + 1
                  ALTDATAPE(IOBS,IALTDT) = ZCOMBUF(ICOUNT) 
                ENDDO
              ENDDO
              DEALLOCATE(ZCOMBUF)

!             how many observations are actually relevant to PE IRANK
              NOBSIN=0
              ALLOCATE(IKALT(NOBSPE(IR)))
              DO IOBS = 1, NOBSPE(IR)
                IJ=IJALTPE(IOBS,1)
                IKALT(IOBS) = KXLT(IJ,1)
                IF(IKALT(IOBS) .GE. KMINLMAX(IRANK) .AND.
     &             IKALT(IOBS) .LE. KMAXLMAX(IRANK) ) THEN
                  NOBSIN=NOBSIN+1
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
                    IF(IKALT(IOBS) .GE. KMINLMAX(IRANK) .AND.
     &                 IKALT(IOBS) .LE. KMAXLMAX(IRANK) ) THEN
                      ICOUNT = ICOUNT + 1
                      IJALT(NOBSPEOLD+ICOUNT,IIJALT) =
     &                      IJALTPE(IOBS,IIJALT)
                    ENDIF
                  ENDDO
                ENDDO

                DO IALTDT = 1, NALTDT 
                  ICOUNT = 0
                  DO IOBS = 1, NOBSPE(IR)
                    IF(IKALT(IOBS) .GE. KMINLMAX(IRANK) .AND.
     &                 IKALT(IOBS) .LE. KMAXLMAX(IRANK) ) THEN
                      ICOUNT = ICOUNT + 1
                      ALTDATA(NOBSPEOLD+ICOUNT,IALTDT) =
     &                        ALTDATAPE(IOBS,IALTDT)
                    ENDIF
                  ENDDO
                ENDDO
 
              ENDIF
              DEALLOCATE(IKALT)
              DEALLOCATE(IJALTPE)
              DEALLOCATE(ALTDATAPE)

            ENDIF
          ENDIF
        ENDDO

        KTAG = KTAG + 2

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
