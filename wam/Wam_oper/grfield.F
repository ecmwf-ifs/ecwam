      SUBROUTINE GRFIELD (WHME, USME)

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

!       *CALL* *GRFIELD (WHME, USME, LDFAILD)*

!       *WHME*   REAL    OUTPUT FIELD OF MEASURED WAVE HEIGHTS.
!       *USME*   REAL    OUTPUT FIELD OF MEASURED SURFACE WINDS.


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

      USE YOWFRED  , ONLY : FR       ,TH       ,COSTH    ,SINTH
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWICE   , ONLY : ICEMASK
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,AMOWEP   ,AMONOP   ,
     &            XDELLA   ,ZDELLO
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,KTAG
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,
     &            XKAPPA2
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : CDTPRO   ,LALTCOR
      USE YOWSPEC  , ONLY : NSTART   ,NEND     ,FL1
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWINTP  , ONLY : WHGTTG   ,CDG      ,U10G     ,WHGTAG   ,
     &            CWHTAG   ,RANGAG

! ----------------------------------------------------------------------

      LOGICAL LLFAILED, LLOPENED
      INTEGER IREAD
      INTEGER,ALLOCATABLE :: IJALT(:,:), IJKAP(:)
      CHARACTER* 3 CLFID
      CHARACTER*12 CTIME
      REAL, DIMENSION(NGX,NGY) :: WHME, USME
      REAL,ALLOCATABLE :: XKAPPA1(:), XKAP(:), DELH_ALT(:), DELH(:)
      REAL,ALLOCATABLE :: FALT(:,:,:)


!*     VARIABLE     TYPE      PURPOSE.
!      --------     ----      --------

!     *WHME*        REAL      GRID WAVE HEIGHT.
!     *USME*        REAL      GRID FRICTION VELOCITY.
!     *IREAD*       INTGER    PROCESSOR WHICH WILL READ THE INPUT FILES.
!     *CLFID        CHARACTER THREE CHARACTER FILE IDENTIFIER.

      DATA CLFID, IREAD /'RFL', 1/

! ----------------------------------------------------------------------

!*    1. OPEN FILE AND READ.
!        -------------------

      IF(LMESSPASS .AND. IRANK .EQ. IREAD .OR. .NOT.LMESSPASS) THEN
        IFIELD = 0
        LLFAILED = .FALSE.
        IFAIL = 1

!*       FIND SPARE UNIT TO ATTACH INPUT FILE.

        UNIT: DO IUME=33,99,1
          INQUIRE ( UNIT=IUME, OPENED=LLOPENED)
          IF ( .NOT. LLOPENED) THEN
            WRITE(IU06,'("  GRFIELD : TAKES UNIT ", I2,
     &       " AS ALTIMETER DATA INPUT UNIT")') IUME
            EXIT UNIT
          ENDIF
        ENDDO UNIT

!*       IF FOUND: SEEK INPUT FILE, READ HEADER AND FIELDS.
!*                 IFIELD MAY BE MISSING INDICATING THERE ARE NO FIELDS.

        IF (IUME.LT.100) THEN
          CALL CONFILE (IU06, IUME, CDTPRO, CLFID, IFAIL)
          IF(IFAIL.EQ.0) THEN
            READ (IUME,ERR=9000) CTIME, IFIELD
            READ (IUME,ERR=9000) WHME, USME
            CLOSE (IUME)
          ELSE
            LLFAILED = .TRUE.
            WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!!!")')
            WRITE(IU06,'("NO FILES FOUND"                            )')
            WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!!!")')
            CALL FLUSH(IU06)
          ENDIF
        ELSE

!*       IF NOT FOUND: ISSUE WARNING MESSAGE.

          LLFAILED = .TRUE.
          WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!!!!!")')
          WRITE(IU06,'("NO SPARE UNIT TO READ SATELLITE DATA FOUND"  )')
          WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!!!!!")')
          CALL FLUSH(IU06)
        ENDIF

!*       IF ALL ATTEMPTS OF THE NEW METHOD FAILED TRY THE OLD CALL TO
!*       GRDATA.

 9000   IF ( LLFAILED .OR. IFIELD .EQ. 0 ) THEN
          WRITE(IU06,'("WARNING IN GRFIELD:    W A R N I N G  !!!!!!")')
          WRITE(IU06,'("NO GRIDDED FIELDS FOUND IN FILE"             )')
          WRITE(IU06,'("TRY THE OLD CALL TO GRDATA                  ")')
          WRITE(IU06,'("WARNING IN GRFIELD:    W A R N I N G  !!!!!!")')
          CALL FLUSH(IU06)
          IF (IFAIL .EQ. 0) CLOSE (IUME)
          CALL GRDATA (WHME, USME, NSTART, NEND) 
        ENDIF
      ENDIF

! ----------------------------------------------------------------------

!*    2. SEND FIELDS TO ALL OTHER PES.
!        -----------------------------

      IF(LMESSPASS) THEN
        CALL MPBCASTGRDFLD(IREAD,KTAG,NGX,NGY,WHME)
        KTAG=KTAG+1
        CALL MPBCASTGRDFLD(IREAD,KTAG,NGX,NGY,USME)
        KTAG=KTAG+1
      ENDIF

! ----------------------------------------------------------------------

!     3. APPLY SCREENING OF DATA BASED ON MODEL FIELD
!        --------------------------------------------

!        EXCLUDE POINTS WHICH ARE ON MODEL ICE

      IF(LMESSPASS) THEN
        DO IJ = NSTART(1), NEND(NPROC)
          IF (ICEMASK(IJ,1).EQ.0) THEN
            IX = IXLG(IJ,1)
            IY = NGY- KXLT(IJ,1) +1
            USME(IX,IY)=-3.
            WHME(IX,IY)=-3.
          ENDIF
        ENDDO
      ELSE
        DO IG = 1, IGL
          DO IJ = IJS(IG),IJL(IG)
            IF (ICEMASK(IJ,IG).EQ.0) THEN
              IX = IXLG(IJ,IG)
              IY = NGY- KXLT(IJ,IG) +1
              USME(IX,IY)=-3.
              WHME(IX,IY)=-3.
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!     SCREEN DATA WHICH ARE NOT OVER MODEL SEA POINTS
      DO J=1,NGY
        DO I=1,NGX
          IF (CDG(I,J) .LT. 0.000) THEN
            USME(I,J)=-3.
            WHME(I,J)=-3.
          ENDIF
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    4. APPLY CORRECTION TO WAVE HEIGHT DATA ACCORDING TO SROKOSZ.
!        ---------------------------------------------------------

      NERS = 0

!     4.1 DETERMINE NUMBER OF OBSERVATIONS ON ONE PE.
!         ------------------------------------------
      
      IF(LMESSPASS) THEN
        DO IJ = NSTART(IRANK), NEND(IRANK)
          IX = IXLG(IJ,1)
          IY = NGY- KXLT(IJ,1) +1
          IF (WHME(IX,IY).GT.0.) NERS = NERS+1
        ENDDO
      ELSE
        DO IG = 1, IGL
          DO IJ = IJS(IG),IJL(IG)
            IX = IXLG(IJ,IG)
            IY = NGY- KXLT(IJ,IG) +1
            IF (WHME(IX,IY).GT.0.) NERS = NERS+1
          ENDDO
        ENDDO
      ENDIF
      
!     4.2 COLLOCATE WAMSPECTRA WITH ALT-HS OBSERVATIONS.
!         ---------------------------------------------
      
      IF (NERS.GT.0) THEN
        ALLOCATE( XKAPPA1(NERS),IJALT(NERS,IGL),DELH_ALT(NERS),
     1            FALT(NERS,NANG,NFRE) )
      
        IERS = 0
        IF(LMESSPASS) THEN
          DO IJ = NSTART(IRANK), NEND(IRANK)
            IX = IXLG(IJ,1)
            IY = NGY- KXLT(IJ,1) +1
            IF (WHME(IX,IY).GT.0.) THEN
              IERS = IERS+1
              IJALT(IERS,1)=IJ
            ENDIF   
          ENDDO
          DO IERS = 1, NERS
            FALT(IERS,:,:)=FL1(IJALT(IERS,1),:,:)
          ENDDO
        ELSE
          DO IG = 1, IGL
            DO IJ = IJS(IG),IJL(IG)
              IX = IXLG(IJ,IG)
              IY = NGY- KXLT(IJ,IG) +1
              IF (WHME(IX,IY).GT.0.) THEN
                IERS = IERS+1
                IJALT(IERS,IG)=IJ
                FALT(IERS,:,:)=FL1(IJ,:,:)
              ENDIF   
            ENDDO
          ENDDO
        ENDIF
      
        CALL SKEWNESS(IU06,FALT,NERS,NANG,NFRE,FR(1),TH,COSTH,SINTH,
     1                XKAPPA1,DELH_ALT)
      
      ENDIF

!     4.3 APPLY SROKOSZ CORRECTION.
!         -------------------------
     
      WHGTAG = ZMISS
      CWHTAG = ZMISS
      RANGAG = ZMISS
 
      IF(LALTCOR.AND.ITEST.GT.0) THEN
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' ALTIMETER WAVE HEIGHT CORRECTION '
        WRITE(IU06,*) ' -------------------------------- ' 
        WRITE(IU06,*) ' '
      ENDIF

      IF(LMESSPASS) THEN
        DO IR = 1, NPROC
          IF(IR.EQ.IRANK) KERS=NERS
          CALL MPBCASTINTFLD(IR,KTAG,1,1,KERS)
          KTAG=KTAG+1
          IF(KERS.GT.0) THEN
            ALLOCATE(XKAP(KERS))
            IF(IR.EQ.IRANK) XKAP=XKAPPA1
            CALL MPBCASTGRDFLD(IR,KTAG,KERS,1,XKAP)
            KTAG=KTAG+1
            ALLOCATE(DELH(KERS))
            IF(IR.EQ.IRANK) DELH=DELH_ALT
            CALL MPBCASTGRDFLD(IR,KTAG,KERS,1,DELH)
            KTAG=KTAG+1
            ALLOCATE(IJKAP(KERS))
            IF(IR.EQ.IRANK) IJKAP=IJALT(:,1)
            CALL MPBCASTINTFLD(IR,KTAG,KERS,1,IJKAP)
            KTAG=KTAG+1
            DO IERS = 1, KERS 
              IJ = IJKAP(IERS)
              IX = IXLG(IJ,1)
              IY = NGY- KXLT(IJ,1) +1
              XK = XKAP(IERS)
              RCOR = DELH(IERS) 
              HS = WHME(IX,IY)
              HS_NEW = SWH_NEW(HS,XK,XKAPPA2)
              
              IF(HS.GT.0.) THEN
                IF(ALLOCATED(WHGTAG)) WHGTAG(IX,IY) = HS 
                IF(ALLOCATED(CWHTAG)) CWHTAG(IX,IY) = HS_NEW
                IF(ALLOCATED(RANGAG)) RANGAG(IX,IY) = RCOR 
              ENDIF

              IF(LALTCOR) THEN
                WHME(IX,IY) = HS_NEW 
                IF(ITEST.GT.0) THEN
                  XLON = AMOWEP + REAL(IX-1)*ZDELLO(KXLT(IJ,1))
                  XLAT = AMONOP - REAL(IY-1)*XDELLA
                  WRITE(IU06,*) XLON, XLAT, XK, HS, WHME(IX,IY),
     &                        WHGTTG(IX,IY), USME(IX,IY), U10G(IX,IY)
                ENDIF
              ENDIF
            ENDDO
            DEALLOCATE(XKAP)
            DEALLOCATE(DELH)
            DEALLOCATE(IJKAP)
          ENDIF
        ENDDO
      ELSE
        DO IG = 1, IGL
          DO IERS = 1, NERS 
            IJ = IJALT(IERS,IG)
            IF(IJ.GE.IJS(IG).AND.IJ.LE.IJL(IG)) THEN
              XK = XKAPPA1(IERS)
              RCOR = DELH_ALT(IERS)
              IX = IXLG(IJ,1)
              IY = NGY- KXLT(IJ,1) +1
              HS = WHME(IX,IY)
              HS_NEW = SWH_NEW(HS,XK,XKAPPA2)

              IF(HS.GT.0.) THEN
                IF(ALLOCATED(WHGTAG)) WHGTAG(IX,IY) = HS 
                IF(ALLOCATED(CWHTAG)) CWHTAG(IX,IY) = HS_NEW 
                IF(ALLOCATED(RANGAG)) RANGAG(IX,IY) = RCOR 
              ENDIF

              IF(LALTCOR) WHME(IX,IY) = HS_NEW 
            ENDIF
          ENDDO
        ENDDO
      ENDIF

   
      IF( ALLOCATED(XKAPPA1) ) DEALLOCATE(XKAPPA1)
      IF( ALLOCATED(IJALT) ) DEALLOCATE(IJALT)
      IF( ALLOCATED(DELH_ALT) ) DEALLOCATE(DELH_ALT)
      IF( ALLOCATED(FALT) ) DEALLOCATE(FALT)
      

!*    5. OBTAIN FRICTION VELOCITIES FROM OBSERVED WINDS AND THE DRAG.
!        ------------------------------------------------------------

      FV: DO J=1,NGY
        DO I=1,NGX
          IF (USME(I,J) .GT. 0) USME(I,J) = SQRT(CDG(I,J)) * USME(I,J)
          IF (CDG(I,J) .LT. 0.000) THEN
            USME(I,J)=-3.
            WHME(I,J)=-3.
          ENDIF
        ENDDO
      ENDDO FV

      RETURN
      END SUBROUTINE GRFIELD
      
      FUNCTION SWH_NEW(SWH,XKAPPA1,K2)
      
      REAL:: SWH, SWH_NEW, XN, XKAPPA1, K2
      
      XN = SWH**2/16. + K2
      
      SWH_NEW = 16.*(XKAPPA1*XN - K2)
      
      IF ( SWH_NEW >  0.0 ) THEN
        SWH_NEW = SQRT(SWH_NEW)
      ELSE
        SWH_NEW = 0.0
      ENDIF

      RETURN
      END FUNCTION SWH_NEW
