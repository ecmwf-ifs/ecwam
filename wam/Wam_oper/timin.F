      SUBROUTINE TIMIN (CDTWIS, CDTWIE,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     &                  ROAIRO, ZIDLOLD, CICOVER, IREAD,
     &                  LWCUR)
! ----------------------------------------------------------------------
!     MODIFIED  J. BIDLOT  FEBRARY 1996  MESSAGE PASSING
!               S. ABDALLA OCTOBER 2001  INCLUSION OF AIR DENSITY & Zi/L
!               J BIDLOT AUGUST 2006 INCLUSION OF ICE MASK.
!               J BIDLOT AUGUST 2008 REMOVE INFORMATION FROM COUPLED RUNS
!                                    AND CLEAN UP OBSOLETE OPTIONS.



!**** *TIMIN* - STEERING MODULE IF TIME INTERPOLATION IS WANTED.

!     H. GUNTHER    ECMWF    MAY 1990     MODIFIED FOR SUB VERSION.
!     H. GUNTHER    ECMWF    DECEMBER 90  MODIFIED FOR CYCLE_4.

!*    PURPOSE.
!     --------

!       TIME INTERPOLATION: PROCESS WINDFIELDS.

!**   INTERFACE.
!     ----------

!       *CALL* *TIMIN (CDTWIS, CDTWIE,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
!    &                 ROAIRO, ZIDLOLD, CICOVER, IREAD,
!    &                 LWCUR)*
!          *CDTWIS* - DATE OF FIRST WIND FIELD.
!          *CDTWIE* - DATE OF LAST WIND FIELD.
!          *U10OLD* - WIND SPEED.
!          *THWOLD* - WIND DIRECTION (RADIANS).
!          *USOLD*  - FRICTION VELOCITY.
!          *TAUW*   - WAVE STRESS.
!          *Z0OLD*  - ROUGHNESS LENGTH IN M.
!          *ROAIRO* - AIR DENSITY IN KG/M3.
!          *ZIDLOLD*- Zi/L
!                     (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!          *CICOVER*- SEA ICE COVER.
!          *IREAD*  - PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!          *LWCUR*  - LOGICAL INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS


!     METHOD.
!     -------

!       WINDFIELDS ARE READ IN EVERY IDELWI SECONDS (U,V),
!       INTERPOLATED IN SPACE, BLOCKED AND SAVED ON SCRATCH UNITS.
!       MAGNITUDE AND DIRECTION ARE INTERPOLATED LINEARLY IN TIME AND
!       WRITTEN TO THE OUTPUT FILES.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.
!       *AIRSEA*    - TOTAL STRESS IN SURFACE LAYER.
!       *CREWFN*    - ASSIGN A FILE NAME.
!       *GETWND*    - READ A WINDFIELD AND COMPUTE WIND
!                     FOR ALL BLOCKS (US,DS).
!       *INCDATE*   - INCREMENT DATE.
!       *LOCINT*    - SPACE INTERPOLATION.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWGRIBHD, ONLY : NKSEK1
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NIBLO    ,NBLO
      USE YOWPCONS , ONLY : PI       ,ZPI
      USE YOWSTAT  , ONLY : IDELPRO  ,IDELWI   ,IDELWO   ,NWAM_BLKS
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWIND  , ONLY : CDA      ,CDTNEXT  ,NSTORE   ,
     &            USNEXT   ,DSNEXT   ,ADSNEXT ,ZIDLNEXT  ,CICOVERNEXT

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      INTEGER :: JKGLO,KIJS,KIJL,NPROMA

      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD,USOLD,Z0OLD,TAUW
     &                                , ROAIRO,ZIDLOLD,CICOVER
      REAL,DIMENSION(NINF:NSUP,NBLO) :: US, DS   ,ADS ,ZIDL, CICR
      REAL,DIMENSION(NINF:NSUP,NBLO) :: US2, DS2 ,ADS2,ZIDL2, CICR2
      REAL,DIMENSION(NINF:NSUP) :: US3, DS3, ADS3, ZIDL3, CICR3
      REAL ZHOOK_HANDLE

      CHARACTER*14 :: CDTWIE, CDTH, CDT1, CDT2, CDTWIS, ZERO
      CHARACTER*80 :: FILENAME

      LOGICAL :: LWNDFILE, LCLOSEWND
      LOGICAL :: LWCUR

! ----------------------------------------------------------------------

!*    1. INITIALIZE TIMECOUNTER.
!        -----------------------

      IF (LHOOK) CALL DR_HOOK('TIMIN',0,ZHOOK_HANDLE)

      ZERO= ' '

      LCLOSEWND=.FALSE.
      IF (LWCOU) THEN
        LWNDFILE=.FALSE.
      ELSE
        LWNDFILE=.TRUE.
      ENDIF
      ILEV=1

      IF (CDA.EQ.ZERO) THEN
        CDT1 = CDTWIS
        CALL GETWND (US, DS, ADS, ZIDL, CICR,
     &               CDT1, LWNDFILE, LCLOSEWND, IREAD,
     &               LWCUR)
        DO IG=1,IGL
          NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1
          CALL GSTATS(1494,0)
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ)
          DO JKGLO=IJS(IG),IJL(IG),NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL(IG))

            DO IJ=KIJS,KIJL
              U10OLD(IJ,IG) = US(IJ,IG)
              THWOLD(IJ,IG) = DS(IJ,IG)
              ROAIRO(IJ,IG) = ADS(IJ,IG)
              ZIDLOLD(IJ,IG)= ZIDL(IJ,IG)
              CICOVER(IJ,IG)= CICR(IJ,IG)
            ENDDO

            CALL AIRSEA (U10OLD(KIJS,IG), TAUW(KIJS,IG),
     &        USOLD(KIJS,IG), Z0OLD(KIJS,IG),
     &        KIJS, KIJL, ILEV)
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1494,1)
        ENDDO
        IF (ITEST.GE.3) THEN
          WRITE(IU06,'(''       SUB. TIMIN: FIRST WIND FIELD '',
     &     ''SAVED IN COMMON WIND'')')
          CALL FLUSH(IU06)
        ENDIF
      ELSE
        CDT1 = CDA
        DO IG=1,IGL
          NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1
          CALL GSTATS(1494,0)
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ)
          DO JKGLO=IJS(IG),IJL(IG),NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
            DO IJ=KIJS,KIJL
              US(IJ,IG) = U10OLD(IJ,IG)
              DS(IJ,IG) = THWOLD(IJ,IG)
              ADS(IJ,IG) = ROAIRO(IJ,IG)
              ZIDL(IJ,IG)= ZIDLOLD(IJ,IG)
              CICR(IJ,IG)= CICOVER(IJ,IG)
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1494,1)
        ENDDO
        CDTH = CDT1
        CALL INCDATE (CDTH,IDELWO)
        IF (CDTWIS.NE.CDTH) THEN
          WRITE(IU06,*) ' *******************************************'
          WRITE(IU06,*) ' *                                         *'
          WRITE(IU06,*) ' *        FATAL ERROR IN --TIMIN--         *'
          WRITE(IU06,*) ' *        =========================        *'
          WRITE(IU06,*) ' * DATES DO NOT MATCH.                     *'
          WRITE(IU06,*) ' * START DATE FOR WIND IS       CDTWIS = ',
     &     CDTWIS
          WRITE(IU06,*) ' * LAST DATE SAVED IN COM WIND IS CDT1 = ',
     &     CDT1
          WRITE(IU06,*) ' * PROCESSING WILL BE ABORTED              *'
          WRITE(IU06,*) ' *                                         *'
          WRITE(IU06,*) ' *******************************************'
          CALL ABORT1
        ENDIF
      ENDIF
      CDT2 = CDT1
      CALL INCDATE(CDT2,IDELWI)

      NTS = IDELWI/IDELWO
      DEL = REAL(IDELWO)/REAL(IDELWI)

!*    2. LOOP OVER INPUT WINDFIELDS.
!        ---------------------------

 2000 CONTINUE

!*    2.1 READ ONE WINDFIELD AND TRANSFORM ALL BLOCKS.
!         -------------------------------------------

      CDT2 = CDT1
      CALL INCDATE(CDT2,IDELWI)
      CALL GETWND (US2, DS2, ADS2, ZIDL2, CICR2,
     &             CDT2, LWNDFILE, LCLOSEWND, IREAD,
     &             LWCUR)

!*    2.2 SAVE BLOCKED WIND FIELD ON SCRATCH UNITS.
!         -----------------------------------------

      NSTORE=NTS

      IF(.NOT.ALLOCATED(CDTNEXT)) ALLOCATE(CDTNEXT(NSTORE))
      IF(LMESSPASS) THEN
        IF(.NOT.ALLOCATED(USNEXT))
     &          ALLOCATE(USNEXT(IJS(1):IJL(1),NSTORE,NBLO))
        IF(.NOT.ALLOCATED(DSNEXT))
     &          ALLOCATE(DSNEXT(IJS(1):IJL(1),NSTORE,NBLO))
        IF(.NOT.ALLOCATED(ADSNEXT))
     &          ALLOCATE(ADSNEXT(IJS(1):IJL(1),NSTORE,NBLO))
        IF(.NOT.ALLOCATED(ZIDLNEXT))
     &          ALLOCATE(ZIDLNEXT(IJS(1):IJL(1),NSTORE,NBLO))
        IF(.NOT.ALLOCATED(CICOVERNEXT))
     &          ALLOCATE(CICOVERNEXT(IJS(1):IJL(1),NSTORE,NBLO))
      ELSE
        IF(.NOT.ALLOCATED(USNEXT))
     &          ALLOCATE(USNEXT(NIBLO,NSTORE,NBLO))
        IF(.NOT.ALLOCATED(DSNEXT))
     &          ALLOCATE(DSNEXT(NIBLO,NSTORE,NBLO))
        IF(.NOT.ALLOCATED(ADSNEXT))
     &          ALLOCATE(ADSNEXT(NIBLO,NSTORE,NBLO))
        IF(.NOT.ALLOCATED(ZIDLNEXT))
     &          ALLOCATE(ZIDLNEXT(NIBLO,NSTORE,NBLO))
        IF(.NOT.ALLOCATED(CICOVERNEXT))
     &          ALLOCATE(CICOVERNEXT(NIBLO,NSTORE,NBLO))
      ENDIF

      DO IG=1,IGL
        CDTH = CDT1
        DO N=1,NTS
          CALL INCDATE(CDTH,IDELWO)
          CDTNEXT(N)=CDTH
          NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1
          CALL GSTATS(1494,0)
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ,D)
          DO JKGLO=IJS(IG),IJL(IG),NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
            DO IJ=KIJS,KIJL
              US3(IJ) = US(IJ,IG)+REAL(N)*DEL*(US2(IJ,IG)-US(IJ,IG))
              D = DS2(IJ,IG)-DS(IJ,IG)
              IF (ABS(D).GT.PI) D = D-ZPI*SIGN(1.,D)
              D = DS(IJ,IG)+REAL(N)*DEL*D
              DS3(IJ) = MOD(D+ZPI,ZPI)
              ADS3(IJ) = ADS(IJ,IG)+
     &                        REAL(N)*DEL*(ADS2(IJ,IG)-ADS(IJ,IG))
              ZIDL3(IJ) = ZIDL(IJ,IG)+
     &                         REAL(N)*DEL*(ZIDL2(IJ,IG)-ZIDL(IJ,IG))
              CICR3(IJ) = CICR(IJ,IG)+
     &                         REAL(N)*DEL*(CICR2(IJ,IG)-CICR(IJ,IG))

              USNEXT(IJ,N,IG)=US3(IJ)
              DSNEXT(IJ,N,IG)=DS3(IJ)
              ADSNEXT(IJ,N,IG)=ADS3(IJ)
              ZIDLNEXT(IJ,N,IG)=ZIDL3(IJ)
              CICOVERNEXT(IJ,N,IG)=CICR3(IJ)
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1494,1)
        ENDDO

        DO IJ=IJS(IG),IJL(IG)
          US(IJ,IG) = US2(IJ,IG)
          DS(IJ,IG) = DS2(IJ,IG)
          ADS(IJ,IG) = ADS2(IJ,IG)
          ZIDL(IJ,IG) = ZIDL2(IJ,IG)
          CICR(IJ,IG) = CICR2(IJ,IG)
        ENDDO

      ENDDO
      CDT1 = CDT2
      IF (ITEST.GE.3) THEN
        WRITE(IU06,*) '       SUB. TIMIN: LAST WIND FIELD AT ',
     &   'CDTH = ', CDTH,' WRITTEN TO SCRATCH UNITS'
          CALL FLUSH (IU06)
      ENDIF

!*    2.3 UPDATE WIND FIELD REQUEST TIME AND READ NEXT IF REQUESTED.
!         ----------------------------------------------------------

      CALL INCDATE (CDTH,IDELWI)
!!! we don't allow to loop anymore
!!!!      IF (CDTH.LE.CDTWIE) GOTO 2000


      IF (ITEST.GE.3) THEN
        WRITE(IU06,*) '       SUB. TIMIN: LAST WIND FIELD AT ',
     &   'CDTH = ', CDTH,' WRITTEN TO OUTPUT UNIT'
          CALL FLUSH (IU06)
      ENDIF

      IF (LHOOK) CALL DR_HOOK('TIMIN',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE TIMIN
