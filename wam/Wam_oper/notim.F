      SUBROUTINE NOTIM (CDTWIS, CDTWIE, IREAD, LWCUR)

! ----------------------------------------------------------------------

!     MODIFIED  J. BIDLOT  FEBRARY 1996  MESSAGE PASSING
!               S. ABDALLA OCTOBER 2001  INCLUSION OF AIR DENSITY & Zi/L
!               J BIDLOT AUGUST 2006 INCLUSION OF ICE MASK.
!               J BIDLOT AUGUST 2008 REMOVE INFORMATION FROM COUPLED RUNS
!                                    AND CLEAN UP OBSOLETE OPTIONS.


!**** *NOTIM* - STEERING MODULE IF NO TIME INTERPOLATION WANTED.

!*    PURPOSE.
!     --------

!       NOTIM NO TIME INTERPOLATION: PROCESS WINDFIELDS.

!**   INTERFACE.
!     ----------

!       *CALL* *NOTIM (CDTWIS, CDTWIE, IREAD, LWCUR)
!          *CDTWIS* - DATE OF FIRST WIND FIELD.
!          *CDTWIE* - DATE OF LAST FIRST WIND FIELD.
!          *IREAD*  - PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!          *LWCUR*  - LOGICAL INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS


!     METHOD.
!     -------

!       NO TIME INTERPOLATION:
!       WINDFIELDS ARE PROCESSED EVERY IDELWI SECONDS (U,V),
!       THE WINDS INTERPOLATED IN SPACE ONLY (US,DS)

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.
!       *AIRSEA*    - TOTAL STRESS IN SURFACE LAYER.
!       *CREWFN*    - ASSIGN A FILE NAME.
!       *GETWND*    - READ A WINDFIELD AND COMPUTE WIND
!                     FOR ALL BLOCKS (US,DS).
!       *INCDATE*   - INCREMENT DATE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWICE   , ONLY : CICOVER  ,CITHICK
      USE YOWPARAM , ONLY : NIBLO    ,NBLO
      USE YOWCOUP  , ONLY : LWCOU
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWSPEC  , ONLY : U10OLD   ,THWOLD   ,USOLD    ,Z0OLD    ,
     &              TAUW   ,ROAIRO   ,ZIDLOLD
      USE YOWSTAT  , ONLY : IDELPRO  ,IDELWO   ,NWAM_BLKS
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWIND  , ONLY : CDA      ,CDTNEXT  ,NSTORE   ,
     &            USNEXT   ,DSNEXT   ,ADSNEXT ,ZIDLNEXT  ,
     &            CICOVERNEXT, CITHICKNEXT

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      INTEGER :: JKGLO,KIJS,KIJL,NPROMA

      REAL,DIMENSION(NINF:NSUP,NBLO) :: US, DS, ADS, ZIDL, CICR, CITH
      REAL :: ZHOOK_HANDLE

      CHARACTER*14 :: CDTWIE, CDTWIH, CDTWIS, CDT, ZERO
      LOGICAL :: LWNDFILE, LCLOSEWND
      LOGICAL :: LWCUR

! ----------------------------------------------------------------------

!*    1. INITIALISE WIND REQUEST DATE AND PROCESS FIRST WINDFIELD
!*       IF COLD START.
!        --------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('NOTIM',0,ZHOOK_HANDLE)

      ZERO = ' '
      CDTWIH = CDTWIS

      LCLOSEWND=.FALSE.
      IF (LWCOU) THEN
        LWNDFILE=.FALSE.
      ELSE
        LWNDFILE=.TRUE.
      ENDIF

      ILEV=1

      IF (CDA.EQ.ZERO) THEN
        CDA = CDTWIS

        CALL GETWND (U10OLD, THWOLD, ROAIRO, ZIDLOLD, CICOVER, CITHICK,
     &               CDA, LWNDFILE, LCLOSEWND, IREAD,
     &               LWCUR)

        DO IG = 1,IGL
          NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1
          CALL GSTATS(1493,0)
!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
          DO JKGLO=IJS(IG),IJL(IG),NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
            CALL AIRSEA (U10OLD(KIJS,IG), TAUW(KIJS,IG),
     &        USOLD(KIJS,IG), Z0OLD(KIJS,IG),
     &        KIJS, KIJL, ILEV)
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1493,1)
        ENDDO
        IF (ITEST.GE.3) THEN
          WRITE(IU06,'(''       SUB. NOTIM: FIRST WIND FIELD '',
     &     ''SAVED IN COMMON WIND'')')
        ENDIF
        IF (CDA.EQ.CDTWIE) THEN
          IF (LHOOK) CALL DR_HOOK('NOTIM',1,ZHOOK_HANDLE)
          RETURN
        ENDIF
        CALL INCDATE (CDTWIH,IDELWO)
      ENDIF

! ----------------------------------------------------------------------

!*    3. LOOP OVER OUTPUT WIND TIMES.
!        ----------------------------

!*    3.1 READ ONE WIND FIELD AND TRANSFORM TO BLOCKS.
!         --------------------------------------------

      CALL GETWND (US, DS, ADS, ZIDL, CICR, CITH,
     &             CDTWIH, LWNDFILE, LCLOSEWND, IREAD,
     &             LWCUR)

!*    3.2 SAVE BLOCKED WIND FIELD.
!         ------------------------

      IF (IDELPRO.GT.IDELWO) THEN
!*      WIND FIELD IS NOT CONSTANT FOR ONE PROPAGATION TIME STEP:
!       THIS OPTION IS NO LONGER AVAILABLE !!!
        WRITE (IU06,*) ' '
        WRITE (IU06,*) ' **********************************'
        WRITE (IU06,*) ' *                                *'
        WRITE (IU06,*) ' *   FATAL ERROR IN SUB. NOTIM:   *'
        WRITE (IU06,*) ' *   ==========================   *'
        WRITE (IU06,*) ' *   THE OPTION IDELPRO.GT.IDELWO *'
        WRITE (IU06,*) ' *      IS LONGER AVAILABLE !!!!! *'
        WRITE (IU06,*) ' *                                *'
        WRITE (IU06,*) ' *          PROGRAM ABORTS        *'
        WRITE (IU06,*) ' *                                *'
        WRITE (IU06,*) ' **********************************'
        CALL ABORT1
      ELSE
!*      WIND FIELD IS CONSTANT FOR ONE PROPAGATION TIME STEP:
!       -----------------------------------------------------
        NSTORE=1
        IF(.NOT.ALLOCATED(CDTNEXT)) ALLOCATE(CDTNEXT(NSTORE))
        IF(LMESSPASS) THEN
          IF(.NOT.ALLOCATED(USNEXT))
     &            ALLOCATE(USNEXT(IJS(1):IJL(1),NSTORE,NBLO))
          IF(.NOT.ALLOCATED(DSNEXT))
     &            ALLOCATE(DSNEXT(IJS(1):IJL(1),NSTORE,NBLO))
          IF(.NOT.ALLOCATED(ADSNEXT))
     &            ALLOCATE(ADSNEXT(IJS(1):IJL(1),NSTORE,NBLO))
          IF(.NOT.ALLOCATED(ZIDLNEXT))
     &            ALLOCATE(ZIDLNEXT(IJS(1):IJL(1),NSTORE,NBLO))
          IF(.NOT.ALLOCATED(CICOVERNEXT))
     &            ALLOCATE(CICOVERNEXT(IJS(1):IJL(1),NSTORE,NBLO))
          IF(.NOT.ALLOCATED(CITHICKNEXT))
     &            ALLOCATE(CITHICKNEXT(IJS(1):IJL(1),NSTORE,NBLO))
        ELSE
          IF(.NOT.ALLOCATED(USNEXT))
     &            ALLOCATE(USNEXT(NIBLO,NSTORE,NBLO))
          IF(.NOT.ALLOCATED(DSNEXT))
     &            ALLOCATE(DSNEXT(NIBLO,NSTORE,NBLO))
          IF(.NOT.ALLOCATED(ADSNEXT))
     &            ALLOCATE(ADSNEXT(NIBLO,NSTORE,NBLO))
          IF(.NOT.ALLOCATED(ZIDLNEXT))
     &            ALLOCATE(ZIDLNEXT(NIBLO,NSTORE,NBLO))
          IF(.NOT.ALLOCATED(CICOVERNEXT))
     &            ALLOCATE(CICOVERNEXT(NIBLO,NSTORE,NBLO))
          IF(.NOT.ALLOCATED(CITHICKNEXT))
     &            ALLOCATE(CITHICKNEXT(NIBLO,NSTORE,NBLO))
        ENDIF

        DO IG=1,IGL
          CDTNEXT(1)=CDTWIH

          NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1
          CALL GSTATS(1493,0)
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ)
          DO JKGLO=IJS(IG),IJL(IG),NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
            DO IJ = KIJS, KIJL
              USNEXT(IJ,1,IG)=US(IJ,IG)
              DSNEXT(IJ,1,IG)=DS(IJ,IG)
              ADSNEXT(IJ,1,IG)=ADS(IJ,IG)
              ZIDLNEXT(IJ,1,IG)=ZIDL(IJ,IG)
              CICOVERNEXT(IJ,1,IG)=CICR(IJ,IG)
              CITHICKNEXT(IJ,1,IG)=CITH(IJ,IG)
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1493,1)
        ENDDO

!*      UPDATE WIND FIELD REQUEST TIME.
        CALL INCDATE (CDTWIH,IDELWO)

!*      IF TIME LEFT BRANCH NOT ALLOWED TO LOOP ANYMORE 
        IF (CDTWIH.LE.CDTWIE) THEN
          WRITE (IU06,*) ' '
          WRITE (IU06,*) ' ********************************************'
          WRITE (IU06,*) ' *                                          *'
          WRITE (IU06,*) ' *       FATAL ERROR IN SUB. NOTIM:         *'
          WRITE (IU06,*) ' *       ==========================         *'
          WRITE (IU06,*) ' * ERROR WHEN WRITTING ON WIND OUTPUT FILE  *'
          WRITE (IU06,*) ' * MORE THAN ONE LOOP NOT ALLOWED ANY LONGER*'
          WRITE (IU06,*) ' * WHEN USING ARRAY INSTEAD OF WRITE TO FILE*'
          WRITE (IU06,*) ' *                                          *'
          WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS         *'
          WRITE (IU06,*) ' *                                          *'
          WRITE (IU06,*) ' ********************************************'
          CALL ABORT1
        ENDIF

        IF (ITEST.GE.3) THEN
          WRITE(IU06,*) '      SUB. NOTIM: NEW FORCING FIELDS GENERATED'
        ENDIF

        IF (LHOOK) CALL DR_HOOK('NOTIM',1,ZHOOK_HANDLE)
        RETURN

      ENDIF

      IF (LHOOK) CALL DR_HOOK('NOTIM',1,ZHOOK_HANDLE)
      RETURN


      END SUBROUTINE NOTIM
