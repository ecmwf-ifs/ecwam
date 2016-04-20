      SUBROUTINE NOTIM (CDTWIS, CDTWIE,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     &                  ROAIRO, ZIDLOLD, CICOVER, CITHICK, IREAD,
     &                  LWCUR)


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

!       *CALL* *NOTIM (CDTWIS, CDTWIE,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
!    &                 ROAIRO, ZIDLOLD, CICOVER, CITHICK, IREAD,
!                      LWCUR)
!          *CDTWIS* - DATE OF FIRST WIND FIELD.
!          *CDTWIE* - DATE OF LAST FIRST WIND FIELD.
!          *U10OLD* - WIND SPEED.
!          *THWOLD* - WIND DIRECTION (RADIANS).
!          *USOLD*  - FRICTION VELOCITY.
!          *TAUW*   - WAVE STRESS.
!          *Z0OLD*  - ROUGHNESS LENGTH IN M.
!          *ROAIRO* - AIR DENSITY IN KG/M3.
!          *ZIDLOLD*- Zi/L 
!                     (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!          *CICOVER*- SEA ICE COVER.
!          *CITHICK*- SEA ICE THICKNESS. 
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
!       *GETWND*    - READ A WINDFIELD AND COMPUTE WIND
!                     FOR ALL BLOCKS (US,DS).
!       *INCDATE*   - INCREMENT DATE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NIBLO    ,NBLO
      USE YOWCOUP  , ONLY : LWCOU
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWSTAT  , ONLY : IDELPRO  ,IDELWO   ,NPROMA_WAM
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWUNIT  , ONLY : IU25     ,IU26
      USE YOWWIND  , ONLY : CDA      ,CDTNEXT  ,NSTORE   ,
     &            USNEXT   ,DSNEXT   ,ADSNEXT ,ZIDLNEXT  ,
     &            CICOVERNEXT, CITHICKNEXT
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IREAD
      
      REAL,DIMENSION(NINF:NSUP,NBLO), INTENT(INOUT) ::
     &               U10OLD,THWOLD,USOLD,Z0OLD,TAUW,
     &               ROAIRO,ZIDLOLD,CICOVER,CITHICK

      CHARACTER(LEN=14), INTENT(IN) :: CDTWIS, CDTWIE

      LOGICAL, INTENT(IN) :: LWCUR


      INTEGER :: JKGLO,KIJS,KIJL,NPROMA
      INTEGER :: IG, ILEV, IJ

      REAL :: ZHOOK_HANDLE
      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10, THW, ADS, ZIDL, CICR, CITH

      CHARACTER(LEN=14) :: CDTWIH, CDT, ZERO

      LOGICAL :: LWNDFILE, LCLOSEWND

! ----------------------------------------------------------------------

!*    1. INITIALISE WIND REQUEST DATE AND PROCESS FIRST WINDFIELD
!*       IF COLD START.
!        --------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('NOTIM',0,ZHOOK_HANDLE)
#endif

      ZERO = ' '
      CDTWIH = CDTWIS

      NPROMA=NPROMA_WAM

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
#ifdef ECMWF
          IF (LHOOK) CALL DR_HOOK('NOTIM',1,ZHOOK_HANDLE)
#endif
          RETURN
        ENDIF
        CALL INCDATE (CDTWIH,IDELWO)
      ENDIF

! ----------------------------------------------------------------------

!*    3. LOOP OVER OUTPUT WIND TIMES.
!        ----------------------------

!*    3.1 READ ONE WIND FIELD AND TRANSFORM TO BLOCKS.
!         --------------------------------------------

      CALL GETWND (U10, THW, ADS, ZIDL, CICR, CITH,
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

          CALL GSTATS(1493,0)
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ)
          DO JKGLO=IJS(IG),IJL(IG),NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
            DO IJ = KIJS, KIJL
              USNEXT(IJ,1,IG)=U10(IJ,IG)
              DSNEXT(IJ,1,IG)=THW(IJ,IG)
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
#ifdef ECMWF
        IF (LHOOK) CALL DR_HOOK('NOTIM',1,ZHOOK_HANDLE)
#endif
        RETURN

      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('NOTIM',1,ZHOOK_HANDLE)
#endif
      RETURN
      END SUBROUTINE NOTIM
