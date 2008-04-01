      SUBROUTINE TIMIN (CDTWIS, CDTWIE,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     &                  ROAIRO, ZIDLOLD, ICEMASK, 
     &                  NFIELDS, NGPTOTG, NC, NR, FCRANGE,
     &                  ISEC1, ISEC2, FIELDS, MASK_IN)

! ----------------------------------------------------------------------
!     MODIFIED  J. BIDLOT  FEBRARY 1996  MESSAGE PASSING
!               S. ABDALLA OCTOBER 2001  INCLUSION OF AIR DENSITY & Zi/L
!               J BIDLOT AUGUST 2006 INCLUSION OF ICE MASK.


!**** *TIMIN* - STEERING MODULE IF TIME INTERPOLATION IS WANTED.

!     H. GUNTHER    ECMWF    MAY 1990     MODIFIED FOR SUB VERSION.
!     H. GUNTHER    ECMWF    DECEMBER 90  MODIFIED FOR CYCLE_4.

!*    PURPOSE.
!     --------

!       TIME INTERPOLATION: PROCESS WINDFIELDS.

!**   INTERFACE.
!     ----------

!       *CALL* *TIMIN (CDTWIS, CDTWIE,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
!    &                 ROAIRO, ZIDLOLD, ICEMASK,
!    &                 NFIELDS, NGPTOTG, NC, NR, FCRANGE,
!    &                 ISEC1, ISEC2, FIELDS)*
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
!          *ICEMASK*- SEA ICE MASK.
!          *NFIELDS*- NUMBER OF FIELDS HOLDING ATMOSPHERIC DATA
!          *NGPTOTG*- NUMBER OF ATMOSPHERIC GRID POINTS
!          *NC*     - NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!          *NR*     - NUMBER OF ATM. ROWS OF LATITUDES
!          *FCRANGE*- FORECAST RANGE IN SECONDS
!          *ISEC1*  - GRIB ISEC1 (TO EXTRACT DATE/TIME OF ATM. FIELDS)
!          *ISEC2*  - GRIB ISEC2 (TO EXTRACT ATM. GRID CHARACTERISTICS)
!          *FIELDS* - ATM. FIELDS (U10, V10, Q2, MSL, SKT, ice, Zi/L)

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
      USE YOWSTAT  , ONLY : IDELPRO  ,IDELWI   ,IDELWO
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWUNIT  , ONLY : IUVELO   ,IUSCR
      USE YOWWIND  , ONLY : CDA      ,CDTNEXT  ,NSTORE   ,
     &            USNEXT   ,DSNEXT   ,ADSNEXT ,ZIDLNEXT  ,ICEMASKNEXT

! ----------------------------------------------------------------------

      CHARACTER*14 CDTWIE, CDTH, CDT1, CDT2, CDTWIS, ZERO
      LOGICAL LWNDFILE, LCLOSEWND

      INTEGER NFIELDS, NGPTOTG, NC, NR
      INTEGER ISEC1(NKSEK1), ISEC2(22+NR)
      REAL FIELDS(NGPTOTG,NFIELDS)
      INTEGER MASK_IN(NGPTOTG)

!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD,USOLD,Z0OLD,TAUW
     &                                , ROAIRO,ZIDLOLD
      INTEGER,DIMENSION(NINF:NSUP,NBLO) :: ICEMASK 

      CHARACTER*80 FILENAME
! ----------------------------------------------------------------------

      REAL,DIMENSION(NINF:NSUP,NBLO) :: US, DS   ,ADS ,ZIDL
      REAL,DIMENSION(NINF:NSUP,NBLO) :: US2, DS2 ,ADS2,ZIDL2
      REAL,DIMENSION(NINF:NSUP) :: US3, DS3, ADS3, ZIDL3
      INTEGER,DIMENSION(NINF:NSUP,NBLO) :: ICMSK, ICMSK2
      INTEGER,DIMENSION(NINF:NSUP) :: ICMSK3 

!        *US*     REAL   OUTPUT WIND FIELD (FRICTION VELOCITIES).
!        *DS*     REAL   OUTPUT WIND FIELD (DIRECTIONS).

! ----------------------------------------------------------------------

!*    1. INITIALIZE TIMECOUNTER.
!        -----------------------

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
        CALL GETWND (US, DS, ADS, ZIDL, ICMSK,
     &               NINF, NSUP, CDT1,
     &               NFIELDS, NGPTOTG, NC, NR, FCRANGE,
     &               ISEC1, ISEC2, FIELDS, MASK_IN,
     &               LWNDFILE, LCLOSEWND)
        DO IG=1,IGL
          DO IJ=IJS(IG),IJL(IG)
            U10OLD(IJ,IG) = US(IJ,IG)
            THWOLD(IJ,IG) = DS(IJ,IG)
            ROAIRO(IJ,IG) = ADS(IJ,IG)
            ZIDLOLD(IJ,IG)= ZIDL(IJ,IG)
            ICEMASK(IJ,IG)= ICMSK(IJ,IG)
          ENDDO
        ENDDO
        DO IG = 1,IGL
          CALL AIRSEA (U10OLD(IJS(IG),IG), TAUW(IJS(IG),IG),
     &     USOLD(IJS(IG),IG), Z0OLD(IJS(IG),IG),
     &     IJS(IG), IJL(IG), ILEV)

        ENDDO
        IF (ITEST.GE.3) THEN
          WRITE(IU06,'(''       SUB. TIMIN: FIRST WIND FIELD '',
     &     ''SAVED IN COMMON WIND'')')
          CALL FLUSH(IU06)
        ENDIF
      ELSE
        CDT1 = CDA
        DO IG=1,IGL
          DO IJ=IJS(IG),IJL(IG)
            US(IJ,IG) = U10OLD(IJ,IG)
            DS(IJ,IG) = THWOLD(IJ,IG)
            ADS(IJ,IG) = ROAIRO(IJ,IG)
            ZIDL(IJ,IG)= ZIDLOLD(IJ,IG)
            ICMSK(IJ,IG)= ICEMASK(IJ,IG)
          ENDDO
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

!!!!  we will instead use arrays
!!!!      CALL CREWFN (IUVELO, CDTWIE)
!!!!      IF (ITEST.GE.3) THEN
!!!!        WRITE(IU06,*) '       SUB. TIMIN: NEW WIND FILE '
!!!!        WRITE(IU06,*) '       UNIT IS IUVELO = ', IUVELO,
!!!!     &   ' WIND FILE DATE IS CDTWIE = ', CDTWIE
!!!!        CALL FLUSH (IU06)
!!!!      ENDIF
!!!!      REWIND (UNIT=IUVELO)

!*    2. LOOP OVER INPUT WINDFIELDS.
!        ---------------------------

 2000 CONTINUE

!*    2.1 READ ONE WINDFIELD AND TRANSFORM ALL BLOCKS.
!         -------------------------------------------

      CDT2 = CDT1
      CALL INCDATE(CDT2,IDELWI)
      CALL GETWND (US2, DS2, ADS2, ZIDL2, ICMSK2,
     &             NINF, NSUP, CDT2,
     &             NFIELDS, NGPTOTG, NC, NR, FCRANGE,
     &             ISEC1, ISEC2, FIELDS, MASK_IN,
     &             LWNDFILE, LCLOSEWND)

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
        IF(.NOT.ALLOCATED(ICEMASKNEXT))
     &          ALLOCATE(ICEMASKNEXT(IJS(1):IJL(1),NSTORE,NBLO))
      ELSE
        IF(.NOT.ALLOCATED(USNEXT))
     &          ALLOCATE(USNEXT(NIBLO,NSTORE,NBLO))
        IF(.NOT.ALLOCATED(DSNEXT))
     &          ALLOCATE(DSNEXT(NIBLO,NSTORE,NBLO))
        IF(.NOT.ALLOCATED(ADSNEXT))
     &          ALLOCATE(ADSNEXT(NIBLO,NSTORE,NBLO))
        IF(.NOT.ALLOCATED(ZIDLNEXT))
     &          ALLOCATE(ZIDLNEXT(NIBLO,NSTORE,NBLO))
        IF(.NOT.ALLOCATED(ICEMASKNEXT))
     &          ALLOCATE(ICEMASKNEXT(NIBLO,NSTORE,NBLO))
      ENDIF

      DO IG=1,IGL
!!!! use arrays instead of writing to files
!!!!        IUNIT =IUSCR(IG)
!!!!        IF(LMESSPASS) THEN
!!!!          FILENAME='wind_temp_file.%p'
!!!!          CALL EXPAND_STRING(IRANK,NPROC,0,0,FILENAME,1)
!!!!          OPEN(IUNIT,FILE=FILENAME,
!!!!     &         FORM='UNFORMATTED',STATUS='NEW')
!!!!        ENDIF
        CDTH = CDT1
        DO N=1,NTS
          CALL INCDATE(CDTH,IDELWO)
          DO IJ=IJS(IG),IJL(IG)
            US3(IJ) = US(IJ,IG)+REAL(N)*DEL*(US2(IJ,IG)-US(IJ,IG))
            D = DS2(IJ,IG)-DS(IJ,IG)
            IF (ABS(D).GT.PI) D = D-ZPI*SIGN(1.,D)
            D = DS(IJ,IG)+REAL(N)*DEL*D
            DS3(IJ) = MOD(D+ZPI,ZPI)
            ADS3(IJ) = ADS(IJ,IG)+
     &                      REAL(N)*DEL*(ADS2(IJ,IG)-ADS(IJ,IG))
            ZIDL3(IJ) = ZIDL(IJ,IG)+
     &                       REAL(N)*DEL*(ZIDL2(IJ,IG)-ZIDL(IJ,IG))
            ICMSK3(IJ) = NINT(ICMSK(IJ,IG)+
     &                       REAL(N)*DEL*(ICMSK2(IJ,IG)-ICMSK(IJ,IG)))

!!!! arrays replacing the files
            CDTNEXT(N)=CDTH
            USNEXT(IJ,N,IG)=US3(IJ)
            DSNEXT(IJ,N,IG)=DS3(IJ)
            ADSNEXT(IJ,N,IG)=ADS3(IJ)
            ZIDLNEXT(IJ,N,IG)=ZIDL3(IJ)
            ICEMASKNEXT(IJ,N,IG)=ICMSK3(IJ)

          ENDDO

!!!!          WRITE (IUNIT,ERR=6100,IOSTAT=IOS)
!!!!     &     CDTH, IG, (US3(I),I=IJS(IG),IJL(IG)),
!!!!     &     (DS3(I),I=IJS(IG),IJL(IG))
!!!!     &     , (ADS3(I),I=IJS(IG),IJL(IG))
!!!!     &     , (ZIDL3(I),I=IJS(IG),IJL(IG))
!!!!     &     , (ICMSK3(I),I=IJS(IG),IJL(IG))


        ENDDO

        DO IJ=IJS(IG),IJL(IG)
          US(IJ,IG) = US2(IJ,IG)
          DS(IJ,IG) = DS2(IJ,IG)
          ADS(IJ,IG) = ADS2(IJ,IG)
          ZIDL(IJ,IG) = ZIDL2(IJ,IG)
          ICMSK(IJ,IG) = ICMSK2(IJ,IG)
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

! ----------------------------------------------------------------------

!*    3. REWIND SCRATCH UNITS.
!        ---------------------

!!!!      DO IG=1,IGL
!!!!        IUNIT = IUSCR(IG)
!!!!        REWIND (UNIT=IUNIT)
!!!!      ENDDO

! ----------------------------------------------------------------------

!*    4. RE-ARRANGE THE BLOCKS.
!        ----------------------

!!!! the use of arrays instead of files will NOT work
!!!! if IDELPRO > IDELWI !!!!

!!!!      MP = MAX(IDELWI,IDELPRO)/IDELWO
!!!!      MSTEP = MAX(1,IDELPRO/IDELWO)

!*    4.1 LOOP OVER WIND FILES.
!         ---------------------


!!!!      DO M=1,MP,MSTEP

!*    4.1.1 LOOP OVER BLOCKS.
!           -----------------

!!!!        DO IG=1,IGL
!!!!          IUNIT = IUSCR(IG)

!*    4.1.1.1 LOOP OVER WIND FIELDS.
!             ----------------------

!!!!          DO MM=1,MSTEP
!!!!            READ (IUNIT,END=6200,IOSTAT=IOS)
!!!!     &       CDTH, IGG, (US(I,IG),I=IJS(IG),IJL(IG)),
!!!!     &       (DS(I,IG),I=IJS(IG),IJL(IG))
!!!!     &       , (ADS(I,IG),I=IJS(IG),IJL(IG))
!!!!     &       , (ZIDL(I,IG),I=IJS(IG),IJL(IG))
!!!!     &       , (ICMSK(I,IG),I=IJS(IG),IJL(IG))
!!!!            IF (CDTH.GE.CDTWIS)
!!!!     &       WRITE (IUVELO,ERR=6300,IOSTAT=IOS)
!!!!     &       CDTH, IGG, (US(I,IG),I=IJS(IG),IJL(IG)),
!!!!     &       (DS(I,IG),I=IJS(IG),IJL(IG))
!!!!     &       , (ADS(I,IG),I=IJS(IG),IJL(IG))
!!!!     &       , (ZIDL(I,IG),I=IJS(IG),IJL(IG))
!!!!     &       , (ICMSK(I,IG),I=IJS(IG),IJL(IG))
!!!!          ENDDO
!!!!        ENDDO
!!!!      ENDDO
!!!!      DO IG=1,IGL
!!!!        CLOSE (UNIT=IUSCR(IG) , STATUS='DELETE')
!!!!      ENDDO
!!!!      CLOSE (UNIT=IUVELO, STATUS='KEEP')

      IF (ITEST.GE.3) THEN
        WRITE(IU06,*) '       SUB. TIMIN: LAST WIND FIELD AT ',
     &   'CDTH = ', CDTH,' WRITTEN TO OUTPUT UNIT'
          CALL FLUSH (IU06)
      ENDIF

      RETURN

! ----------------------------------------------------------------------

!*    6. ERROR MESSAGES

 6100 CONTINUE
      WRITE (IU06,*) ' '
      WRITE (IU06,*) ' ********************************************'
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' *       FATAL ERROR IN SUB. TIMIN:         *'
      WRITE (IU06,*) ' *       ==========================         *'
      WRITE (IU06,*) ' * ERROR WHEN WRITTING ON SCRATCH UNIT      *'
      WRITE (IU06,*) ' * FOR BLOCK NUMBER        IG = ', IG
      WRITE (IU06,*) ' * NUMBER OF WIND FIELD IS  N = ', N
      WRITE (IU06,*) ' * DATE OF WIND FIELD IS CDTH = ', CDTH
      WRITE (IU06,*) ' * MESSAGE IS             IOS = ', IOS
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS         *'
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' ********************************************'
      CALL ABORT1
 6200 CONTINUE
      WRITE (IU06,*) ' '
      WRITE (IU06,*) ' ********************************************'
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' *       FATAL ERROR IN SUB. TIMIN:         *'
      WRITE (IU06,*) ' *       ==========================         *'
      WRITE (IU06,*) ' * ERROR WHEN READING FROM SCRATCH UNIT     *'
      WRITE (IU06,*) ' * FOR BLOCK NUMBER        IG = ', IG
      WRITE (IU06,*) ' * DATE OF WIND FIELD IS CDTH = ', CDTH
      WRITE (IU06,*) ' * MESSAGE IS             IOS = ', IOS
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS         *'
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' ********************************************'
      CALL ABORT1
 6300 CONTINUE
      WRITE (IU06,*) ' '
      WRITE (IU06,*) ' ********************************************'
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' *       FATAL ERROR IN SUB. TIMIN:         *'
      WRITE (IU06,*) ' *       ==========================         *'
      WRITE (IU06,*) ' * ERROR WHEN WRITTING ON WIND OUTPUT FILE  *'
      WRITE (IU06,*) ' * FOR BLOCK NUMBER       IG  = ', IG
      WRITE (IU06,*) ' * DATE OF WIND FIELD IS CDTH = ', CDTH
      WRITE (IU06,*) ' * MESSAGE IS            IOS  = ', IOS
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS         *'
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' ********************************************'
      CALL ABORT1

      END SUBROUTINE TIMIN
