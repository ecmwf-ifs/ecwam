      SUBROUTINE TIMIN (CDTWIS, CDTWIE,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     &                  IUGB, IVGB, ILEN)

! ----------------------------------------------------------------------
!     MODIFIED  J. BIDLOT  FEBRARY 1996  MESSAGE PASSING


!**** *TIMIN* - STEERING MODULE IF TIME INTERPOLATION IS WANTED.

!     H. GUNTHER    ECMWF    MAY 1990     MODIFIED FOR SUB VERSION.
!     H. GUNTHER    ECMWF    DECEMBER 90  MODIFIED FOR CYCLE_4.

!*    PURPOSE.
!     --------

!       TIME INTERPOLATION: PROCESS WINDFIELDS.

!**   INTERFACE.
!     ----------

!       *CALL* *TIMIN (CDTWIS, CDTWIE,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
!    &                 IUGB, IVGB, ILEN)*
!          *CDTWIS* - DATE OF FIRST WIND FIELD.
!          *CDTWIE* - DATE OF LAST FIRST WIND FIELD.
!          *U10OLD* - WIND SPEED.
!          *THWOLD* - WIND DIRECTION (RADIANS).
!          *USOLD*  - FRICTION VELOCITY.
!          *TAUW*   - WAVE STRESS.
!          *Z0OLD*  - ROUGHNESS LENGTH IN.
!          *IUGB*   - INPUT GRIB WIND FIELD U COMPONENT.
!          *IVGB*   - INPUT GRIB WIND FIELD U COMPONENT.
!          *ILEN*   - LENGTH OG IUGB AND IVGB

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
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NBLO
      USE YOWPCONS , ONLY : PI       ,ZPI
      USE YOWSTAT  , ONLY : IDELPRO  ,IDELWI   ,IDELWO
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWUNIT  , ONLY : IUVELO   ,IUSCR
      USE YOWWIND  , ONLY : CDA      ,UWND     ,VWND

! ----------------------------------------------------------------------

      CHARACTER*12 CDTWIE, CDTH, CDT1, CDT2, CDTWIS, ZERO
      LOGICAL LWNDFILE, LCLOSEWND
      INTEGER IUGB(ILEN), IVGB(ILEN)
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD,USOLD,Z0OLD,TAUW
! ----------------------------------------------------------------------

      REAL,DIMENSION(NINF:NSUP,NBLO) :: US,DS
      REAL,DIMENSION(NINF:NSUP,NBLO) :: US2,DS2
      REAL,DIMENSION(NINF:NSUP) :: US3,DS3

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
        CALL GETWND (US, DS, NINF, NSUP, NSUP, CDT1,
     &   IUGB, IVGB, ILEN, LWNDFILE, LCLOSEWND)
        DO IG=1,IGL
          DO IJ=IJS(IG),IJL(IG)
            U10OLD(IJ,IG) = US(IJ,IG)
            THWOLD(IJ,IG) = DS(IJ,IG)
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
      CALL CREWFN (IUVELO, CDTWIE)
      IF (ITEST.GE.3) THEN
        WRITE(IU06,*) '       SUB. TIMIN: NEW WIND FILE '
        WRITE(IU06,*) '       UNIT IS IUVELO = ', IUVELO,
     &   ' WIND FILE DATE IS CDTWIE = ', CDTWIE
      ENDIF
      REWIND (UNIT=IUVELO)

!*    2. LOOP OVER INPUT WINDFIELDS.
!        ---------------------------

 2000 CONTINUE

!*    2.1 READ ONE WINDFIELD AND TRANSFORM ALL BLOCKS.
!         -------------------------------------------

      CDT2 = CDT1
      CALL INCDATE(CDT2,IDELWI)
      CALL GETWND (US2, DS2, NINF, NSUP, NSUP, CDT2,
     &             IUGB, IVGB, ILEN, LWNDFILE, LCLOSEWND)

!*    2.2 SAVE BLOCKED WIND FIELD ON SCRATCH UNITS.
!         -----------------------------------------

      DO IG=1,IGL
        IUNIT =IUSCR(IG)
        CDTH = CDT1
        DO N=1,NTS
          CALL INCDATE(CDTH,IDELWO)
          DO IJ=IJS(IG),IJL(IG)
            US3(IJ) = US(IJ,IG)+REAL(N)*DEL*(US2(IJ,IG)-US(IJ,IG))
            D = DS2(IJ,IG)-DS(IJ,IG)
            IF (ABS(D).GT.PI) D = D-ZPI*SIGN(1.,D)
            D = DS(IJ,IG)+REAL(N)*DEL*D
            DS3(IJ) = MOD(D+ZPI,ZPI)
          ENDDO
          WRITE (IUNIT,ERR=6100,IOSTAT=IOS)
     &     CDTH, IG, (US3(I),I=IJS(IG),IJL(IG)),
     &     (DS3(I),I=IJS(IG),IJL(IG))
        ENDDO
        DO IJ=IJS(IG),IJL(IG)
          US(IJ,IG) = US2(IJ,IG)
          DS(IJ,IG) = DS2(IJ,IG)
        ENDDO
      ENDDO
      CDT1 = CDT2
      IF (ITEST.GE.3) THEN
        WRITE(IU06,*) '       SUB. TIMIN: LAST WIND FIELD AT ',
     &   'CDTH = ', CDTH,' WRITTEN TO SCRATCH UNITS'
      ENDIF

!*    2.3 UPDATE WIND FIELD REQUEST TIME AND READ NEXT IF REQUESTED.
!         ----------------------------------------------------------

      CALL INCDATE (CDTH,IDELWI)
      IF (CDTH.LE.CDTWIE) GOTO 2000

! ----------------------------------------------------------------------

!*    3. REWIND SCRATCH UNITS.
!        ---------------------

      DO IG=1,IGL
        IUNIT = IUSCR(IG)
        REWIND (UNIT=IUNIT)
      ENDDO

! ----------------------------------------------------------------------

!*    4. RE-ARRANGE THE BLOCKS.
!        ----------------------

      MP = MAX(IDELWI,IDELPRO)/IDELWO
      MSTEP = MAX(1,IDELPRO/IDELWO)

!*    4.1 LOOP OVER WIND FILES.
!         ---------------------

      DO M=1,MP,MSTEP

!*    4.1.1 LOOP OVER BLOCKS.
!           -----------------

        DO IG=1,IGL
          IUNIT = IUSCR(IG)

!*    4.1.1.1 LOOP OVER WIND FIELDS.
!             ----------------------

          DO MM=1,MSTEP
            READ (IUNIT,END=6200,IOSTAT=IOS)
     &       CDTH, IGG, (US(I,IG),I=IJS(IG),IJL(IG)),
     &       (DS(I,IG),I=IJS(IG),IJL(IG))
            IF (CDTH.GE.CDTWIS)
     &       WRITE (IUVELO,ERR=6300,IOSTAT=IOS)
     &       CDTH, IGG, (US(I,IG),I=IJS(IG),IJL(IG)),
     &       (DS(I,IG),I=IJS(IG),IJL(IG))
          ENDDO
        ENDDO
      ENDDO
      DO IG=1,IGL
        CLOSE (UNIT=IUSCR(IG) , STATUS='DELETE')
      ENDDO
      CLOSE (UNIT=IUVELO, STATUS='KEEP')
      IF (ITEST.GE.3) THEN
        WRITE(IU06,*) '       SUB. TIMIN: LAST WIND FIELD AT ',
     &   'CDTH = ', CDTH,' WRITTEN TO OUTPUT UNIT'
      ENDIF

      IF(ALLOCATED(UWND)) DEALLOCATE(UWND)
      IF(ALLOCATED(VWND)) DEALLOCATE(VWND)

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
