      SUBROUTINE NOTIM  (CDTWIS, CDTWIE,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     &                   IUGB, IVGB, ILEN)

! ----------------------------------------------------------------------

!     MODIFIED  J. BIDLOT  FEBRARY 1996  MESSAGE PASSING

!**** *NOTIM* - STEERING MODULE IF NO TIME INTERPOLATION WANTED.

!*    PURPOSE.
!     --------

!       NOTIM NO TIME INTERPOLATION: PROCESS WINDFIELDS.

!**   INTERFACE.
!     ----------

!       *CALL* *NOTIM (CDTWIS, CDTWIE,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
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

      USE YOWPARAM , ONLY : NIBLO    ,NBLO
      USE YOWCOUP  , ONLY : LWCOU
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWSTAT  , ONLY : IDELPRO  ,IDELWO
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWUNIT  , ONLY : IUVELO   ,IUSCR
      USE YOWWIND  , ONLY : CDA      ,CDTNEXT  ,UWND     ,VWND     ,
     &            USNEXT   ,DSNEXT

! ----------------------------------------------------------------------
      CHARACTER*12 CDTWIE, CDTWIH, CDTWIS, CDT, ZERO
      LOGICAL LWNDFILE, LCLOSEWND
      INTEGER IUGB(ILEN), IVGB(ILEN)
! ----------------------------------------------------------------------

!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD,USOLD,Z0OLD,TAUW
! ----------------------------------------------------------------------

      REAL,DIMENSION(NINF:NSUP,NBLO) :: US,DS

!        *US*      REAL   OUTPUT WIND FIELD ARRAY (FRICTION VELOCITY).
!        *DS*      REAL   OUTPUT WIND FIELD ARRAY (DIRECTION).

! ----------------------------------------------------------------------

!*    1. INITIALISE WIND REQUEST DATE AND PROCESS FIRST WINDFIELD
!*       IF COLD START.
!        --------------------------------------------------------


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
        CALL GETWND (U10OLD, THWOLD, NINF, NSUP, NSUP, CDA,
     &               IUGB, IVGB, ILEN, LWNDFILE, LCLOSEWND)

        DO IG = 1,IGL
          CALL AIRSEA (U10OLD(IJS(IG),IG), TAUW(IJS(IG),IG),
     &     USOLD(IJS(IG),IG), Z0OLD(IJS(IG),IG),
     &     IJS(IG), IJL(IG), ILEV)
        ENDDO
        IF (ITEST.GE.3) THEN
          WRITE(IU06,'(''       SUB. NOTIM: FIRST WIND FIELD '',
     &     ''SAVED IN COMMON WIND'')')
        ENDIF
        IF (CDA.EQ.CDTWIE) THEN
          IF(ALLOCATED(UWND)) DEALLOCATE(UWND)
          IF(ALLOCATED(VWND)) DEALLOCATE(VWND)
          RETURN
        ENDIF
        CALL INCDATE (CDTWIH,IDELWO)
      ENDIF

! ----------------------------------------------------------------------

!*    2. CREATE WIND FILE NAME FROM END DATE.
!        ------------------------------------

!!!!      CALL CREWFN (IUVELO, CDTWIE)
!!!!      IF (ITEST.GE.3) THEN
!!!!         WRITE(IU06,*) '       SUB. NOTIM: NEW WIND FILE '
!!!!         WRITE(IU06,*) '       UNIT IS IUVELO = ', IUVELO,
!!!!     1                 ' WIND FILE DATE IS CDTWIE = ', CDTWIE
!!!!      ENDIF
!!!!      REWIND (UNIT=IUVELO)
      MP = 0

! ----------------------------------------------------------------------

!*    3. LOOP OVER OUTPUT WIND TIMES.
!        ----------------------------

 3000 CONTINUE

!*    3.1 READ ONE WIND FIELD AND TRANSFORM TO BLOCKS.
!         --------------------------------------------

      CALL GETWND (US, DS, NINF, NSUP, NSUP, CDTWIH,
     &             IUGB, IVGB, ILEN, LWNDFILE, LCLOSEWND)
      MP = MP + 1

!*    3.2 SAVE BLOCKED WIND FIELD.
!         ------------------------

      IF (IDELPRO.GT.IDELWO) THEN

!*        MORE THAN ONE FIELD PER PROPAGATION TIME STEP.
!*        SAVE ON SCRATCH UNITS:

        DO IG=1,IGL
          IUNIT =IUSCR(IG)
          WRITE (IUNIT,ERR=6100,IOSTAT=IOS)
     &     CDTWIH, IG, (US(I,IG),I=IJS(IG),IJL(IG)),
     &     (DS(I,IG),I=IJS(IG),IJL(IG))
        ENDDO
        IF (ITEST.GE.3) THEN
          WRITE(IU06,*) '       SUB. NOTIM: NEW WIND FILES AT ',
     &     'CDTWIH = ', CDTWIH,' WRITTEN TO SCRATCH UNIT'
          WRITE(IU06,'(''+'', A8)') IUNIT
        ENDIF

!*        UPDATE WIND FIELD REQUEST TIME.

        CALL INCDATE (CDTWIH,IDELWO)

!*        IF TIME LEFT BRANCH BACK TO 3.0 FOR NEXT FIELD.

        IF (CDTWIH.LE.CDTWIE) GOTO 3000
      ELSE

!*        WIND FIELD IS CONSTANT FOR ONE PROPAGATION TIME STEP:
!*        WRITE TO OUTPUT.
!         -----------------------------------------------------

        IF(.NOT.ALLOCATED(USNEXT)) ALLOCATE(USNEXT(NIBLO,NBLO))
        IF(.NOT.ALLOCATED(DSNEXT)) ALLOCATE(DSNEXT(NIBLO,NBLO))

        DO IG=1,IGL

!!!!            WRITE (IUVELO,ERR=6400,IOSTAT=IOS)
!!!     1            CDTWIH, IG, (US(I,IG),I=IJS(IG),IJL(IG)),
!!     2            (DS(I,IG),I=IJS(IG),IJL(IG))
!! 
!!save to array instead
          CDTNEXT=CDTWIH
          DO I=IJS(IG),IJL(IG)
            USNEXT(I,IG)=US(I,IG)
            DSNEXT(I,IG)=DS(I,IG)
          ENDDO

        ENDDO

!*        UPDATE WIND FIELD REQUEST TIME.
!         -------------------------------

        CALL INCDATE (CDTWIH,IDELWO)

!*        IF TIME LEFT BRANCH BACK TO 3.0 FOR NEXT FIELD.

!!!!!
!!!!         IF (CDTWIH.LE.CDTWIE) GOTO 3000
!!! we don't allow to loop anymore
        IF (CDTWIH.LE.CDTWIE) GOTO 6500

        IF (ITEST.GE.3) THEN
          WRITE(IU06,*) '       SUB. NOTIM: NEW WIND FILE GENERATED '
        ENDIF

!!!!!!!!!!        CLOSE (UNIT=IUVELO, STATUS='KEEP')

        IF(ALLOCATED(UWND)) DEALLOCATE(UWND)
        IF(ALLOCATED(VWND)) DEALLOCATE(VWND)

        RETURN

      ENDIF

! ----------------------------------------------------------------------

!*    4. RE-ARRANGE THE BLOCKS.
!        ----------------------

      DO IG=1,IGL
        REWIND (UNIT=IUSCR(IG))
      ENDDO

      MSTEP = IDELPRO/IDELWO
      DO M=1,MP,MSTEP
        DO IG=1,IGL
          IUNIT = IUSCR(IG)
          DO MM=1,MSTEP
            READ (IUNIT,ERR=6200,IOSTAT=IOS)
     &       CDT,IGG, (US(I,IG),I=IJS(IG),IJL(IG)),
     &       (DS(I,IG),I=IJS(IG),IJL(IG))

            WRITE (IUVELO,ERR=6300,IOSTAT=IOS)
     &       CDT, IGG, (US(I,IG),I=IJS(IG),IJL(IG)),
     &       (DS(I,IG),I=IJS(IG),IJL(IG))
          ENDDO
        ENDDO
      ENDDO
      REWIND (UNIT=IUVELO)
      DO IG=1,IGL
        CLOSE (UNIT=IUSCR(IG), STATUS='DELETE')
      ENDDO
      CLOSE (UNIT=IUVELO, STATUS='KEEP')
      IF (ITEST.GE.3) THEN
        WRITE(IU06,*) '       SUB. NOTIM: NEW WIND FILE ',
     &   'GENERATED AS COPY FROM  SCRATCH UNITS'
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
      WRITE (IU06,*) ' *       FATAL ERROR IN SUB. NOTIM:         *'
      WRITE (IU06,*) ' *       ==========================         *'
      WRITE (IU06,*) ' * ERROR WHEN WRITTING ON SCRATCH UNIT      *'
      WRITE (IU06,*) ' * FOR BLOCK NUMBER          IG = ', IG
      WRITE (IU06,*) ' * NUMBER OF WIND FIELD IS   MP = ', MP
      WRITE (IU06,*) ' * DATE OF WIND FIELD IS CDTWIH = ', CDTWIH
      WRITE (IU06,*) ' * MESSAGE IS               IOS = ', IOS
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS         *'
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' ********************************************'
      CALL ABORT1
 6200 CONTINUE
      WRITE (IU06,*) ' '
      WRITE (IU06,*) ' ********************************************'
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' *       FATAL ERROR IN SUB. NOTIM:         *'
      WRITE (IU06,*) ' *       ==========================         *'
      WRITE (IU06,*) ' * ERROR WHEN READING FROM SCRATCH UNIT     *'
      WRITE (IU06,*) ' * FOR BLOCK NUMBER       IG = ', IG
      WRITE (IU06,*) ' * DATE OF WIND FIELD IS CDT = ', CDT
      WRITE (IU06,*) ' * MESSAGE IS            IOS = ', IOS
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS         *'
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' ********************************************'
      CALL ABORT1
 6300 CONTINUE
      WRITE (IU06,*) ' '
      WRITE (IU06,*) ' ********************************************'
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' *       FATAL ERROR IN SUB. NOTIM:         *'
      WRITE (IU06,*) ' *       ==========================         *'
      WRITE (IU06,*) ' * ERROR WHEN WRITTING ON WIND OUTPUT FILE  *'
      WRITE (IU06,*) ' * AFTER READING FROM SCRATCH UNIT          *'
      WRITE (IU06,*) ' * FOR BLOCK NUMBER       IG = ', IG
      WRITE (IU06,*) ' * DATE OF WIND FIELD IS CDT = ', CDT
      WRITE (IU06,*) ' * MESSAGE IS            IOS = ', IOS
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS         *'
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' ********************************************'
      CALL ABORT1
      WRITE (IU06,*) ' '
      WRITE (IU06,*) ' ********************************************'
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' *       FATAL ERROR IN SUB. NOTIM:         *'
      WRITE (IU06,*) ' *       ==========================         *'
      WRITE (IU06,*) ' * ERROR WHEN WRITTING ON WIND OUTPUT FILE  *'
      WRITE (IU06,*) ' * SCRATCH UNITS ARE NOT USED               *'
      WRITE (IU06,*) ' * FOR BLOCK NUMBER          IG = ', IG
      WRITE (IU06,*) ' * DATE OF WIND FIELD IS CDTWIH = ',CDTWIH
      WRITE (IU06,*) ' * MESSAGE IS               IOS = ', IOS
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS         *'
      WRITE (IU06,*) ' *                                          *'
      WRITE (IU06,*) ' ********************************************'
      CALL ABORT1
 6500 CONTINUE
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

      END SUBROUTINE NOTIM
