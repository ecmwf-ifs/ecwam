      SUBROUTINE GETWND (U10, THW, ADS, ZIDL, CICVR, CITH,
     &                   CDTWIS, LWNDFILE, LCLOSEWND, IREAD,
     &                   LWCUR)

! ----------------------------------------------------------------------

!**** *GETWND* - ROUTINE TO READ AND PROCESS ONE WINDFIELD.

!*    PURPOSE.
!     --------

!        READ A WINDFIELD FROM THE WINDFILE (SEARCH FOR IT)
!        AND CALCULATES THE WIND VELOCITY  AND DIRECTION
!        FOR ALL WAM BLOCKS.
!        ALSO INPUT THE OTHER FORCING FIELD.

!**   INTERFACE.
!     ----------

!       *CALL* *GETWND (U10, THW, ADS, ZIDL, CICVR, CITH,
!                       CDTWIS, LWNDFILE, LCLOSEWND,
!                       LWCUR)*
!         *U10*    - MAGNITUDE OF 10m WIND AT EACH POINT AND BLOCK.
!         *THW*    - DIRECTION OF 10m WIND AT EACH POINT AND BLOCK.
!         *ADS*    - AIR DENSITY AT EACH POINT AND BLOCK.
!         *ZIDL*   - Zi/L  AT EACH POINT AND BLOCK.
!         *CICVR*  - SEA ICE COVER.
!         *CITH*   - SEA ICE THICKNESS.
!         *CDTWIS* - DATE OF WIND FIELD TO BE LOOKED FOR.
!         *LWNDFILE - FLAG USED TO DETERMINE WHETHER WINDS ARE READ FROM
!                     FILE OR ARE AVAILABLE IN ARRAY FIELDG (SEE *IFSTOWAM).
!         *LCLOSEWND* IF TRUE THE INPUT FILE WILL BE CLOSED AND
!                     THE UNIT RESET
!         *IREAD*  - PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!                    (IF NEEDED)
!         *LWCUR*  -  LOGICAL INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS
!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.
!       *READWIND*   - READING WINDS.
!       *WAMWND*    - CALCULATE WIND IN WAM POINTS.

!     REFERENCE.
!     ----------

!       NONE.


!    MODIFIED BY:
!    ------------
!    B. HANSEN    ECMWF 1997
!                 RESTRUCTURE CALL TO READWIND.
!
!    S. ABDALLA   ECMWF OCTOBER 2001
!                 MODIFICATION THE CALL TO READWIND WAMWND; AND 
!                 INCLUSION OF AIR DENSITY AND Zi/L.
!    J. BIDLOT    ECMWF NOVEMEBR 2003
!                 INTRODUCE openMP
!    J. BIDLOT    ECMWF AUGUST 2006
!                 SEA ICE FRACTION IF COUPLED.
!    J. BIDLOT    ECMWF AUGUST 2008
!                 READWIND WAS SPLIT BETWEEN READWIND AND IFSTOWAM:
!                 IN COUPLED RUNS, READWIND IS NOT CALLED 
! ----------------------------------------------------------------------
      USE YOWCOUP  , ONLY : LWCOU
      USE YOWICE   , ONLY : IPARAMCI
      USE YOWMPP   , ONLY : IRANK    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NBLO     ,NIBLO    ,
     &             CLDOMAIN,LWDINTS
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWSTAT  , ONLY : NPROMA_WAM 
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWIND  , ONLY : WSPMIN   ,IUNITW
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWUNPOOL, ONLY : LLUNSTR
      USE UNWAM, ONLY : USE_DIRECT_WIND_FILE
      USE UNSTRUCT_WIND, ONLY : SET_WIND_UNSTRUCTURED
      USE GRIB_API_INTERFACE

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IREAD

      REAL, DIMENSION (NINF:NSUP,NBLO), INTENT(INOUT) :: 
     &                 U10,THW,ADS,ZIDL,CICVR,CITH

      CHARACTER(LEN=14), INTENT(IN) :: CDTWIS

      LOGICAL, INTENT(IN) :: LWNDFILE, LCLOSEWND, LWCUR


      INTEGER :: IG, MIJS, MIJL
      INTEGER :: JKGLO,KIJS,KIJL,NPROMA
      INTEGER :: I, J

      REAL :: ZHOOK_HANDLE

      CHARACTER(LEN=14), INTENT(IN) :: CDTWIS

      LOGICAL, INTENT(IN) :: LWNDFILE, LCLOSEWND, LWCUR

      LOGICAL :: LLNOTOPENED
      LOGICAL :: LONLYONCE
      LOGICAL :: IsAssigned

      CHARACTER(LEN=14) :: CDTWIR
      CHARACTER(LEN=24) :: FILNM

      DATA LONLYONCE /.TRUE./

      SAVE CDTWIR
      SAVE LLNOTOPENED

! ----------------------------------------------------------------------

!*    1. WIND DATA ARE READ
!        ------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('GETWND',0,ZHOOK_HANDLE)
#endif

 1000 CONTINUE

      IF(IUNITW.EQ.0) THEN
        LLNOTOPENED=.TRUE.
      ELSE
        LLNOTOPENED=.FALSE.
      ENDIF

      NPROMA=NPROMA_WAM

!     GET FORCING FIELDS FROM INPUT FILES (if needed)
!     -----------------------------------
      IsAssigned=.FALSE.
      IF (LLUNSTR .and. USE_DIRECT_WIND_FILE) THEN
        IsAssigned=.TRUE.
        CALL SET_WIND_UNSTRUCTURED
      END IF


      IF(LWNDFILE .and. (IsAssigned .eqv. .FALSE.)) THEN
        CALL READWIND (CDTWIR, FILNM, LLNOTOPENED, IREAD)

!       CHECK WIND FIELD DATE

        IF (CDTWIR.LT.CDTWIS) THEN
!         DATE OF INPUT FIELD IS BEFORE REQUESTED DATE
!         TRY AGAIN
          IF(LWNDFILE) THEN
            IF (ITEST.GT.1) THEN
              WRITE(IU06,*) ' SUB. GETWND - BEFORE REQUESTED DATE '
              WRITE(IU06,*) ' CDTWIR= ',CDTWIR
              WRITE(IU06,*) ' CDTWIS= ',CDTWIS
              WRITE(IU06,*) ' SUB. GETWND - CALLING READWIND AGAIN'
              CALL FLUSH(IU06)
            ENDIF
            GOTO 1000
          ELSE
            WRITE (IU06,*) ' ****************************************'
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *      FATAL ERROR SUB. GETWND          *'
            WRITE (IU06,*) ' *      =======================          *'
            WRITE (IU06,*) ' * WIND DATE IS EARLIER THAN EXPECTED    *'
            WRITE (IU06,*) ' * DECODED DATE IS  CDTWIR = ', CDTWIR
            WRITE (IU06,*) ' * DATE EXPECTED IS CDTWIS = ', CDTWIS
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' ****************************************'
            CALL ABORT1
          ENDIF
        ELSEIF (CDTWIR.GT.CDTWIS) THEN

!         DATE OF INPUT FIELD IS LATER THAN REQUESTED DATE
          WRITE (IU06,*) ' ****************************************'
          WRITE (IU06,*) ' *                                      *'
          WRITE (IU06,*) ' *      FATAL ERROR SUB. GETWND         *'
          WRITE (IU06,*) ' *      =======================         *'
          WRITE (IU06,*) ' * WIND DATE IS LATER THAN EXPECTED     *'
          IF(LWNDFILE) THEN
            WRITE (IU06,*) ' * DATE READ IS    CDTWIR = ', CDTWIR
          ELSE
            WRITE (IU06,*) ' * DECODED DATE IS CDTWIR = ', CDTWIR
          ENDIF
          WRITE (IU06,*) ' * DATE EXPECTED IS CDTWIS = ', CDTWIS
          WRITE (IU06,*) ' *                                      *'
          WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS     *'
          WRITE (IU06,*) ' *                                      *'
          WRITE (IU06,*) ' ****************************************'
          CALL ABORT1
        ENDIF

        IF (LCLOSEWND .AND. LWNDFILE .AND. 
     &     .NOT.(CLDOMAIN.EQ.'s' .OR. LWDINTS) ) THEN
          IF(IRANK.EQ.IREAD) THEN
            CALL IGRIB_CLOSE_FILE(IUNITW)
            LLNOTOPENED = .TRUE.
            IUNITW=0
            IF (ITEST.GT.1)
     &        WRITE(IU06,*) ' SUB. GETWND - CLOSE ',FILNM
          ENDIF
        ENDIF

      ENDIF
! ----------------------------------------------------------------------

!*    3. INTERPOLATE AND BLOCK WINDFIELD
!        -------------------------------

      DO IG=1,IGL
        MIJS=IJS(IG)
        MIJL=IJL(IG)
! Mod for OPENMP
        CALL GSTATS(1444,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=MIJS,MIJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,MIJL)
          CALL WAMWND (U10(KIJS,IG), THW(KIJS,IG), ADS(KIJS,IG),
     &                 ZIDL(KIJS,IG), CITH(KIJS,IG), KIJS, KIJL, IG,
     &                 LWCUR)
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1444,1)

        IF(LONLYONCE) THEN
          WRITE (IU06,*) ' '
          WRITE (IU06,*) ' SUB. GETWND : '
          WRITE (IU06,*) ' '
          WRITE (IU06,*) ' WIND SPEEDS LOWER THAN ',WSPMIN, ' M/S'
          WRITE (IU06,*) ' WERE RESET TO  ',WSPMIN, ' M/S'
          WRITE (IU06,*) ' '
          CALL FLUSH(IU06)
          LONLYONCE=.FALSE.
        ENDIF

!       USE THE SEA ICE FRACTION TO DEFINE THE SEA ICE BOUNDARY

        CALL GSTATS(1444,0)
! Mod for OPENMP
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=MIJS,MIJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,MIJL)
          CALL MICEP(IPARAMCI,CICVR, CITH, KIJS,KIJL,IG)
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1444,1)

      ENDDO ! END LOOP ON IG

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('GETWND',1,ZHOOK_HANDLE)
#endif

      RETURN
      END SUBROUTINE GETWND
