      SUBROUTINE GETWND (US, DS, ADS, ZIDL, CICVR, CITH,
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

!       *CALL* *GETWND (US, DS, ADS, ZIDL, CICVR, CITH,
!                       CDTWIS, LWNDFILE, LCLOSEWND,
!                       LWCUR)*
!         *US*     - MAGNITUDE OF USTAR AT EACH POINT AND BLOCK.
!         *DS*     - DIRECTION OF USTAR AT EACH POINT AND BLOCK.
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
      USE YOWSTAT  , ONLY : NWAM_BLKS
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWIND  , ONLY : FIELDG   ,WSPMIN   ,IUNITW
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE GRIB_API_INTERFACE

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: IREAD
      INTEGER :: IG, MIJS, MIJL
      INTEGER :: JKGLO,KIJS,KIJL,NPROMA

      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION (NINF:NSUP,NBLO) :: US,DS,ADS,ZIDL,CICVR,CITH

      LOGICAL :: LWCUR
      LOGICAL :: LLNOTOPENED, LWNDFILE, LCLOSEWND
      LOGICAL :: LONLYONCE

      CHARACTER(LEN=14) :: CDTWIS, CDTWIR
      CHARACTER(LEN=24) :: FILNM

      DATA LONLYONCE /.TRUE./

      SAVE CDTWIR
      SAVE LLNOTOPENED

! ----------------------------------------------------------------------

!*    1. WIND DATA ARE READ
!        ------------------

      IF (LHOOK) CALL DR_HOOK('GETWND',0,ZHOOK_HANDLE)

 1000 CONTINUE

      IF(IUNITW.EQ.0) THEN
        LLNOTOPENED=.TRUE.
      ELSE
        LLNOTOPENED=.FALSE.
      ENDIF

!     GET FORCING FIELDS FROM INPUT FILES (if needed)
!     -----------------------------------
      IF(LWNDFILE) THEN
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
        NPROMA=(MIJL-MIJS+1)/NWAM_BLKS+1
        CALL GSTATS(1444,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=MIJS,MIJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,MIJL)
          CALL WAMWND (US(KIJS,IG), DS(KIJS,IG), ADS(KIJS,IG),
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
        NPROMA=(MIJL-MIJS+1)/NWAM_BLKS+1
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=MIJS,MIJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,MIJL)
          CALL MICEP(IPARAMCI,CICVR, CITH, KIJS,KIJL,IG)
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1444,1)

      ENDDO ! END LOOP ON IG

      IF (LHOOK) CALL DR_HOOK('GETWND',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE GETWND
