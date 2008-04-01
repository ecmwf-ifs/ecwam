      SUBROUTINE GETWND (US, DS, ADS, ZIDL, ICMSK,
     &                   KINF, KSUP, CDTWIS,
     &                   NFIELDS, NGPTOTG, NC, NR, FCRANGE,
     &                   ISEC1, ISEC2, FIELDS, MASK_IN,
     &                   LWNDFILE, LCLOSEWND)

! ----------------------------------------------------------------------

!**** *GETWND* - ROUTINE TO READ AND PROCESS ONE WINDFIELD.

!*    PURPOSE.
!     --------

!        READ A WINDFIELD FROM THE WINDFILE (SEARCH FOR IT)
!        AND CALCULATES THE WIND VELOCITY  AND DIRECTION
!        FOR ALL WAM BLOCKS.

!**   INTERFACE.
!     ----------

!       *CALL* *GETWND (US, DS, ADS, ZIDL, KINF, KSUP, CDTWIS,
!                       NFIELDS, NGPTOTG, NC, NR, FCRANGE,
!                       ISEC1, ISEC2, FIELDS,
!                       LWNDFILE, LCLOSEWND)*
!         *US*     - MAGNITUDE OF USTAR IN EACH POINT AND BLOCK.
!         *DS*     - DIRECTION OF USTAR IN EACH POINT AND BLOCK.
!         *ADS*    - AIR DENSITY IN EACH POINT AND BLOCK.
!         *ZIDL*   - Zi/L  IN EACH POINT AND BLOCK.
!         *KINF*   - LOWER BOUND ON FIRST DIEMENSION OF FIELD
!         *KSUP*   - UPPER BOUND ON FIRST DIEMENSION OF FIELD
!         *CDTWIS* - DATE OF WIND FIELD TO BE LOOKED FOR.
!         *NFIELDS*- NUMBER OF FIELDS HOLDING ATMOSPHERIC DATA
!         *NGPTOTG*- NUMBER OF ATMOSPHERIC GRID POINTS
!         *NC*     - NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!         *NR*     - NUMBER OF ATM. ROWS OF LATITUDES
!         *FCRANGE*- FORECAST RANGE IN SECONDS
!         *ISEC1*  - GRIB ISEC1 (TO EXTRACT DATE/TIME OF ATM. FIELDS)
!         *ISEC2*  - GRIB ISEC2 (TO EXTRACT ATM. GRID CHARACTERISTICS)
!         *FIELDS* - ATMOSPHERIC FIELDS (U10, V10, AIR DENSITY, Zi/L)
!         *LWNDFILE - FLAG USED TO DETERMINE WHETHER WINDS ARE READ FROM
!                     FILE OR PASSED AS ARGUMENTS IN IUGB AND IVGB
!         *LCLOSEWND* IF TRUE THE INPUT FILE WILL BE CLOSED AND
!                     THE UNIT RESET

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
! ----------------------------------------------------------------------

      USE YOWICE   , ONLY : IPARAMCI ,CI
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NBLO     ,NIBLO    ,
     &             CLDOMAIN,LWDINTS
      USE YOWGRIBHD, ONLY : NKSEK1 
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWSTAT  , ONLY : NWAM_BLKS
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWIND  , ONLY : UWND     ,VWND     ,AIRD     ,PZIDL     ,
     &               WSPMIN,IUNITW

! ----------------------------------------------------------------------
      INTEGER NFIELDS, NGPTOTG, NC, NR
      INTEGER ISEC1(NKSEK1), ISEC2(22+NR)
      REAL FIELDS(NGPTOTG,NFIELDS)
      INTEGER MASK_IN(NGPTOTG)
! ----------------------------------------------------------------------

      REAL, DIMENSION (KINF:KSUP,NBLO) :: US, DS, ADS, ZIDL
      INTEGER, DIMENSION (KINF:KSUP,NBLO) :: ICMSK 

! Mod for OPENMP
      INTEGER JKGLO,KIJS,KIJL,NPROMA
! End Mod for OPENMP
      LOGICAL LLNOTOPENED, LWNDFILE, LCLOSEWND
      LOGICAL LONLYONCE

      CHARACTER*14 CDTWIS, CDTWIR
      CHARACTER*24 :: FILNM

      DATA LONLYONCE /.TRUE./

      SAVE CDTWIR
      SAVE LLNOTOPENED

! ----------------------------------------------------------------------

!*    0. GLOBAL FIELD FOR THE INPUTS
!        ---------------------------

      ALLOCATE(UWND(NGX,NGY))
      ALLOCATE(VWND(NGX,NGY))
      ALLOCATE(AIRD(NGX,NGY))
      ALLOCATE(PZIDL(NGX,NGY))
      ALLOCATE(CI(NGX,NGY))

!*    1. WIND DATA ARE READ
!        ------------------

 1000 CONTINUE

      IF(IUNITW.EQ.0) THEN
        LLNOTOPENED=.TRUE.
      ELSE
        LLNOTOPENED=.FALSE.
      ENDIF

      CALL READWIND (CDTWIR, NFIELDS, NGPTOTG, NC, NR, FCRANGE,
     &               ISEC1, ISEC2, FIELDS, MASK_IN,
     &               LWNDFILE, LCLOSEWND,
     &               KINF, KSUP,
     &               FILNM, LLNOTOPENED )

!*    2. CHECK WIND FIELD DATE
!        ----------------------

      IF (CDTWIR.LT.CDTWIS) THEN

!     2.1 DATE OF INPUT FIELD IS BEFORE REQUESTED DATE
!         TRY AGAIN
!         --------------------------------------------
        IF(LWNDFILE) THEN
          GOTO 1000
        ELSE
          WRITE (IU06,*) ' ******************************************'
          WRITE (IU06,*) ' *                                         *'
          WRITE (IU06,*) ' *      FATAL ERROR SUB. GETWND            *'
          WRITE (IU06,*) ' *      =======================            *'
          WRITE (IU06,*) ' * WIND DATE IS EARLIER THAN EXPECTED      *'
          WRITE (IU06,*) ' * DECODED DATE IS  CDTWIR = ', CDTWIR
          WRITE (IU06,*) ' * DATE EXPECTED IS CDTWIS = ', CDTWIS
          WRITE (IU06,*) ' *                                         *'
          WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS        *'
          WRITE (IU06,*) ' *                                         *'
          WRITE (IU06,*) ' ******************************************'
          CALL ABORT1
        ENDIF
      ELSEIF (CDTWIR.GT.CDTWIS) THEN

!     2.2 DATE OF INPUT FIELD IS LATER THAN REQUESTED DATE
!         ------------------------------------------------
        WRITE (IU06,*) ' ******************************************'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' *      FATAL ERROR SUB. GETWND           *'
        WRITE (IU06,*) ' *      =======================           *'
        WRITE (IU06,*) ' * WIND DATE IS LATER THAN EXPECTED       *'
        IF(LWNDFILE) THEN
          WRITE (IU06,*) ' * DATE READ IS     CDTWIR = ', CDTWIR
        ELSE
          WRITE (IU06,*) ' * DECODED DATE IS  CDTWIR = ', CDTWIR
        ENDIF
        WRITE (IU06,*) ' * DATE EXPECTED IS CDTWIS = ', CDTWIS
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS       *'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' ******************************************'
        CALL ABORT1
      ENDIF


      IF (LCLOSEWND .AND. LWNDFILE .AND. 
     &   .NOT.(CLDOMAIN.EQ.'s' .OR. LWDINTS) ) THEN
        CALL PBCLOSE(IUNITW,KRET)
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '****************************************'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO PBCLOSE     *'
          WRITE (IU06,*) '*  IN GETWND                           *'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '****************************************'
          CALL ABORT1
        ENDIF
        LLNOTOPENED = .TRUE.
        IUNITW=0
        IF (ITEST.GT.1)
     &    WRITE(IU06,*) ' SUB. GETWND - PBCLOSE ON ',FILNM
      ENDIF

! ----------------------------------------------------------------------

!*    3. INTERPOLATE AND BLOCK WINDFIELD
!        -------------------------------

      DO IG=1,IGL
        IF(KINF.EQ.1.AND.KSUP.EQ.NIBLO.AND.LMESSPASS) THEN
          MIJS=KINF
          MIJL=KSUP
        ELSE
          MIJS=IJS(IG)
          MIJL=IJL(IG)
        ENDIF

! Mod for OPENMP
        NPROMA=(MIJL-MIJS+1)/NWAM_BLKS+1
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=MIJS,MIJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,MIJL)
          CALL WAMWND (US, DS, ADS, ZIDL, KINF, KSUP, KIJS, KIJL, IG)
        ENDDO
!$OMP   END PARALLEL DO

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
! Mod for OPENMP
        NPROMA=(MIJL-MIJS+1)/NWAM_BLKS+1
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=MIJS,MIJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,MIJL)
          CALL MICEP (CI, IPARAMCI, ICMSK, KINF, KSUP, KIJS, KIJL, IG)
        ENDDO
!$OMP   END PARALLEL DO

      ENDDO ! END LOOP ON IG

!*    4. DEALLOCATE GRID ARRAYS FOR INPUT FORCING FIELDS
!        -----------------------------------------------

      DEALLOCATE(UWND)
      DEALLOCATE(VWND)
      DEALLOCATE(AIRD)
      DEALLOCATE(PZIDL)
      DEALLOCATE(CI)

      RETURN
      END SUBROUTINE GETWND
