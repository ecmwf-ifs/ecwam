      SUBROUTINE MSTART (IU12, IU14, IU15, IOPTI, FETCH, FRMAX,
     &                   FL1,U10OLD,THWOLD,USOLD,Z0OLD)
! ----------------------------------------------------------------------

!**** *MSTART* - MAKES START FIELDS FOR WAMODEL.

!      H. GUNTHER    ECMWF    MAY 1990
!      H. GUNTHER    ECMWF    DECEMBER 90  MODIFIED FOR CYCLE_4.
!      J. BIDLOT     ECMWF    FEBRUARY 96  MESSAGE PASSING
!      J. BIDLOT     ECMWF    MARCH 97     REMOVE ALL OUTPUT TO IU12
!                                          and IU15 

!*    PURPOSE.
!     --------

!       TO GENERATE WAMODEL START FIELDS.

!**   INTERFACE.
!     ----------

!   *CALL* *MSTART (IU12, IU14, IU15, IOPTI, FETCH, FRMAX,
!    1              FL1,U10OLD,THWOLD,USOLD,Z0OLD)*
!      *IU12*   INTEGER    OUTPUT UNIT BLOCKS OF SPECTRA.
!      *IU14*   INTEGER    OUTPUT UNIT SECOND LAT OF BLOCKS.
!      *IU15*   INTEGER    OUTPUT UNIT LAST WINDFIELDS.
!      *IOPTI*  INTEGER    START FIELD OPTION
!                          = 0 FROM PARAMTERS.
!                          = 1 FROM WINDS CALM ENERGY=0.
!                          = 2 FROM WINDS CALM FROM PARAMETERS.
!      *FETCH*  REAL       FETCH IN METERS.
!      *FRMAX*  REAL       MAXIMUM PEAK FREQUENCY IN HERTZ.
!      *FL1*      REAL      2-D SPECTRUM FOR EACH GRID POINT 
!      *U10NEW*    NEW WIND SPEED IN M/S.
!      *U10OLD*    INTERMEDIATE STORAGE OF MODULUS OF WIND
!                  VELOCITY.
!      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!      *THWOLD*    INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
!                  WIND VELOCITY.
!      *USNEW*     NEW FRICTION VELOCITY IN M/S.
!      *USOLD*     INTERMEDIATE STORAGE OF MODULUS OF FRICTION
!                  VELOCITY.
!      *Z0OLD*     INTERMEDIATE STORAGE OF ROUGHNESS LENGTH IN
!                  M.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       *ABORT*     - TERMINATES PROCESSING.
!       *PEAK*      - COMPUTE PARAMETERS FROM WIND FOR A BLOCK.
!       *SPECTRA*   - COMPUTES SPECTRA OF A BLOCK.
!       *SPLITBL*   - STORING OVERLAPPING LATITUDES OF BLOCKS.

!    REFERENCE.
!    ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLO    ,NBLO
      USE YOWPCONS , ONLY : DEG
      USE YOWCOUT  , ONLY : NGOUT    ,IGAR     ,IJAR
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL2     ,IJL
      USE YOWJONS  , ONLY : FP       ,ALPHJ    ,THES     ,FM       ,
     &            ALFA     ,GAMMA    ,SA       ,SB
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL

! ----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD,USOLD,Z0OLD
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL1 

! ----------------------------------------------------------------------

      CHARACTER*12 ZERO
      PARAMETER(ZERO=' ')

! ----------------------------------------------------------------------

!     0. ALLOCATE ARRAYS
!        ---------------

      ALLOCATE(FP(NIBLO))
      ALLOCATE(ALPHJ(NIBLO))
      ALLOCATE(THES(NIBLO))

!*    1. DEFINE SPECTRUM FOR LAND POINTS AND WRITE OUTPUT.
!        -------------------------------------------------

      DO M=1,NFRE
        DO K=1,NANG
          FL1(NINF-1,K,M) = 0.
        ENDDO
      ENDDO

      WRITE (IU06,'(1X,/,1X,''  PARAMETER AT OUTPUT SIDES:'')')
      WRITE (IU06,'(1X,''  NGOU    IG    IJ     U10    UDIR'',
     &            ''      FP   ALPHAJONS   GAMMA      SA      SB'')')

! ----------------------------------------------------------------------

!*    2. LOOP FOR BLOCKS.
!        ----------------

      DO IG=1,IGL

!*    2.1 COMPUTE PEAK FREQUENCIES AND ALPHAJONS PARAMETERS.
!         ----------------------------------------------


!*    2.1.1 INITIAL VALUES DUE TO OPTION.
!           -----------------------------

        IF (IOPTI.EQ.1) THEN
          DO IJ = IJS(IG), IJL(IG)
            FP(IJ) = 0.
            ALPHJ(IJ) = 0.
            THES(IJ) = THWOLD(IJ,IG)
          ENDDO
        ELSE
          DO IJ = IJS(IG), IJL(IG)
            FP(IJ) = FM
            ALPHJ(IJ) = ALFA
            IF (U10OLD(IJ,IG) .GT. 0.1E-08) THEN
              THES(IJ) = THWOLD(IJ,IG)
            ELSE
              THES(IJ) = THES(IJ)
            ENDIF
          ENDDO
        ENDIF

!*    2.1.2 PEAK FREQUENCY AND ALPHAJONS FROM FETCH LAW.
!           --------------------------------------------

        IF (IOPTI.NE.0) THEN
          CALL PEAK (IJS(IG), IJL(IG), IG, FETCH, FRMAX,U10OLD)
        ENDIF
        IF (ITEST.GT.1) THEN
          IF (IG.LE.ITESTB) WRITE (IU06,*) '    SUB. PEAK DONE'
        ENDIF

!*    2.1.3 PRINT PARAMETERS AT OUTPUT POINTS.
!           ----------------------------------

        DO NGOU = 1, NGOUT
          IF (IG.EQ.IGAR(NGOU)) THEN
            IJ = IJAR(NGOU)
            WRITE (IU06,'(1X,3I6,F8.2,F8.2,5F8.4)')  NGOU, IG, IJ,
     &       U10OLD(IJ,IG), THWOLD(IJ,IG)*DEG, FP(IJ), ALPHJ(IJ),
     &       GAMMA, SA, SB
          ENDIF
        ENDDO

!*    2.2 COMPUTE SPECTRA FROM PARAMETERS.
!         --------------------------------

        CALL SPECTRA (IJS(IG), IJL(IG), IG, FL1)
        IF (ITEST.GT.1) THEN
          IF (IG.LE.ITESTB) WRITE (IU06,*) '    SUB. SPECTRA DONE'
        ENDIF

!*    2.3 OUTPUT OF OVERLAPPING LATITUDES.
!         --------------------------------

        IF (IG.GT.1) CALL SPLITBL (FL1, IJS(IG), IJL2(IG), IU14)
        IF (ITEST.GT.1) THEN
          IF (IG.LE.ITESTB) WRITE (IU06,*) '    SUB. SPLITBL DONE'
        ENDIF


!*    BRANCHING BACK TO 2. FOR NEXT BLOCK.

      ENDDO

! ----------------------------------------------------------------------

!*    3. PREPARE OUTPUT OF WIND.
!        ----------------------

      CDTPRO  = ZERO
      CDATEWO = ZERO
      CDAWIFL = ZERO
      CDATEFL = ZERO

      RETURN
      END SUBROUTINE MSTART
