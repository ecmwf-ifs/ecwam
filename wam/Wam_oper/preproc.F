      PROGRAM preproc 

! ----------------------------------------------------------------------

!**** *PROGRAM PREPROC* - PREPARE DATA (BUT NOT WINDS) FOR INPUT
!                         TO WAM WAVE MODELS.

!     SUSANNE HASSELMANN  MPI     JUNE 1986.

!     ANNEGRET SPEIDEL    MPI  OCTOBER 1988. MODFIED FOR CYCLE_2.

!     K. HUBBERT          POL     JUNE 1989  DEPTH AND CURRENT
!                                            REFRACTION.

!     H. GUNTHER   ECMWF/GKSS    APRIL 1990  LAND POINTS ARE REMOVED
!                                            FROM BLOCKS AND THE CODE
!                                            HAS BEEN RESTRUCTURED.

!     R. PORTZ     MPI         JANUARY 1991  NESTED GRID OPTION.

!     H. GUNTHER   ECMWF/GKSS    APRIL 1991  CYCLE_4 MODIFICATIONS.
!                                            MULTI-PART REMOVED.
!                                            NEW SOURCE FUNCTIONS.
!                                            LOG. DEPTH TABLE.
!     J. BIDLOT    ECMWF         SEPTEMBER 1996 REDUCED GRID.

!     J. BIDLOT    ECMWF OCTOBER 1998 : INTRODUCE USE OF MODULES
!                                       IN PLACE COMMON BLOCKS 

!*    PURPOSE.
!     --------

!       TO ARRANGE A GRID FOR THE WAM WAVE MODEL AND COMPUTE
!       ALL FIXED MODEL PARAMETERS WHICH ARE STORED IN DIFFERENT
!       COMMON BLOCKS.

!     METHOD.
!     -------

!       A REPRESENTATIVE TOPOGRAPHIC DATA SET ON LAT-LONG
!       COORDINATES CONTAINING THE MODEL SQUARE BOX REGION IS
!       READ IN.THE MODEL REGION IS EXTRACTED AND INTERPOLATED
!       ONTO GIVEN LAT-LONG GRID INCREMENTS (SEE SUB TOPOAR).
!       THE PROGRAM CHECKS FOR A PERIODIC LATITUDE GRID. IF THE
!       GRID IS NOT PERIODIC A CLOSED BASIN IS ASSUMED.
!       THE PROGRAM DOES NOT DISTINGUISH BETWEEN DEEP AND SHALLOW
!       WATER.

!       -BLOCK STRUCTURE :
!        GRID POINTS ARE COLLECTED INTO A 1-DIMENSIONAL ARRAY,
!        BLOCKS OF MAXIMALLY NIBLO ELEMENTS,  GRID POINTS
!        (ONLY SEAPOINTS) ARE COUNTED ALONG LINES OF LATITUDES
!        FROM WEST TO EAST WORKING FROM SOUTH TO NORTH.
!        BLOCKS OVERLAP OVER TWO LATITUDE LINES,TO COMPUTE NORTH
!        -SOUTH ADVECTION TERMS.

!       -NESTED GRIDS: THE GRID GENERATED CAN BE A
!         - COARSE GRID WHICH MEANS OUTPUT OF SPECTRA
!                       FOR A FOLLOW UP FINE GRID RUN.
!         - FINE   GRID WHICH MEANS INPUT OF SPECTRA
!                       FROM  AN EARLIER COARSE GRID RUN.
!         - COARSE AND FINE GRID

!       - REFRACTION: CONTROLLED BY THE REFRACTION OPTION
!         A CURRENT FIELD IS READ, INTERPOLATED TO THE MODEL
!         GRID AND STORED IN THE GRID OUTPUT FILE.

!       - PARAMETERS FOR ARRAY DIMENSIONS: THE PRORAM CHECKS
!         ALL DIMENSIONS INTERNALLY. ONLY THE BLOCK LENGTH
!         (NIBLO) IS USED FOR THE SET UP OF THE GRID, ALL
!         THE OTHER PARAMETERS HAVE TO BE LARGE ENOUGH TO
!         GET A SUCCESFULL RUN OF THE JOB. AT THE END OF
!         THE OUTPUT PROTOCOLL A LIST IS PRINTED FOR THE
!         OPTIMAL SETTINGS OF THE DIMENSION.

!**   INTERFACE.
!     ----------

!       *PROGRAM* *PREPROC*

!       *IU01*   - LOGICAL UNIT FOR INPUT OF TOPOGRAPHIC DATA.
!                  (SEE SUB TOPOAR).
!       *IU02*   - LOGICAL UNIT FOR INPUT OF CURRENTS.
!                  (SEE SUB READCUR).
!       *IU03*   - LOGICAL UNIT FOR INPUT OF COARSE GRID
!                  BOUNDARY ORGANISATION (COMMON CBOUND).
!                  IF THIS IS A FINE GRID PREPROC.
!                  FORMATED IF IFORM = 2 OTHERWISE UNFORMATED.
!                  (SEE SUB MBOUNF).
!       *IU05*   - LOGICAL UNIT FOR USER INPUT. (SEE SUB UIPREP).
!       *IU06*   - LOGICAL UNIT FOR PRINTER OUTPUT UNIT
!       *IU07*   - LOGICAL UNIT FOR OUTPUT OF GRID ORGANISATION
!                  AND COMPUTED CONSTANTS. (UNFORMATED)
!                  (SEE SUB OUTCOM).
!       *IU08*   - LOGICAL UNIT FOR OUTPUT OF COMMON UBUF.
!                  (UNFORMATED) (SEE SUB OUTUBUF).
!       *IU09*   - LOGICAL UNIT FOR UNFORMATED OUTPUT OF COARSE
!                  GRID BOUNDARY ORGANISATION (COMMON CBOUND),
!                  IF THIS IS A COARSE GRID PREPROC.
!                  (SEE SUB MBOUNC).
!       *IU10*   - LOGICAl UNIT FOR UNFORMATED OUTPUT OF FINE
!                  GRID BOUNDARY ORGANISATION (COMMON CBOUND).
!                  IF THIS IS A FINE GRID PREPROC.
!                  (SEE SUB MBOUNF).
!       *IU17*   - SAME AS IU07 BUT FORMATED.
!       *IU18*   - SAME AS IU08 BUT FORMATED.
!       *IU19*   - SAME AS IU09 BUT FORMATED.
!       *IU20*   - SAME AS IU10 BUT FORMATED.

!       ALL UNITS ARE DEFINE IN SECTION 1. OF THIS PROGRAM.

!       MODULES YOWPARAM, YOWCOUPL, YOWCURR, YOWFRED, YOWINDNL, YOWGRID,
!       YOWMAP, YOWCOUT, YOWTABL, AND YOWSHAL ARE WRITTEN TO UNIT
!       IU07 AND/OR IU17.
!       ALL FREQUENCY AND DIRECTION DEPENDENT ARRAYS ARE WRITTEN FROM
!       1 TO THE USED NUMBER OF FREQUENCIES AND THE USED NUMBER OF
!       DIRECTIONS.
!       OTHER ARRAYS ARE WRITTEN ACCORDING TO THEIR DIMENSIONS.

!     EXTERNALS.
!     ----------

!       *ABORT*     - TERMINATES PROCESSING.
!       *ADJUST*    - CORRECTS LONGITUDE INPUT.
!       *AKI*       - COMPUTES WAVE NUMBER.
!       *CHECK*     - CHECKS CONSISTENCY OF BLOCK OVERLAPS.
!       *FINDB*     - FIND BLOCK AND GRID POINT NUMBERS.
!       *JAFU*      - ANGULAR INDEX OF NON LINEAR INTERACTION
!       *LOCINT*    - INTERPOLATE TO BLOCKS.
!       *MBLOCK*    - PREPARES ONE BLOCK
!       *MBOUNC*    - MAKE COARSE GRID BOUNDARY.
!       *MBOUNF*    - MAKE FINE   GRID BOUNDARY.
!       *MBOXB*     - MAKE BOX FOR FINE GRID IN COARSE GRID.
!       *MCOUT*     - PREPARES OUTPUT COMMON COUT
!       *MFREDIR*   - COMPUTES FREQUENCY/DIRECTION COMMON FREDIR
!       *MGRID*     - ARRANGES GRID FOR MODEL.
!       *MINTF*     - MAKE INTERPOLATION TABLES FOR BOUNDARY INPUT.
!       *MTABS*     - COMPUTES TABLES USED FOR SHALLOW WATER
!       *MUBUF*     - COMPUTES COMMON UBUF.
!       *NLWEIGT*   - COMPUTES NON LINEAR WEIGHTS IN COMMON INDNL
!       *OUTCOM*    - OUTPUT OF COMPUTED COMMONS
!       *OUTPP*     - OUTPUT OF GRID.
!       *OUTUBUF*   - OUTPUT OF BLOCK COMMON UBUF
!       *PACKI*     - PACKS AN INTEGER ARRAY.
!       *PACKR*     - PACKS A REAL ARRAY.
!       *READCUR*   - READ A CURRENT FIELD.
!       *STRESS*    - STRESS TABLE.
!       *TAUHF*     - HIGH FREQUENCY STRESS TABLE.
!       *TOPOAR*    - PREPARE TOPOGRAPHY FOR GRID.
!       *UIPREP*    - READS USER INPUT
!       *WAMCUR*    - PREPARE CURRENT FIELD.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NIBLO    ,NBLO     ,NIBLD    ,NBLD     ,
     &            NIBLC    ,NBLC
      USE YOWCPBO  , ONLY : IBOUNC   ,NBOUNC
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,XKAPPA   ,
     &            XNLEV
      USE YOWFPBO  , ONLY : IBOUNF   ,NBOUNF
      USE YOWFRED  , ONLY : FR       ,DFIM     ,GOM      ,C        ,
     &            TH       ,COSTH    ,SINTH
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,SINPH    ,COSPH    ,
     &            NLONRGG  ,IGL
      USE YOWMAP   , ONLY : NX       ,NY       ,IPER     ,IRGG     ,
     &            AMOWEP   ,AMOSOP   ,AMOEAP   ,XDELLA   ,XDELLO   ,
     &            ZDELLO
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWPCONS , ONLY : CIRC     ,RAD

! ----------------------------------------------------------------------

      INTEGER, ALLOCATABLE :: IA1(:,:)

!          *IA1*   INTEGER  TOPOGRAPHIC DATA IN MODEL GRID.

! ----------------------------------------------------------------------

      LMESSPASS=.FALSE.

!*    1. INITIALISATION OF INPUT/OUTPUT UNITS.
!        -------------------------------------

      IU01  = 1
      IU02  = 2
      IU03  = 3
      IU05  = 5
      IU06  = 6
      IU07  = 7
      IU08  = 8
      IU09  = 9
      IU10  = 10
      IU17  = 17
      IU18  = 18
      IU19  = 19
      IU20  = 20

! ----------------------------------------------------------------------

!*    2. USER INPUT AND LINEPRINTER PROTOCOL.
!        ------------------------------------

      CALL UIPREP (IU05, IFORM, IREFRA, ML, KL)

! ----------------------------------------------------------------------

!*    2.1 ALLOCATE NECESSARY ARRAYS
!     -----------------------------

      ALLOCATE(IA1(NX, NY))

!     FR(KL) IS ALLOCATED IN UIPREP 
      ALLOCATE(DFIM(ML)) 
      ALLOCATE(GOM(ML)) 
      ALLOCATE(C(ML)) 
      ALLOCATE(TH(KL)) 
      ALLOCATE(COSTH(KL)) 
      ALLOCATE(SINTH(KL)) 

      ALLOCATE(DELLAM(NY)) 
      ALLOCATE(SINPH(NY)) 
      ALLOCATE(COSPH(NY)) 
      ALLOCATE(NLONRGG(NY)) 
      ALLOCATE(ZDELLO(NY))

! ----------------------------------------------------------------------

!*    3. INITIALISE TOTAL NUMBER OF BLOCKS,
!*       AND GRID INCREMENTS IN RADIANS AND METRES.
!        ------------------------------------------

      IGL=0
      DELPHI = XDELLA*CIRC/360.
      DO K=1,NY
        XLAT       = (AMOSOP + REAL(K-1)*XDELLA)*RAD
        SINPH(K)   = SIN(XLAT)
        COSPH(K)   = COS(XLAT)
        IF (IRGG.EQ.1) THEN
          NLONRGG(K) = NINT(NX*COSPH(K))
          IF(MOD(NLONRGG(K),2).eq.1) NLONRGG(K) = NLONRGG(K)+1
        ELSE
          NLONRGG(K) = NX
        ENDIF      

        PLONS=(AMOEAP-AMOWEP) + IPER*XDELLO
        IF(IPER.EQ.1) THEN
          ZDELLO(K)  = PLONS/REAL(NLONRGG(K))
        ELSE
          ZDELLO(K)  = PLONS/REAL(NLONRGG(K)-1)
        ENDIF
        DELLAM(K)  = ZDELLO(K)*CIRC/360.
      ENDDO

! ----------------------------------------------------------------------

!*    4. COMPUTE GRID INDEPENDENT COMMON BLOCKS.
!        ---------------------------------------


!*    4.1 COMMON FREDIR (FREQUENCY/DIRECTION CONST).
!         ------------------------------------------

      CALL MFREDIR (ML, KL)

!*    4.2 COMMON INDNL (WEIGHT OF NON-LINEAR INTERACTION).
!         ------------------------------------------------

      CALL NLWEIGT (ML, KL)

!*    4.3 COMMON SHALLOW (SHALLOW WATER TABLES).
!         --------------------------------------

      CALL MTABS (ML, KL)

!*    4.4 COMMON COUPLE.
!         --------------

      BETAMAX = 1.20
      ZALP    = 0.0110
      ALPHA   = 0.0090
      XKAPPA  = 0.40

!*    UP TO JPLEVC LEVELS CAN BE SPECIFIED.

      XNLEV(1)= 10.0

!*    4.4 COMMON TABLE (STRESS TABLES).
!         -----------------------------

      CALL STRESS(IU06,ITEST)
      CALL TAUHF (FR(ML))

! ----------------------------------------------------------------------

!*    5. GENERATE OUTPUT GRID INFORMATION.
!        ---------------------------------


!*    5.1 READ IN TOPOGRAPHY AND ARRANGE ON REQUESTED MODEL AREA OF
!*        REQUESTED RESOLUTION.
!         ---------------------------------------------------------

      CALL TOPOAR (IU01, IA1)

!*    5.2 COMPUTATION OF BLOCKS.
!         ----------------------

      CALL MGRID (IA1)
      IF (ITEST.GT.0) WRITE (IU06,*) ' SUB MGRID DONE'

      IF(IREFRA.EQ.2) THEN
        NIBLD=NIBLO
        NBLD=NBLO
        NIBLC=NIBLO
        NBLC=NBLO
      ELSEIF(IREFRA.EQ.1) THEN
        NIBLD=NIBLO
        NBLD=NBLO
        NIBLC=1
        NBLC=1
      ELSE
        NIBLD=1
        NBLD=1
        NIBLC=1
        NBLC=1
      ENDIF

!*    5.3 COMPUTE OUTPUT POINT INDICES (MODULE YOWCOUT).
!         -------------------------------------------

      CALL MCOUT
      IF (ITEST.GT.0) WRITE (IU06,*) ' SUB MCOUT DONE'

! ----------------------------------------------------------------------

!*    6. COMPUTE NEST INFORMATION.
!        -------------------------


!*    6.1 COMPUTE FINE GRID NEST INFORMATION (MODULE YOWFBOU).
!         ---------------------------------------------------

      IF (IBOUNF.EQ.1) THEN
        CALL MBOUNF (IU03, IU10, IU20, IFORM, IINPC)
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB MBOUNF DONE'
      ELSE
        IINPC  = 0
        NBOUNF = 0
      ENDIF

!*    6.2 COMPUTE COARSE GRID NEST INFORMATION (MODULE YOWCBOU).
!         -----------------------------------------------------

      IF (IBOUNC.EQ.1) THEN
        CALL MBOUNC (IU09, IU19, IFORM)
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB MBOUNC DONE'
      ELSE
        NBOUNC = 0
      ENDIF

! ----------------------------------------------------------------------

!*    7. GENERATE CURRENT FIELD.
!        -----------------------

      IF (IREFRA.EQ.2) THEN
        CALL WAMCUR (IU02)
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB WAMCUR DONE'
      ENDIF

! ----------------------------------------------------------------------

!*    8. GENERATE AND WRITE COMMON UBUF.
!        -------------------------------

      DO IG=1,IGL
        CALL MUBUF (IA1, IG, IU08, IU18, IFORM)
        IF (ITEST.GT.0) THEN
          IF (IG.LE.ITESTB) WRITE (IU06,*)
     &     ' SUB MUBUF DONE FOR BLOCK ',IG
        ENDIF
      ENDDO

! ----------------------------------------------------------------------

!*    9. OUTPUT OF COMMON BLOCKS.
!        ------------------------

      CALL OUTCOM (IU07, IU17, IFORM, IREFRA, ML, KL)
      IF (ITEST.GT.0) WRITE (IU06,*) ' SUB OUTCOM DONE'

! ----------------------------------------------------------------------

!*    10. CONSISTENCY CHECK OF COMPUTED BLOCK PARAMETERS AND
!*        OUTPUT OF NECESSARY DIMENSIONS.
!         --------------------------------------------------

      CALL CHECK (IREFRA, ML, KL, IINPC)

      STOP
      END




