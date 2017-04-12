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
!                                       IN PLACE OF COMMON BLOCKS 
!     P.A.E.M. JANSSEN  ECMWF JULY 2010 : INTRODUCE 2ND ORDER INTERACTION
!                                         COEFFICIENTS
!
!*    PURPOSE.
!     --------

!       TO ARRANGE A GRID FOR THE WAM WAVE MODEL AND COMPUTE
!       ALL FIXED MODEL PARAMETERS WHICH ARE STORED IN DIFFERENT
!       MODULE.

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
!                  (SEE SUB TOPOAR AND MUBUF).
!       *IU02*   - LOGICAL UNIT FOR INPUT OF CURRENTS.
!                  (SEE SUB READCUR).
!       *IU03*   - LOGICAL UNIT FOR INPUT OF COARSE GRID
!                  BOUNDARY ORGANISATION (MODULE CBOUND).
!                  IF THIS IS A FINE GRID PREPROC.
!                  FORMATED IF IFORM = 2 OTHERWISE UNFORMATED.
!                  (SEE SUB MBOUNF).
!       *IU06*   - LOGICAL UNIT FOR PRINTER OUTPUT UNIT
!       *IU07*   - LOGICAL UNIT FOR OUTPUT OF GRID ORGANISATION
!                  AND COMPUTED CONSTANTS. (UNFORMATED)
!                  (SEE SUB OUTCOM).
!       *IU08*   - LOGICAL UNITS FOR OUTPUT OF MODULE UBUF.
!                  (UNFORMATED) (SEE SUB MUBUF).
!       *IU09*   - LOGICAL UNIT FOR UNFORMATED OUTPUT OF COARSE
!                  GRID BOUNDARY ORGANISATION (MODULE CBOUND),
!                  IF THIS IS A COARSE GRID PREPROC.
!                  (SEE SUB MBOUNC).
!       *IU10*   - LOGICAl UNIT FOR UNFORMATED OUTPUT OF FINE
!                  GRID BOUNDARY ORGANISATION (MODULE CBOUND).
!                  IF THIS IS A FINE GRID PREPROC.
!                  (SEE SUB MBOUNF).
!       *IU17*   - SAME AS IU07 BUT FORMATED.
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

!       *ADJUST*    - CORRECTS LONGITUDE INPUT.
!       *AKI*       - COMPUTES WAVE NUMBER.
!       *CHECK*     - CHECKS CONSISTENCY OF BLOCK OVERLAPS.
!       *FINDB*     - FIND BLOCK AND GRID POINT NUMBERS.
!       *JAFU*      - ANGULAR INDEX OF NON LINEAR INTERACTION
!       *INIWCST*   - SETS GLOBAL CONSTANTS.
!       *LOCINT*    - INTERPOLATE TO BLOCKS.
!       *MBLOCK*    - PREPARES ONE BLOCK
!       *MBOUNC*    - MAKE COARSE GRID BOUNDARY.
!       *MBOUNF*    - MAKE FINE   GRID BOUNDARY.
!       *MBOXB*     - MAKE BOX FOR FINE GRID IN COARSE GRID.
!       *MCOUT*     - PREPARES OUTPUT MODULE COUT
!       *MFREDIR*   - COMPUTES FREQUENCY/DIRECTION MODULE FREDIR
!       *MGRID*     - ARRANGES GRID FOR MODEL.
!       *MINTF*     - MAKE INTERPOLATION TABLES FOR BOUNDARY INPUT.
!       *MTABS*     - COMPUTES TABLES USED FOR SHALLOW WATER
!       *MUBUF*     - COMPUTES MODULE UBUF.
!       *NLWEIGT*   - COMPUTES NON LINEAR WEIGHTS IN MODULE INDNL
!       *OUTCOM*    - OUTPUT OF COMPUTED MODULES
!       *OUTPP*     - OUTPUT OF GRID.
!       *PACKI*     - PACKS AN INTEGER ARRAY.
!       *PACKR*     - PACKS A REAL ARRAY.
!       *READCUR*   - READ A CURRENT FIELD.
!       *SECONDHH*  - COMPUTES TABLES FOR SECOND HARMONICS EVALUATION
!       *STRESS*    - STRESS TABLE.
!       *TAUHF*     - HIGH FREQUENCY STRESS TABLE.
!       *TOPOAR*    - PREPARE TOPOGRAPHY FOR GRID.
!       *UIPREP*    - READS USER INPUT

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWALTAS , ONLY : EGRCRV   ,AGRCRV   ,BGRCRV   ,AFCRV    ,
     &            BFCRV    ,ESH      ,ASH      ,BSH      ,ASWKM    ,
     &            BSWKM
      USE YOWPARAM , ONLY : NIBLO    ,NBLO     ,NIBLD    ,NBLD     ,
     &            NIBLC    ,NBLC
      USE YOWCPBO  , ONLY : IBOUNC   ,NBOUNC
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,XKAPPA   ,
     &            XNLEV    ,TAUWSHELTER, TAILFACTOR, TAILFACTOR_PM
      USE YOWFPBO  , ONLY : IBOUNF   ,NBOUNF
      USE YOWFRED  , ONLY : FR       ,DFIM     ,GOM      ,C        ,
     &            TH       ,COSTH    ,SINTH
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,SINPH    ,COSPH    ,
     &            NLONRGG  ,IGL
      USE YOWMAP   , ONLY : NX       ,NY       ,IPER     ,IRGG     ,
     &            AMOWEP   ,AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLA   ,
     &            XDELLO   ,ZDELLO   ,LAQUA
      USE YOWPHYS  , ONLY : ALPHAPMAX
      USE YOWSTAT  , ONLY : IPHYS
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWPCONS , ONLY : OLDPI    ,CIRC     ,RAD
      USE YOWUNIT  , ONLY : NPROPAGS ,IU07     ,IU08
      USE YOWUNPOOL ,ONLY : LLUNSTR  ,LPREPROC
      USE UNWAM

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IU01, IU02, IU03, IU09, IU10, IU17, IU19, IU20 
      INTEGER :: IG
      INTEGER :: K, IX, ICL, IFORM, ML, KL, LNAME, IINPC, LFILE
      INTEGER :: I_GET_UNIT

      REAL :: PRPLRADI
      REAL :: OLDRAD
      REAL :: XLAT, XLATD, XLATMAX,PLONS, COSPHMIN
      REAL, ALLOCATABLE :: BATHY(:,:)

      CHARACTER(LEN=1) :: C1 
      CHARACTER(LEN=80) :: FILENAME

      LOGICAL :: LLEXIST
      LOGICAL :: LLGRID

! ----------------------------------------------------------------------

      LMESSPASS=.FALSE.
      PRPLRADI=1.0

      CALL INIWCST(PRPLRADI)

!*    1. INITIALISATION OF INPUT/OUTPUT UNITS.
!        -------------------------------------

      IU02  = 2
      IU06  = 6

! ----------------------------------------------------------------------

!*    2. USER INPUT AND LINEPRINTER PROTOCOL.
!        ------------------------------------

      CALL UIPREP (IFORM, ML, KL, LLGRID)

! ----------------------------------------------------------------------

!*    2.1 ALLOCATE NECESSARY ARRAYS
!     -----------------------------

      ALLOCATE(BATHY(NX, NY))

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
!     NLONRGG IS ALLOCATED IN UIPREP
      ALLOCATE(ZDELLO(NY))

      IU01 = I_GET_UNIT(IU06, 'wam_topo', 'r', 'f', 0)

      FILENAME='wam_topo'
      LLEXIST=.FALSE.
      LNAME = LEN_TRIM(FILENAME)
      INQUIRE(FILE=FILENAME(1:LNAME),EXIST=LLEXIST)
      IF(.NOT. LLEXIST) THEN
        WRITE (IU06,*) '*************************************'
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
        WRITE (IU06,*) '*  IN PREPROC:                     *'
        WRITE (IU06,*) '*  COULD NOT FIND FILE ',FILENAME
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*************************************'
        CALL ABORT1
      ENDIF
      IU01 = I_GET_UNIT(IU06, FILENAME(1:LNAME), 'r', 'f', 0)

      IF (IFORM.NE.2) THEN
        IU07 = I_GET_UNIT(IU06, 'wam_grid_tables', 'w', 'u', 0)
      ELSE
        IU17 = I_GET_UNIT(IU06, 'wam_grid_tables_form', 'w', 'f', 0)
      ENDIF

      DO ICL=0,NPROPAGS
        WRITE(C1,'(I1)') ICL
        FILENAME='wam_subgrid_'//C1
        LFILE=0
        IF (FILENAME.NE. ' ') LFILE=LEN_TRIM(FILENAME)
        IU08(ICL) = I_GET_UNIT(IU06,FILENAME(1:LFILE) , 'w', 'u', 0)
      ENDDO

      IF (IBOUNC.EQ.1) THEN
!       Information of the nested grid(s) that will be produce by a coarse grid run
        IF (IFORM.NE.2) THEN
          IU09=I_GET_UNIT(IU06,'wam_nested_grids_info','w', 'u', 0)
        ELSE
          IU19=I_GET_UNIT(IU06,'wam_nested_grids_info_form','w','f',0)
        ENDIF
      ENDIF

      IF (IBOUNF.EQ.1) THEN
!       Information of the nested grid(s) that were produced by a coarse grid run
        IF (IFORM.NE.2) THEN
          IU03=I_GET_UNIT(IU06,'wam_nested_grids_from_coarse_info',
     &                    'r', 'u', 0)
        ELSE
          IU03=I_GET_UNIT(IU06,'wam_nested_grids_info_from_coarse_form',
     &                    'r','f',0)
        ENDIF

!       Information about the boundary points that will be needed for a fine
!       grid run.
        IF (IFORM.NE.2) THEN
          IU10=I_GET_UNIT(IU06,'wam_boundary_grid_info', 'w', 'u', 0)
        ELSE
          IU20=I_GET_UNIT(IU06,'wam_boundary_grid_info_form','w','f',0)
        ENDIF
      ENDIF

! ----------------------------------------------------------------------

!*    3. INITIALISE TOTAL NUMBER OF BLOCKS,
!*       AND GRID INCREMENTS IN RADIANS AND METRES.
!        ------------------------------------------

      IGL=0
      DELPHI = XDELLA*CIRC/360.
      OLDRAD=OLDPI/180.
      DO K=1,NY
        XLAT       = (AMOSOP + REAL(K-1)*XDELLA)*RAD
        XLATD      = (AMOSOP + REAL(K-1)*XDELLA)
        SINPH(K)   = SIN(XLAT)
        COSPH(K)   = COS(XLAT)
        IF(.NOT.LLGRID) THEN
          IF (IRGG.EQ.1) THEN
            IF(XDELLA.EQ.0.25.AND.AMOWEP.EQ.-98..AND.AMOSOP.EQ.9..AND.
     &         AMOEAP.EQ.42..AND.AMONOP.EQ.81.) THEN
!         the old value for pi has to be taken in order to reproduce
!         exactly the irregular grid of the operational LAW model 0.25
              NLONRGG(K)=NINT(NX*COS((AMOSOP+REAL(K-1)*XDELLA)*OLDRAD))
            ELSE
#ifdef _CRAYFTN
              NLONRGG(K) = MAX(NINT(NX*COSD(XLATD)),2)
#else
              NLONRGG(K) = MAX(NINT(NX*COS(RAD*XLATD)),2)
#endif
            ENDIF
            IF(MOD(NLONRGG(K),2).EQ.1) NLONRGG(K) = NLONRGG(K)+1
          ELSE
            NLONRGG(K) = NX
          ENDIF      
        ENDIF      

        PLONS=(AMOEAP-AMOWEP) + IPER*XDELLO
        IF(NX.EQ.1.AND.NY.EQ.1) THEN
          NLONRGG(K) = NX
          ZDELLO(K)  = XDELLO
          DELLAM(K)  = ZDELLO(K)*CIRC/360.
          EXIT
        ENDIF
        IF(IPER.EQ.1) THEN
          ZDELLO(K)  = 360./REAL(NLONRGG(K))
        ELSE
          ZDELLO(K)  = PLONS/REAL(NLONRGG(K)-1)
        ENDIF
        DELLAM(K)  = ZDELLO(K)*CIRC/360.
      ENDDO

!     IF THE POLES ARE INCLUDED, ARTIFICIALLY REMOVE THE SINGULARITY
      XLATMAX=87.5
      COSPHMIN=COS(XLATMAX*RAD)
      DO K=1,NY
        IF(COSPH(K).LE.COSPHMIN) THEN
          COSPH(K)=COS(XLATMAX*RAD)
          SINPH(K)=SIN(XLATMAX*RAD)
        ENDIF
      ENDDO 

! ----------------------------------------------------------------------

!*    4. COMPUTE GRID INDEPENDENT MODULE.
!        -------------------------------


!*    4.1 MODULE FREDIR (FREQUENCY/DIRECTION CONST).
!         ------------------------------------------

      CALL MFREDIR (ML, KL)

!*    4.2 MODULE INDNL (WEIGHT OF NON-LINEAR INTERACTION).
!         ------------------------------------------------

      CALL NLWEIGT

!*    4.3 MODULE SHALLOW (SHALLOW WATER TABLES).
!         --------------------------------------

      CALL MTABS (ML, KL)

!*    4.4 MODULE COUPLE.
!         --------------
      IF (IPHYS.EQ.0) THEN
!       ECMWF PHYSICS:
        IF(ML.GT.30) THEN
          ALPHA   = 0.0060
        ELSE
          ALPHA   = 0.0075
        ENDIF
!!!        ALPHA   = 0.00525
        BETAMAX = 1.20
        ZALP    = 0.008
        ALPHAPMAX = 0.03
        TAUWSHELTER=0.0
        TAILFACTOR=2.5
!!!        TAILFACTOR_PM=3.4
        TAILFACTOR_PM=0.   ! i.e. not used

!!!     EMPIRICAL CONSTANCE FOR  SPECTRAL UPDATE FOLLOWING DATA ASSIMILATION
        EGRCRV = 1108.0
        AGRCRV = 0.06E+6
        BGRCRV = 9.70
        AFCRV = 4.0E-4
        BFCRV = -3.0
        ESH = 1711.
        ASH = 8.0E-4 
        BSH = 0.96 
        ASWKM=0.0981
        BSWKM=0.425

      ELSE IF (IPHYS.EQ.1) THEN
!       METEO FRANCE PHYSICS WITH HIGH FREQUENCy TAIL FROM ECMWF PHYSICS
!        ALPHA   = 0.0065
!        BETAMAX = 1.52
!        ZALP    = 0.0060
!        ALPHAPMAX = 0.03
!        TAUWSHELTER=0.5
!        TAILFACTOR=2.5
!        TAILFACTOR_PM=0.0
!       TEST 471 adapted to ECWAM
        ALPHA   = 0.006
        BETAMAX = 1.38
        ZALP    = 0.006
        ALPHAPMAX = 0.04
        TAUWSHELTER=0.3
        TAILFACTOR=2.5
        TAILFACTOR_PM=3.5

!!!     EMPIRICAL CONSTANCE FOR  SPECTRAL UPDATE FOLLOWING DATA ASSIMILATION
        EGRCRV = 1000.0
        AGRCRV = 0.113E+6
        BGRCRV = 4.97 
        AFCRV = 4.5E-4
        BFCRV = -3.0
        ESH = 1711.
        ASH = 8.0E-4 
        BSH = 0.96 
        ASWKM=0.0981
        BSWKM=0.425

      ELSE
        WRITE (IU06,*) '*************************************'
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*  ERROR IN PREPROC                 *'
        WRITE (IU06,*) '*  UKNOWN PHYSICS SELECTION :       *'
        WRITE (IU06,*) '*  IPHYS =' , IPHYS
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*************************************'
        CALL ABORT1
      ENDIF

      XKAPPA  = 0.40

!*    UP TO JPLEVC LEVELS CAN BE SPECIFIED.

      XNLEV(1)= 10.0

!*    4.4 MODULE TABLE (STRESS TABLES).
!         -----------------------------

      CALL STRESS(IU06,ITEST)


!*    4.5 HIGHER HARMONICS TABLES
!         -----------------------
!
!         FIRST ROUTINE IS FOR DEEP WATER ONLY AND IS APPLIED FOR 
!         ALTIMETRIC CORRECTIONS
!         SECOND ROUTINE IS MORE GENERAL, HENCE APPLIES TO DEEP AND 
!         SHALLOW WATER
!     

      CALL SECONDHH

      CALL SECONDHH_GEN

! ----------------------------------------------------------------------

!*    5. GENERATE OUTPUT GRID INFORMATION.
!        ---------------------------------


!*    5.1 READ IN TOPOGRAPHY AND ARRANGE ON REQUESTED MODEL AREA OF
!*        REQUESTED RESOLUTION.
!         ---------------------------------------------------------

      IF(.NOT.LAQUA) THEN
        CALL TOPOAR (IU01, BATHY)
        IF (ITEST.GT.0) THEN
          WRITE (IU06,*) ' SUB TOPOAR DONE'
          CALL FLUSH(IU06)
        ENDIF
      ELSE
!       AQUA PLANET SET TO DEEP EVERYWHERE
!       EXCEPT AT THE POLES THAT ARE EXCLUDED AS LAND.
        BATHY(:,:)=999.
        DO IX=1,NX
          BATHY(IX,1)=-999.
          BATHY(IX,NY)=-999.
        ENDDO
      ENDIF

!*    5.2 COMPUTATION OF BLOCKS.
!         ----------------------

      CALL MGRID (BATHY)
      IF (ITEST.GT.0) THEN 
        WRITE (IU06,*) ' SUB MGRID DONE'
        CALL FLUSH(IU06)
      ENDIF

      IF (LLUNSTR .AND. LPREPROC) THEN

        CALL INIT_UNWAM

      END IF ! LLUNSTR

      NIBLD=0
      NBLD=0
      NIBLC=0
      NBLC=0

!*    5.3 COMPUTE OUTPUT POINT INDICES (MODULE YOWCOUT).
!         ----------------------------------------------

      CALL MCOUT
      IF (ITEST.GT.0) WRITE (IU06,*) ' SUB MCOUT DONE'

! ----------------------------------------------------------------------

!*    6. COMPUTE NEST INFORMATION.
!        -------------------------


!*    6.1 COMPUTE FINE GRID NEST INFORMATION (MODULE YOWFPBO).
!         ---------------------------------------------------
!AR: We just omit the nest informations in this module since it is anyway broken ...
      IF (.NOT. LLUNSTR) THEN

        IF (IBOUNF.EQ.1) THEN
          CALL MBOUNF (IU03, IU10, IU20, IFORM, IINPC)
          IF (ITEST.GT.0) WRITE (IU06,*) ' SUB MBOUNF DONE'
        ELSE
          IINPC  = 0
          NBOUNF = 0
        ENDIF

!*    6.2 COMPUTE COARSE GRID NEST INFORMATION (MODULE YOWCPBO).
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

!     THE READING OF THE INPUT CURRENTS IN NOW DONE IN INITMDL 

! ----------------------------------------------------------------------

!*    8. GENERATE AND WRITE MODULE UBUF.
!        -------------------------------

        DO IG=1,IGL
          CALL MUBUF (IU01, BATHY, IG, ML, IU08, NPROPAGS)
          IF (ITEST.GT.0) THEN
            IF (IG.LE.ITESTB) WRITE (IU06,*)
     &       ' SUB MUBUF DONE FOR BLOCK ',IG
          ENDIF
        ENDDO
!
      END IF ! .NOT. LLUNSTR


! ----------------------------------------------------------------------

!*    9. OUTPUT OF MODULES.
!        ------------------

      CALL OUTCOM (IU07, IU17, IFORM)
      IF (ITEST.GT.0) WRITE (IU06,*) ' SUB OUTCOM DONE'

! ----------------------------------------------------------------------

!*    10. CONSISTENCY CHECK OF COMPUTED BLOCK PARAMETERS AND
!*        OUTPUT OF NECESSARY DIMENSIONS.
!         --------------------------------------------------

      CALL CHECK (ML, KL, IINPC)
  
      END PROGRAM
