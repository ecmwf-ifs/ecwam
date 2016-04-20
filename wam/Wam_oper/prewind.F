      SUBROUTINE PREWIND (U10OLD, THWOLD, USOLD, TAUW, Z0OLD,
     &                    ROAIRO, ZIDLOLD,
     &                    CICOVER, CITHICK, CIWA,
     &                    LLINIT, LLALLOC_FIELDG_ONLY,
     &                    IREAD,
     &                    NFIELDS, NGPTOTG, NC, NR,
     &                    FIELDS, LWCUR, MASK_IN)

! ----------------------------------------------------------------------

!**** *PREWIND* - PREPARES WIND DATA FOR WAVE MODEL.

!     P.GROENWOUD     DELFT HYDRAULICS LABORATORY  OKTOBER 1986

!     E. BAUER        MPI       FEB 1987   VERSION FOR CDC 205 HAMBURG.

!     S. HASSELMANN   MPI       MAY 1987   COMBINED CDC 205 AND CRAY
!                                          VERSIONS.
!     W. BRUEGGEMANN  MPI      AUGUST 1988   SIMPLIFIED PROGRAM.

!     L. ZAMBRESKY    ECMWF    JUNE 1988   MODIFIED EXTENSIVELY FOR
!                                          COUPLING TO SPECTRAL MODEL.

!     H. GUNTHER      ECMWF    JUNE 1990   MODIFIED FOR CYCLE_4.

!     J. BIDLOT       ECMWF    FEBRUARY 1996 MESSAGE PASSING

!     S. ABDALLA      ECMWF    OCTOBER 2001  AIR DENSITY AND Zi/L

!*    PURPOSE.
!     --------

!       EVALUATE WIND SPEED AND DIRECTION AT WAVE MODEL GRID POINTS.

!**   INTERFACE.
!     ----------

!     *CALL* *PREWIND (U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
!    &                 ROAIRO, ZIDLOLD, CICOVER,
!    &                 LLINIT,
!    &                 IREAD,
!    &                 NFIELDS, NGPTOTG, NC, NR,
!    &                 FIELDS, LWCUR, MASK_IN)*

!      *U10OLD*    REAL      INTERMEDIATE STORAGE OF MODULUS OF WIND
!                            VELOCITY.
!                            CLOCKWISE FROM NORTH).
!      *THWOLD*    REAL      INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
!                            WIND VELOCITY.
!      *USOLD*     REAL      INTERMEDIATE STORAGE OF MODULUS OF FRICTION
!                            VELOCITY.
!      *Z0OLD*     REAL      INTERMEDIATE STORAGE OF ROUGHNESS LENGTH IN
!                            M.
!      *TAUW*      REAL      WAVE STRESS IN (M/S)**2
!      *ROAIRO*    REAL      AIR DENSITY IN KG/M3.
!      *ZIDLOLD*   REAL      Zi/L (Zi: INVERSION HEIGHT,
!                                   L: MONIN-OBUKHOV LENGTH).
!      *CICOVER*   REAL      SEA ICE COVER.
!      *CITHICK*   REAL      SEA ICE THICKNESS.
!      *CIWA*      REAL      SEA ICE WAVE ATTENUATION FACTOR.
!      *LLINIT*    LOGICAL   TRUE IF WAMADSZIDL NEEDS TO BE CALLED.   
!      *LLALLOC_FIELDG_ONLY  LOGICAL IF TRUE THEN FIELDG DATA STRUCTURE WILL
!                            ONLY BE ALLOCATED BUT NOT INITIALISED.
!      *IREAD*     INTEGER   PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!                            (IF NEEDED).
!      *NFIELDS*   INTEGER   NUMBER OF FIELDS HOLDING ATMOSPHERIC DATA
!      *NGPTOTG*   INTEGER   NUMBER OF ATMOSPHERIC GRID POINTS
!      *NC*        INTEGER   NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!      *NR*        INTEGER   NUMBER OF ATM. ROWS OF LATITUDES
!      *FIELDS*    REAL      ATM. FIELDS (U10, V10, AIR DENSITY, Zi/L, U and V CURRENTS)
!      *LWCUR*     LOGICAL   INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS
!      *MASK_IN*   INTEGER   MASK TO INDICATE WHICH PART OF FIELDS IS RELEVANT.

!       *UNIT* *DESCRIPTION*

!          IU01    INPUT WIND DATA (SUB READWIND).
!          IU06    PRINTER OUTPUT (SUB INITMDL).
!          IUVELO  OUTPUT OF BLOCKED WIND FIELDS.
!          IUSCR   SCRATCH UNITS FOR ALL BLOCKS (INTERMEDIATE STORAGE,
!                  INPUT/OUTPUT) (SUB INITMDL).

!     METHOD.
!     -------

!       INPUT WIND FIELDS WHICH CAN BE YOWPONENTS OF
!                USTAR, U10, USTRESS
!       ARE TRANSFORMED TO FRICTION VELOCITIES.
!       THE INPUT FIELDS HAVE TO BE ON A LAT /LONG GRID.
!       SEE SUB READWIND FOR FORMATS AND HEADER INFORMATION,
!       WHICH HAVE TO BE GIVEN TO THE PROGRAM.

!       A DOUBLE LINEAR INTERPOLATION IN SPACE IS PERFORMED
!       ONTO THE MODEL BLOCKS.
!       IF THE WIND OUTPUT TIMSTEP IS LESS THAN THE INPUT TIMESTEP
!       A LINEAR INTERPOLATION IN TIME IS PERFORMED.

!       THERE ARE TWO POSSIBILITIES WITH RESPECT TO THE WIND
!       OUTPUT FILES:

!           1. PROPAGATION TIMESTEP >= WIND INPUT STEP
!              ONE OUTPUT FILE CONTAINS IDELPRO/IDELWO WINDFIELDS
!              I.E. INFORMATION FOR ONE PROPAGATION TIMESTEP.
!              TIME FILE(I+1)= TIME FILE(I)+ IDELPRO

!           2. PROPAGATION TIMESTEP < INPUT WIND TIMESTEP
!              ONE OUTPUT FILE CONTAINS IDELWI/IDELWO WINDFIELDS
!              I.E. INFORMATION FOR ONE WIND INPUT TIMESTEP.
!              TIME FILE(I+1)= TIME FILE(I) + IDELWI

!     EXTERNALS.
!     ----------

!       *ABORT*     - TERMINATES PROCESSING.
!       *AIRSEA*    - SURFACE LAYER STRESS.
!       *INCDAT*    - INCREMENTS DATE TIME GROUP.
!       *LOCINT*    - INTERPOLATES IN SPACE.
!       *NOTIM*     - STEERING SUB FOR INTERPOLATION IN SPACE ONLY.
!       *READWIND*   - READS A WIND FIELD.
!       *TIMIN*     - STEERING SUB FOR INTERPOLATION IN SPACE AND TIME.
!       *WAMWND*    - BLOCKS A WIND FIELD AND CONVERTS TO USTAR.

!     REFERENCE.
!     -----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWNEMOCOU,LWNEMOCOURECV
      USE YOWGRID  , ONLY : IJS      ,IJL
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NBLO     ,NFRE
      USE YOWSTAT  , ONLY : CDATEA   ,CDATEE   ,IDELPRO  ,IDELWI   ,
     &            IDELWO   ,LANAONLY, IDELT
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWTEXT  , ONLY : LRESTARTED
      USE YOWUNPOOL ,ONLY : LLUNSTR
      USE YOWWIND  , ONLY : CDA      ,CDAWIFL  ,FIELDG

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
#ifdef MODEL_COUPLING_ATM_WAV
      USE pgmcl_lib_WAM, only : WAV_coupl_prewind
#endif

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: NFIELDS, NGPTOTG, NC, NR
      INTEGER :: IDELWH, IREAD
      INTEGER :: ISTORE, IJ, IG, IIG
      INTEGER :: MASK_IN(NGPTOTG)

      REAL :: ZHOOK_HANDLE
      REAL :: FIELDS(NGPTOTG,NFIELDS)
      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD,USOLD,Z0OLD,TAUW,
     &                                  ROAIRO,ZIDLOLD,
     &                                  CICOVER,CITHICK
      REAL,DIMENSION(NINF:NSUP, NFRE, NBLO) :: CIWA 

      CHARACTER(LEN=14) :: CDTWIE, CDTWIS, ZERO

      LOGICAL :: LLFRSTNEMO
      LOGICAL :: LWCUR
      LOGICAL :: LLFRST
      LOGICAL :: LLINIT
      LOGICAL :: LLALLOC_FIELDG_ONLY, LLINIALL, LLOCAL

      DATA LLFRST / .TRUE. /
      DATA LLFRSTNEMO / .TRUE. /

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('PREWIND',0,ZHOOK_HANDLE)
#endif

!*    1. BEGIN AND END DATES OF WIND FIELDS TO BE PROCESSED.
!        ---------------------------------------------------

      ZERO = ' '

      IF (CDA.EQ.ZERO.OR.LANAONLY) THEN

!        IF START FROM PRESET FIELDS DO FIRST FIELD IN ADDITION.

        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) "  PREWIND AT CDATEA = ", CDATEA 
          FLUSH (IU06)
        ENDIF
        CDTWIS = CDATEA
      ELSE
        CDTWIS = CDAWIFL
        IDELWH = -MAX(IDELPRO,IDELWI)+IDELWO
        CALL INCDATE (CDTWIS,IDELWH)
      ENDIF

      IF(CDAWIFL.LT.CDATEE) THEN
        CDTWIE = CDAWIFL
        LLFRST = .FALSE.
      ELSE
        CDTWIE = CDATEE
      ENDIF
! ----------------------------------------------------------------------

!*    2. PROCESS FORCING FIELDS.
!        -----------------------

!*    2.0 GLOBAL FIELD FOR THE INPUTS
!         ---------------------------
      LLINIALL=.TRUE.
      LLOCAL=.TRUE.
      CALL INIT_FIELDG(LLALLOC_FIELDG_ONLY,LLINIALL,LLOCAL)

!     2.1 IN COUPLED RUNS, TRANSFORM INPUT FORCING FIELDS TO WAM GRID.
!         -----------------------------------------------------------

      CALL IFSTOWAM (NFIELDS, NGPTOTG, NC, NR,
     &               FIELDS, LWCUR, MASK_IN)

!     2.1.1 GET DATA FROM NEMO.
!           ------------------
!
!     IF RESTART DO NOT CALL RECVNEMOFIELDS AT THE FIRST TIME STEP
!     BUT USE THE VALUES FROM THE RESTART.
!      
      IF(LWNEMOCOU.AND.LWNEMOCOURECV) THEN
        CALL RECVNEMOFIELDS(LLFRSTNEMO.AND.LRESTARTED,
     &                      LLFRSTNEMO.AND..NOT.LRESTARTED)
        LLFRSTNEMO=.FALSE.
      ENDIF

!!!   WHEN COUPLED, IT IS ALSO NEEDED TO INITIALISE THE AIR DENSITY AND
!!!   THE CONVECTIVE VELOCITY SCALE ARRAYS SINCE THESE ARE NOT (YET) PROVIDED
!!!   AS PART OF THE GRIB RESTART FILES.
      IF (LLINIT) CALL WAMADSZIDL(ROAIRO,ZIDLOLD)

#if !defined MODEL_COUPLING_ATM_WAV && !defined MODEL_COUPLING_OCN_WAV

!     2.2 GET SURFACE CURRENTS TO WAM BLOCK STRUCTURE (if needed) 
!         -------------------------------------------

      CALL GETCURR(LWCUR, IREAD)


!*    PROCESS THE OTHER FORCING FIELDS.
!     ---------------------------------
      IF (IDELWO.GE.IDELWI) THEN

!*    2.2 NO TIME INTERPOLATION.
!         ----------------------

        IF (ITEST.GE.2) THEN
          WRITE (IU06,*) '   SUB. PREWIND: WIND REQUEST'
          WRITE (IU06,*) '     NO TIME INTERPOLATION'
          WRITE (IU06,*) '     START OF PERIOD IS    CDTWIS = ',CDTWIS
          WRITE (IU06,*) '     END   OF PERIOD IS    CDTWIE = ',CDTWIE
          WRITE (IU06,*) '     WIND INPUT TIME STEP  IDELWI = ',IDELWI
          WRITE (IU06,*) '     WIND OUTPUT TIME STEP IDELWO = ',IDELWO
          FLUSH (IU06)
        ENDIF
        CALL NOTIM (CDTWIS, CDTWIE,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     &              ROAIRO, ZIDLOLD, CICOVER, CITHICK, IREAD,
     &              LWCUR)

      ELSE

!*    2.3 TIME INTERPOLATION.
!         -------------------

        IF (ITEST.GE.2) THEN
          WRITE (IU06,*) '   SUB. PREWIND: WIND REQUEST'
          WRITE (IU06,*) '     TIME INTERPOLATION'
          WRITE (IU06,*) '     START OF PERIOD IS    CDTWIS = ',CDTWIS
          WRITE (IU06,*) '     END   OF PERIOD IS    CDTWIE = ',CDTWIE
          WRITE (IU06,*) '     WIND INPUT TIME STEP  IDELWI = ',IDELWI
          WRITE (IU06,*) '     WIND OUTPUT TIME STEP IDELWO = ',IDELWO
          FLUSH (IU06)
        ENDIF
        CALL TIMIN (CDTWIS, CDTWIE,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     &              ROAIRO, ZIDLOLD, CICOVER, CITHICK, IREAD,
     &              LWCUR)

      ENDIF

#else
      CALL WAV_coupl_prewind(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,          
     &                       ROAIRO, ZIDLOLD,                         
     &                       CICOVER, CITHICK, CIWA,                  
     &                       LLINIT,                                  
     &                       IREAD,                                   
     &                       NFIELDS, NGPTOTG, NC, NR,                
     &                       FIELDS, LWCUR, MASK_IN)
#endif


!*    2.5 DEALLOCATE GRID ARRAYS FOR INPUT FORCING FIELDS
!         -----------------------------------------------

      DEALLOCATE(FIELDG)

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('PREWIND',1,ZHOOK_HANDLE)
#endif

      RETURN

      END SUBROUTINE PREWIND
