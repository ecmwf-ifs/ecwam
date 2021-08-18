      SUBROUTINE PREWIND (U10OLD, THWOLD, USOLD, Z0OLD,                &
     &                    ROAIRO, WSTAROLD,                            &
     &                    CICOVER, CITHICK,                            &
     &                    LLINIT, LLALLOC_FIELDG_ONLY,                 &
     &                    IREAD,                                       &
     &                    NFIELDS, NGPTOTG, NC, NR,                    &
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

!     *CALL* *PREWIND (U10OLD,THWOLD,USOLD,Z0OLD,
!    &                 ROAIRO, WSTAROLD, CICOVER,
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
!      *ROAIRO*    REAL      AIR DENSITY IN KG/M3.
!      *WSTAROLD*  REAL      CONVECTIVE VELOCITY M/S.
!      *CICOVER*   REAL      SEA ICE COVER.
!      *CITHICK*   REAL      SEA ICE THICKNESS.
!      *LLINIT*    LOGICAL   TRUE IF WAMADSWSTAR NEEDS TO BE CALLED.   
!      *LLALLOC_FIELDG_ONLY  LOGICAL IF TRUE THEN FIELDG DATA STRUCTURE WILL
!                            ONLY BE ALLOCATED BUT NOT INITIALISED.
!      *IREAD*     INTEGER   PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!                            (IF NEEDED).
!      *NFIELDS*   INTEGER   NUMBER OF FIELDS HOLDING ATMOSPHERIC DATA
!      *NGPTOTG*   INTEGER   NUMBER OF ATMOSPHERIC GRID POINTS
!      *NC*        INTEGER   NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!      *NR*        INTEGER   NUMBER OF ATM. ROWS OF LATITUDES
!      *FIELDS*    REAL      ATM. FIELDS (U10, V10, AIR DENSITY, w*, U and V CURRENTS)
!      *LWCUR*     LOGICAL   INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS
!      *MASK_IN*   INTEGER   MASK TO INDICATE WHICH PART OF FIELDS IS RELEVANT.

!       *UNIT* *DESCRIPTION*

!          IU01    INPUT WIND DATA (SUB READWIND).
!          IU06    PRINTER OUTPUT (SUB INITMDL).
!          IUVELO  OUTPUT OF BLOCKED WIND FIELDS.

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
!       *WAMWND*    - BLOCKS A WIND FIELD AND CONVERTS TO USTAR.

!     REFERENCE.
!     -----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWNEMOCOU,LWNEMOCOURECV
      USE YOWGRID  , ONLY : IJS      ,IJL
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWSTAT  , ONLY : CDATEA   ,CDATEE   ,IDELPRO  ,IDELWI   ,    &
     &            IDELWO   ,LANAONLY, IDELT
      USE YOWTEST  , ONLY : IU06
      USE YOWTEXT  , ONLY : LRESTARTED
!      USE YOWWIND  , ONLY : CDATEWL   ,CDAWIFL  ,FIELDG  ,FF_NEXT
      USE YOWWIND  , ONLY : FORCING_FIELDS, CDATEWL   ,CDAWIFL  ,FIELDG

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "getcurr.intfb.h"
#include "getfrstwnd.intfb.h"
#include "ifstowam.intfb.h"
#include "incdate.intfb.h"
#include "init_fieldg.intfb.h"
#include "notim.intfb.h"
#include "recvnemofields.intfb.h"
#include "wamadswstar.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: NFIELDS
      INTEGER(KIND=JWIM), INTENT(IN) :: NGPTOTG
      INTEGER(KIND=JWIM), INTENT(IN) :: NC
      INTEGER(KIND=JWIM), INTENT(IN) :: NR
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      LOGICAL, INTENT(IN) :: LWCUR
      LOGICAL, INTENT(IN) :: LLINIT
      LOGICAL, INTENT(IN) :: LLALLOC_FIELDG_ONLY
      INTEGER(KIND=JWIM),DIMENSION(NGPTOTG), INTENT(INOUT)  :: MASK_IN
      REAL(KIND=JWRB),DIMENSION(NGPTOTG,NFIELDS), INTENT(IN) :: FIELDS
      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(INOUT) :: U10OLD
      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(INOUT) :: THWOLD
      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(INOUT) :: USOLD
      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(INOUT) :: Z0OLD
      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(INOUT) :: ROAIRO
      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(INOUT) :: WSTAROLD
      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(INOUT) :: CICOVER
      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(INOUT) :: CITHICK


      INTEGER(KIND=JWIM) :: IDELWH
      INTEGER(KIND=JWIM) :: ISTORE, IJ

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      TYPE(FORCING_FIELDS), DIMENSION(IJS:IJL) :: FF_NEXT

      CHARACTER(LEN=14) :: CDTWIE, CDTWIS, ZERO

      LOGICAL, SAVE :: LLFRSTNEMO
      LOGICAL :: LLINIALL, LLOCAL, LLMORE

      DATA LLFRSTNEMO / .TRUE. /

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('PREWIND',0,ZHOOK_HANDLE)

!*    1. BEGIN AND END DATES OF WIND FIELDS TO BE PROCESSED.
!        ---------------------------------------------------

      ZERO = ' '

      IF (CDATEWL == ZERO .OR. LANAONLY) THEN
!       IF START FROM PRESET FIELDS DO FIRST FIELD IN ADDITION.
        CDTWIS = CDATEA
      ELSE
        CDTWIS = CDAWIFL
        IDELWH = -MAX(IDELPRO,IDELWI)+IDELWO
        CALL INCDATE (CDTWIS,IDELWH)
      ENDIF

      IF (CDAWIFL < CDATEE) THEN
        CDTWIE = CDAWIFL
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

      CALL IFSTOWAM (NFIELDS, NGPTOTG, NC, NR, FIELDS, LWCUR, MASK_IN)

!     2.1.1 GET DATA FROM NEMO (OR BINARY RESTART).
!           -------------------------------------
 
      IF (LWNEMOCOU .AND. LWNEMOCOURECV) THEN
        CALL RECVNEMOFIELDS(LLFRSTNEMO.AND.LRESTARTED,                  &
     &                      LLFRSTNEMO.AND..NOT.LRESTARTED)
        LLFRSTNEMO=.FALSE.
      ENDIF

!!!   WHEN COUPLED, IT IS ALSO NEEDED TO INITIALISE THE AIR DENSITY AND
!!!   THE CONVECTIVE VELOCITY SCALE ARRAYS SINCE THESE ARE NOT (YET) PROVIDED
!!!   AS PART OF THE GRIB RESTART FILES.
      IF (LLINIT) CALL WAMADSWSTAR(IJS, IJL, ROAIRO, WSTAROLD)

!     2.2 GET SURFACE CURRENTS TO WAM BLOCK STRUCTURE (if needed) 
!         -------------------------------------------

      CALL GETCURR(LWCUR, IREAD)

!*    PROCESS THE OTHER FORCING FIELDS.
!     ---------------------------------
      IF (IDELWO >= IDELWI) THEN

!*      2.2 NO TIME INTERPOLATION.
!       ----------------------

        IF (CDATEWL == ZERO) THEN
!         Initialisation (either first time or following a restart)
          CALL GETFRSTWND (CDTWIS, CDTWIE,                       &
     &                     IJS, IJL,                             &
     &                     U10OLD, THWOLD, USOLD, Z0OLD,         &
     &                     ROAIRO, WSTAROLD, CICOVER, CITHICK,   &
     &                     IREAD, LWCUR, LLMORE)
        ELSE
          LLMORE = .TRUE.
        ENDIF


        IF (LLMORE) THEN
!         Update forcing
          CALL NOTIM (CDTWIS, CDTWIE,                       &
     &                IJS, IJL, FF_NEXT,                    &
     &                IREAD, LWCUR)
        ENDIF
      ELSE

!*      2.3 TIME INTERPOLATION.
!       -------------------

              WRITE (IU06,*) ' *********************************'
              WRITE (IU06,*) ' *                               *'
              WRITE (IU06,*) ' * PROBLEM IN PREWIND .........  *'
              WRITE (IU06,*) ' * NOT AN OPTION ANY MORE        *'
              WRITE (IU06,*) ' *                               *'
              WRITE (IU06,*) ' *********************************'
              CALL FLUSH(IU06)
              CALL ABORT1

      ENDIF


!*    2.5 DEALLOCATE GRID ARRAYS FOR INPUT FORCING FIELDS
!         -----------------------------------------------

      DEALLOCATE(FIELDG)

      IF (LHOOK) CALL DR_HOOK('PREWIND',1,ZHOOK_HANDLE)

      END SUBROUTINE PREWIND
