      MODULE YOWPARAM

      IMPLICIT NONE

      INTEGER, PARAMETER :: IMDLGRDID=215

!*    ** *PARAM*  GENERAL GRID SIZE INFORMATION.

      INTEGER NANG 
      INTEGER NFRE
      INTEGER NGX
      INTEGER NGY
      INTEGER NBLO
      INTEGER NIBLO
      INTEGER NOVER 
      INTEGER NIBL1 
      INTEGER NIBLD 
      INTEGER NBLD 
      INTEGER NIBLC 
      INTEGER NBLC 

      REAL    XKAPPA2
      REAL    SWAMPWIND 
      REAL    SWAMPWIND2 
      REAL    DTNEWWIND
      REAL    HSCOEFCOR
      REAL    HSCONSCOR

      CHARACTER*1 CLDOMAIN

      LOGICAL LTURN90
      LOGICAL :: LL1D

!*     VARIABLE.   TYPE.     PURPOSE.
!       ---------   -------   --------
!      *IMDLGRDID* INTEGER   MODEL IDENTIFICATION NUMBER.
!                            IT SHOULD BE UPDATED FOR EACH NEW CYCLE!!!!
!      *NANG*      INTEGER   NUMBER OF ANGLES.
!      *NFRE*      INTEGER   NUMBER OF FREQUENCIES.
!      *NGX*       INTEGER   NUMBER OF LONGITUDES IN GRID.
!      *NGY*       INTEGER   NUMBER OF LATITUDES IN GRID.
!      *NBLO*      INTEGER   NUMBER OF BLOCKS.
!      *NIBLO*     INTEGER   NUMBER OF SEA POINTS IN BLOCK.
!      *NOVER*     INTEGER   MAXIMUM NUMBER POINTS IN FIRST LATITUDE
!                            OF BLOCKS.
!      *NIBL1*     INTEGER   = NIBLO IF MULTI BLOCK VERSION.
!                            =     1 IF ONE BLOCK VERSION.
!      *NIBLD*     INTEGER   = NIBLO IF DEPTH OR CURRENT REFRACTION.
!                            = 1     ELSE.
!      *NBLD*      INTEGER   = NBLO  IF DEPTH OR CURRENT REFRACTION.
!                            = 1     ELSE.
!      *NIBLC*     INTEGER   = NIBLO IF CURRENT REFRACTION.
!                            = 1     ELSE.
!      *NBLC*      INTEGER   = NBLO  IF CURRENT REFRACTION.
!                            = 1     ELSE.
!      *XKAPPA2*   REAL      KAPPA2 PARAMETER USED IN THE DETERMINATION
!                            OF ALTIMETER WAVE HEIGHTS (SEE GRFIELD)
!      *SWAMPWIND* REAL      CONSTANT WIND SPEED USED TO RUN SWAMP CASE.
!                            FIRST VALUE
!      *SWAMPWIND2*REAL      CONSTANT WIND SPEED USED TO RUN SWAMP CASE.
!                            SECOND VALUE PROVIDED IT'S NOT ZERO
!      *DTNEWWIND* REAL      TIME AFTER WHICH SWAMPWIND2 IS APPLIED
!                            IN HOURS
!      *LTURN90*   LOGICAL   IF TRUE THE NEW SWAMP WIND WILL BE TURNED
!                            BY 90 DEGREES.
!      *HSCOEFCOR* REAL      COEFFICIENT OF THE CORRECTIVE LINEAR
!                            REGRESSION FOR ALTIMETER WAVE HEIGHTS.
!      *HSCONSCOR* REAL      CONSTANT COEFFICIE$NT  OF THE CORRECTIVE
!                            LINEAR REGRESSION FOR ALTIMETER WAVE
!                            HEIGHTS.
!      *CLDOMAIN*  CHARACTER DEFINES THE DOMAIN OF THE MODEL (for the
!                            FDB and for selection of some variables)
!      *LL1D*      LOGICAL   IF TRUE THEN THE DOMAIN DECOMPOSITION IS ONLY
!                            DONE IN LATITUNAL BANDS
!                            (like it used to be done).
! ----------------------------------------------------------------------

      END MODULE YOWPARAM
