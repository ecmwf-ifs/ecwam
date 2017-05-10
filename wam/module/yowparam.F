      MODULE YOWPARAM

      IMPLICIT NONE

      INTEGER, PARAMETER :: IMDLGRDID=218
      INTEGER, PARAMETER :: KWAMVER=6

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

      REAL    SWAMPWIND 
      REAL    SWAMPWIND2 
      REAL    DTNEWWIND
      REAL    SWAMPCIFR
      REAL    SWAMPCITH

      CHARACTER*1 CLDOMAIN

      LOGICAL LTURN90
      LOGICAL :: LL1D
      LOGICAL :: LWDINTS

!*     VARIABLE.   TYPE.     PURPOSE.
!       ---------   -------   --------
!      *IMDLGRDID* INTEGER   MODEL GRID IDENTIFICATION NUMBER.
!                            IT SHOULD BE UPDATED ONLY IF NEW gridglou
!                            AND ubufglou FILES ARE PRODUCED!!!!
!      *KWAMVER*   INTEGER   WAM SOFTWARE VERSION NUMBER.
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
!      *SWAMPWIND* REAL      CONSTANT WIND SPEED USED TO RUN SWAMP CASE.
!                            FIRST VALUE
!      *SWAMPWIND2*REAL      CONSTANT WIND SPEED USED TO RUN SWAMP CASE.
!                            SECOND VALUE PROVIDED IT'S NOT ZERO
!      *DTNEWWIND* REAL      TIME AFTER WHICH SWAMPWIND2 IS APPLIED
!                            IN HOURS
!      *SWAMPCIFR* REAL      SEA ICE COVER FOR THE NORTHERN HALF OF THE SWAMP DOMAIN 
!      *SWAMPCITH* REAL      SEA ICE THICKNESS FOR THE NORTHERN HALF OF THE SWAMP DOMAIN 
!      *LTURN90*   LOGICAL   IF TRUE THE NEW SWAMP WIND WILL BE TURNED
!                            BY 90 DEGREES.
!      *LWDINTS*   LOGICAL   IF TRUE A FILE CONTAINING WIND SPEED (m/s) AND DIRECTION
!                            (in degrees following the meteorological convention)
!                            TIME SERIES WILL BY USED TO SET THE WIND FORCING
!                            OVER THE ALL DOMAIN.
!                            THE PRESENCE OF THE FILE, windforcing_time_series
!                            WILL DETERMINE WHETHER OR NOT LWDINTS IS TRUE OR NOT
!      *CLDOMAIN*  CHARACTER DEFINES THE DOMAIN OF THE MODEL (for the
!                            FDB and for selection of some variables)
!      *LL1D*      LOGICAL   IF TRUE THEN THE DOMAIN DECOMPOSITION IS ONLY
!                            DONE IN LATITUNAL BANDS
!                            (like it used to be done).
! ----------------------------------------------------------------------

      END MODULE YOWPARAM
