      MODULE YOWPARAM

      IMPLICIT NONE

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

      CHARACTER*1 CLDOMAIN

!*     VARIABLE.   TYPE.     PURPOSE.
!       ---------   -------   --------
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
!      *CLDOMAIN*  CHARACTER DEFINES THE DOMAIN OF THE MODEL (for the
!                            FDB and for selection of some variables)

C ----------------------------------------------------------------------

      END MODULE YOWPARAM

