      MODULE YOWGRID

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *GRIDPAR*  GENERAL GRID INFORMATION.

      REAL(KIND=JWRB)              :: DELPHI 
      REAL(KIND=JWRB), ALLOCATABLE :: DELLAM(:) 
      REAL(KIND=JWRB), ALLOCATABLE :: SINPH(:) 
      REAL(KIND=JWRB), ALLOCATABLE :: COSPH(:) 
      REAL(KIND=JWRB), ALLOCATABLE :: CDR(:,:) 
      REAL(KIND=JWRB), ALLOCATABLE :: SDR(:,:) 
      REAL(KIND=JWRB), ALLOCATABLE :: PRQRT(:)

      INTEGER(KIND=JWIM) :: IJS
      INTEGER(KIND=JWIM) :: IJL
      INTEGER(KIND=JWIM) :: IJSLOC, IJLLOC, IJGLOBAL_OFFSET
      INTEGER(KIND=JWIM) :: NPROMA_WAM
      INTEGER(KIND=JWIM) :: NBLOC
      INTEGER(KIND=JWIM), ALLOCATABLE :: IJFROMBLOC

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *DELPHI*    REAL      GRID INCREMENT FOR LATITUDE (METRES).
!      *DELLAM*    REAL      GRID INCREMENT FOR LONGITUDE AT EQUATOR
!                            IN METRES.
!      *SINPH*     REAL      SIN OF LATITUDE.
!      *COSPH*     REAL      COS OF LATITUDE.
!      *CDR*       REAL      ADVECTION VELOCITY IN SE-NW DIRECTION WHEN
!                            GRID ROTATED BY 45 DEGREES.
!      *SDR*       REAL      ADVECTION VELOCITY IN SW-NE DIRECTION WHEN
!                            GRID ROTATED BY 45 DEGREES.
!      *PRQRT*     REAL      PERCENTAGE OF THE ROTATED (45 degree) GRID
!                            USED FOR THE PROPAGATION (see PROPAGS1).
!                            IT IS 50% EXCEPT IF FOR STABILITY REASON
!                            IT HAS TO BE DECREASED ABOVE ~75 DEGREE

!      *IJS*       INTEGER   INDEX OF FIRST POINT ON A GIVEN PROCESSOR
!      *IJL*       INTEGER   INDEX OF LAST POINT ON A GIVEN PROCESSOR

!!!!!! for the unstructured part of the code, it was coded so that halo points
!!!!!! are included as part of the points on a given processor (at the end of the vector of points).
!!!!!! For output, it was necessary to introduce IJSLOC and IJLLOC to point to the points that are purely local
!      *IJSLOC*    INTEGER   INDEX OF FIRST LOCAL POINT
!      *IJLLOC*    INTEGER   INDEX OF LAST LOCAL POINT
!      *IJGLOBAL_OFFSET INTEGER OFFSET TO PLACE FIRST LOCAL POINT IN GLOBAL ARRAY

!!!!! if not unstructured then IJSLOC = IJS and IJLLOC = IJL

!      FOR OPENMP:
!      *NPROMA_WAM* INTEGER  MAX NUMBERS OF GRID POINTS PER THREAD
!      *NBLOC*      INTEGER  NUMBER OF BLOCKS OF MAXIMUM NPROMA_WAM POINTS

! ----------------------------------------------------------------------
      END MODULE YOWGRID
