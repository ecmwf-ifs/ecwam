      MODULE YOWSHAL

      IMPLICIT NONE

!*    ** *SHALLOW*   SHALLOW WATER TABLES.

      INTEGER,PARAMETER :: NDEPTH=52 

      REAL, ALLOCATABLE :: DEPTH(:,:) 
      REAL DEPTHA 
      REAL DEPTHD 
      REAL, ALLOCATABLE :: TCGOND(:,:) 
      REAL, ALLOCATABLE :: TFAK(:,:) 
      REAL, ALLOCATABLE :: TSIHKD(:,:) 
      INTEGER, ALLOCATABLE :: INDEP(:) 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *NDEPTH*    INTEGER   LENGTH OF SHALLOW WATER TABLES.
!      *DEPTH*     REAL      WATER DEPTH IN METRES.
!      *DEPTHA*    REAL      MINIMUM DEPTH FOR TABLES (METRES).
!      *DEPTHD*    REAL      DEPTH INCREMENT (METRES).
!      *TCGOND*    REAL      SHALLOW WATER GROUP VELOCITY TABLE.
!      *TFAK*      REAL      WAVE NUMBER TABLE.
!      *TSIHKD*    REAL      TABLE FOR OMEGA/SINH(2KD).
!      *INDEP*     INTEGER   DEPTH INDEX FOR A BLOCK.

! ----------------------------------------------------------------------
      END MODULE YOWSHAL
