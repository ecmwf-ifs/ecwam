      MODULE YOWINDN

      IMPLICIT NONE

!*    ** *INDNL* - INDICES AND WEIGHTS USED IN THE COMPUTATION
!                  OF THE NONLINEAR TRANSFER RATE.

      INTEGER, PARAMETER :: KFRH=30 

      INTEGER, ALLOCATABLE :: IKP(:) 
      INTEGER, ALLOCATABLE :: IKP1(:) 
      INTEGER, ALLOCATABLE :: IKM(:) 
      INTEGER, ALLOCATABLE :: IKM1(:) 
      INTEGER, ALLOCATABLE :: K1W(:,:) 
      INTEGER, ALLOCATABLE :: K2W(:,:) 
      INTEGER, ALLOCATABLE :: K11W(:,:) 
      INTEGER, ALLOCATABLE :: K21W(:,:) 
      REAL, ALLOCATABLE :: AF11(:) 
      REAL, ALLOCATABLE :: FKLAP(:) 
      REAL, ALLOCATABLE :: FKLAP1(:) 
      REAL, ALLOCATABLE :: FKLAM(:) 
      REAL, ALLOCATABLE :: FKLAM1(:) 
      REAL ACL1 
      REAL ACL2 
      REAL CL11 
      REAL CL21 
      REAL DAL1 
      REAL DAL2 
      REAL, ALLOCATABLE :: FRH(:)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   -------
!      *KFRH*      INTEGER   SIZE OF FRH
!      *IKP*       INTEGER   FREQUENCY INDEX ARRAY FOR STORING ENERGY
!                            TRANSFER INCREMENTS INTO BINS, WAVE NO. 3.
!      *IKP1*      INTEGER   IKP+1.
!      *IKM*       INTEGER   FREQUENCY INDEX ARRAY FOR STORING ENERGY
!                            TRANSFER INCREMENTS INTO BINS, WAVE NO. 4.
!      *IKM1*      INTEGER   IKM+1
!      *K1W*       INTEGER   ANGULAR INDEX ARRAY FOR STORING ENERGY
!                            TRANSFER INCREMENTS INTO BINS, WAVE NO. 3.
!      *K11W*      INTEGER   K1W(.,1)-1, K1W(.,2)+1.
!      *K2W*       INTEGER   ANGULAR INDEX ARRAY FOR STORING ENERGY
!                            TRANSFER INCREMENTS INTO BINS, WAVE NO. 4.
!      *K21W*      INTEGER   K2W(.,1)+1, K2W(.,2)-1.
!      *AF11*      REAL      WEIGHTS FOR DISCRETE APPROXIMATION OF NONL
!                            TRANSFER (AT PRESENT ONE TERM ONLY SET TO
!                            3000). MULTIPLIED BY FREQUENCIES **11.
!      *FKLAP*     REAL      WEIGHT IN FREQUENCY GRID FOR INTERPOLATION,
!                            WAVE NO. 3 ("1+LAMBDA" TERM).
!      *FKLAP1*    REAL      1-FKLAP.
!      *FKLAM*     REAL      WEIGHT IN FREQUENCY GRID FOR INTERPOLATION,
!                            WAVE NO. 4 ("1-LAMBDA" TERM).
!      *ACL1*      REAL      WEIGHT IN ANGULAR GRID FOR INTERPOLATION,
!                            WAVE NO. 3 ("1+LAMBDA" TERM).
!      *ACL2*      REAL      WEIGHT IN ANGULAR GRID FOR INTERPOLATION,
!                            WAVE NO. 4 ("1-LAMBDA" TERM).
!      *CL11*      REAL      1.-ACL1.
!      *CL21*      REAL      1.-ACL2.
!      *DAL1*      REAL      1./ACL1.
!      *DAL2*      REAL      1./ACL2.
!      *FRH*       REAL      TAIL FREQUENCY RATION **5

! ----------------------------------------------------------------------
      END MODULE YOWINDN
