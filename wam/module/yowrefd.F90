      MODULE YOWREFD

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *REFDOT* - DEPTH AND CURRENT PART OF THETA DOT.

      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:)   :: THDD
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:)   :: THDC
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:,:) :: SDOT

      LOGICAL :: LLUPDTTD

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *THDD*      REAL      DEPTH GRADIENT PART OF THETA DOT.
!      *THDC*      REAL      CURRENT GRADIENT PART OF THETA DOT.
!      *SDOT*      SIGMA DOT ARRAY

!      *LLUPDTTD*  LOGICAL  IF TRUE THETA DOT ARRAYS NEED UPDATING.

! ----------------------------------------------------------------------
      END MODULE YOWREFD
