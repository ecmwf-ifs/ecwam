      MODULE YOWCOER

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *COERS* - CONTROLS INFORMATION OF ERS OUTPUT.          

      INTEGER(KIND=JWIM)              :: NERS
      INTEGER(KIND=JWIM)              :: IDELERS 
      INTEGER(KIND=JWIM)              :: IERS 
      INTEGER(KIND=JWIM), ALLOCATABLE :: IJERS(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: IGERS(:) 

      CHARACTER(LEN=14)               :: CDTERS 

      END MODULE YOWCOER