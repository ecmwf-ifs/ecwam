      MODULE YOWTEST

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *TESTO* - PRINTER OUTPUT UNIT AND TEST FLAGS.

      INTEGER(KIND=JWIM) :: IU06 
      INTEGER(KIND=JWIM) :: ITEST 
      INTEGER(KIND=JWIM) :: ITESTB 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IU06*      INTEGER   UNIT FOR PRINTER OUTPUT.
!      *ITEST*     INTEGER   TEST OUTPUT LEVEL:
!                             .LE. 0  NO OUTPUT
!                             .GE. I  OUTPUT TO SUB. LEVEL I
!      *ITESTB*    INTEGER   MAX BLOCK NUMBER FOR OUTPUT IN BLOCK LOOPS

! ----------------------------------------------------------------------
      END MODULE YOWTEST
