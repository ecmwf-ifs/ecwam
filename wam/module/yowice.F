      MODULE YOWICE

      IMPLICIT NONE

!*     ** *ICE* ICE POINTS

      LOGICAL LICERUN 
      INTEGER, ALLOCATABLE :: ICEMASK(:,:) 

!--------------------------------------------------------------------

!*    VARIABLE     TYPE      PURPOSE
!     --------     ----      -------
!     ICEMASK      INTEGER   SET TO 0 FOR ALL LOCAL SEA POINTS WHICH
!                            ARE ICE, 1 FOR ALL OTHER SEA POINTS.
!     LICERUN      LOGICAL   SET TO TRUE IF ICE IS TAKEN INTO ACCOUNT

!--------------------------------------------------------------------

      END MODULE YOWICE
