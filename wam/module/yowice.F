      MODULE YOWICE

      IMPLICIT NONE

!*     ** *ICE* ICE POINTS

      LOGICAL LICERUN 
      INTEGER, ALLOCATABLE :: ICEMASK(:,:) 
      REAL, PARAMETER :: FLMIN = 0.0001

!--------------------------------------------------------------------

!*    VARIABLE     TYPE      PURPOSE
!     --------     ----      -------
!     ICEMASK      INTEGER   SET TO 0 FOR ALL LOCAL SEA POINTS WHICH
!                            ARE ICE, 1 FOR ALL OTHER SEA POINTS.
!     LICERUN      LOGICAL   SET TO TRUE IF ICE IS TAKEN INTO ACCOUNT
!     FLMIN        REAL      SETS THE MINIMUM VALUE OF EACH SPECTRAL
!                            COMPONENTS.
!                            THIS MAKES HS=3cm WHEN U10 = 0.0
!--------------------------------------------------------------------

      END MODULE YOWICE
