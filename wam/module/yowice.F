      MODULE YOWICE

      IMPLICIT NONE

!*     ** *ICE* ICE POINTS

      LOGICAL LICERUN 
      INTEGER, ALLOCATABLE :: ICEMASK(:,:) 
      REAL, PARAMETER :: FLMIN = 0.0001
      REAL, ALLOCATABLE :: FLMINFR(:,:)

!--------------------------------------------------------------------

!*    VARIABLE     TYPE      PURPOSE
!     --------     ----      -------
!     ICEMASK      INTEGER   SET TO 0 FOR ALL LOCAL SEA POINTS WHICH
!                            ARE ICE, 1 FOR ALL OTHER SEA POINTS.
!     LICERUN      LOGICAL   SET TO TRUE IF ICE IS TAKEN INTO ACCOUNT
!     FLMIN        REAL      SETS THE MAXIMUM OF THE MINIMUM VALUE
!                            OF EACH SPECTRAL COMPONENTS.
!     FLMINFR      REAL      ARRAY CONTAINING THE MINIMUM VALUE
!                            OF EACH SPECTRAL COMPONENTS FOR A GIVEN
!                            FREQUENCY.
!--------------------------------------------------------------------

      END MODULE YOWICE
