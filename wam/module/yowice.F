      MODULE YOWICE

      IMPLICIT NONE

!*     ** *ICE* ICE POINTS

      LOGICAL LICERUN 
      INTEGER IPARAMCI
      INTEGER, ALLOCATABLE :: ICEMASK(:,:) 
      REAL :: FLMIN
      REAL, ALLOCATABLE :: FLMINFR(:,:)
      REAL, ALLOCATABLE :: CI(:,:)

!--------------------------------------------------------------------

!*    VARIABLE     TYPE      PURPOSE
!     --------     ----      -------
!     LICERUN      LOGICAL   SET TO TRUE IF ICE IS TAKEN INTO ACCOUNT
!     IPARAMCI     INTEGER   GRIB PARAMETER VALUE OF SEA ICE OR SST.
!     ICEMASK      INTEGER   SET TO 0 FOR ALL LOCAL SEA POINTS WHICH
!                            ARE ICE, 1 FOR ALL OTHER SEA POINTS.
!     FLMIN        REAL      SETS THE MAXIMUM OF THE MINIMUM VALUE
!                            OF EACH SPECTRAL COMPONENTS.
!     FLMINFR      REAL      ARRAY CONTAINING THE MINIMUM VALUE
!                            OF EACH SPECTRAL COMPONENTS FOR A GIVEN
!                            FREQUENCY.
!     CI           REAL      INPUT SEA ICE (in COUPLED MODEL)
!--------------------------------------------------------------------

      END MODULE YOWICE
