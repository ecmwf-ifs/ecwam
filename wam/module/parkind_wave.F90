MODULE PARKIND_WAVE
!
!     *** Define usual kinds ***
!
USE PARKIND1, ONLY : JPIM, JPRB, JPRD

IMPLICIT NONE
SAVE
!
!     Integer Kinds
!     -------------
!
INTEGER, PARAMETER :: JWIM = JPIM

!
!     Real Kinds
!     ----------
!
#ifdef ECMWF
INTEGER, PARAMETER :: JWRB = JPRB
INTEGER, PARAMETER :: JWRU = JPRD 
#else
INTEGER, PARAMETER :: JWRB = 4 
INTEGER, PARAMETER :: JWRU = 8 
#endif
!

END MODULE PARKIND_WAVE
