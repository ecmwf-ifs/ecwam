MODULE PARKIND_WAVE
!
!     *** Define usual kinds ***
!
USE PARKIND1, ONLY : JPIM, JPRB

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
#else
INTEGER, PARAMETER :: JWRB = 4 
#endif
INTEGER, PARAMETER :: JWRU = 8 
!

END MODULE PARKIND_WAVE
