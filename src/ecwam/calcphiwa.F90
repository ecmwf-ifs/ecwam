FUNCTION CALCPHIWA(SPOS,SNEG,DSII,SIG) RESULT(PHIWA)

! ----------------------------------------------------------------------------
!
!  1. Purpose :
!
!      Calculate energy flux from wind into waves, obtained from wind-energy-input (Sin).
!
!                            / FRMAX
!      tau = g * rho_water * | Sin(f) df
!                            /

!----------------------------------------------------------------------
!
!     INTERFACE VARIABLES.
!     --------------------

!     ORIGIN.
!     ----------
!     Adapted from TAUWINDS
!     Implementation into ECWAM DECEMBER 2021 by J. Kousal 

! ----------------------------------------------------------------------------
!

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWPCONS , ONLY : G        ,ROWATER, ZPI
      USE YOWFRED  , ONLY : DELTH, FRATIO
      USE YOWPARAM , ONLY : NANG     ,NFRE    
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

!----------------------------------------------------------------------

      IMPLICIT NONE
#include "irange.intfb.h"

      REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN)  :: SPOS ! POS Sin(sigma) in [m2/rad-Hz]
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN)  :: SNEG ! NEG Sin(sigma) in [m2/rad-Hz]
      REAL(KIND=JWRB), DIMENSION(NFRE),      INTENT(IN)  :: DSII ! freq. bandwidths in [radians]
      REAL(KIND=JWRB), DIMENSION(NFRE),      INTENT(IN)  :: SIG  
      
      REAL(KIND=JWRB), PARAMETER :: FRQMAX  = 10.0_JWRB ! Upper freq. limit to extrap. to 
      REAL, ALLOCATABLE :: IK10Hz(:), SIG10Hz(:)
      REAL, ALLOCATABLE :: SPOSDENS10Hz(:), SNEGDENS10Hz(:), DSII10Hz(:)

      INTEGER(KIND=JWIM) :: NK10Hz
      INTEGER(KIND=JWIM) :: NK, NTH
      INTEGER(KIND=JWIM), DIMENSION(NANG) :: ITHN

      REAL(KIND=JWRB) :: PHIWA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
! ----------------------------------------------------------------------------
!

      IF (LHOOK) CALL DR_HOOK('CALCPHIWA',0,ZHOOK_HANDLE)


      NTH   = NANG  ! NUMBER OF DIRS , SAME AS KL
      NK    = NFRE  ! NUMBER OF FREQS, SAME AS ML

!/ 0) --- Find the number of frequencies required to extend arrays
!/        up to f=10Hz and allocate arrays --------------------------- /
      NK10Hz = CEILING(LOG(FRQMAX/(SIG(1)/ZPI))/LOG(FRATIO))+1
      NK10Hz = MAX(NK,NK10Hz)
!
      ALLOCATE(IK10Hz(NK10Hz))
      IK10Hz = REAL( IRANGE(1,NK10Hz,1) )
!
      ALLOCATE(SPOSDENS10Hz(NK10Hz))
      ALLOCATE(SNEGDENS10Hz(NK10Hz))
      ALLOCATE(DSII10Hz(NK10Hz))
      ALLOCATE(SIG10Hz(NK10Hz))
!
      ITHN   = IRANGE(1,NTH,1)    ! Index vector 1:NTH
!
!/ 1) --- Either extrapolate arrays up to 10Hz or use discrete spectral
!         grid per se. Limit the constraint to the positive part of the
!         wind input only. ---------------------------------------------- /
      IF (NK .LT. NK10Hz) THEN
            SPOSDENS10Hz(1:NK)      = SUM(SPOS,1) * DELTH
            SNEGDENS10Hz(1:NK)      = SUM(SNEG,1) * DELTH
            SIG10Hz                 = SIG(1)*FRATIO**(IK10Hz-1.0_JWRB)
            DSII10Hz                = 0.5_JWRB * SIG10Hz * (FRATIO-1.0_JWRB/FRATIO)
!        The first and last frequency bin:
            DSII10Hz(1)             = 0.5_JWRB * SIG10Hz(1) * (FRATIO-1.0_JWRB)
            DSII10Hz(NK10Hz)        = 0.5_JWRB * SIG10Hz(NK10Hz) * (FRATIO-1.0_JWRB) / FRATIO
!
!        --- Spectral slope for S_IN(F) is proportional to F**(-2) ------ /
            SPOSDENS10Hz(NK+1:NK10Hz)  = SPOSDENS10Hz(NK)  * (SIG10Hz(NK)/SIG10Hz(NK+1:NK10Hz))**2
            SNEGDENS10Hz(NK+1:NK10Hz)  = SNEGDENS10Hz(NK)  * (SIG10Hz(NK)/SIG10Hz(NK+1:NK10Hz))**2
      ELSE
            SIG10Hz          = SIG
            DSII10Hz         = DSII
            SPOSDENS10Hz(1:NK)  = SUM(SPOS,1) * DELTH
            SNEGDENS10Hz(1:NK)  = SUM(SNEG,1) * DELTH
      END IF

!/ 2) --- Calculate PHIWA from the extended arrays ------------------- /
      PHIWA      = G * ROWATER * ( SUM(SPOSDENS10Hz*DSII10Hz) + SUM(SNEGDENS10Hz*DSII10Hz) )

      IF (LHOOK) CALL DR_HOOK('CALCPHIWA',1,ZHOOK_HANDLE)

      END FUNCTION CALCPHIWA
