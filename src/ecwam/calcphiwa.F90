FUNCTION CALCPHIWA(SPOS,SNEG) RESULT(PHIWA)

! ----------------------------------------------------------------------------
!
!  1. Purpose :
!
!      Calculate energy flux from wind into waves, obtained from wind-energy-input (Sin).
!
!                              / FRMAX
!      phiwa = g * rho_water * | Sin(f) df
!                              /

!----------------------------------------------------------------------
!
!     INTERFACE VARIABLES.
!     --------------------

!     ORIGIN.
!     ----------
!     Adapted from TAUWINDSXY
!     Implementation into ECWAM DECEMBER 2021 by J. Kousal 

! ----------------------------------------------------------------------------
!

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWPCONS , ONLY : G        ,ROWATER, ZPI
      USE YOWFRED  , ONLY : FR       ,DELTH  , DFIM
      USE YOWPARAM , ONLY : NANG     ,NFRE    
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

!----------------------------------------------------------------------

      IMPLICIT NONE

REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN)  :: SPOS ! POS Sin(sigma) in [m2/Hz]
REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN)  :: SNEG ! NEG Sin(sigma) in [m2/Hz]


REAL(KIND=JWRB), DIMENSION(NANG)      :: ZA_SPOS, ZA_SNEG

REAL(KIND=JWRB)      :: SPOS_LF, SNEG_LF
REAL(KIND=JWRB)      :: SPOS_HF, SNEG_HF
REAL(KIND=JWRB)      :: PHIWA_HF, PHIWA_LF

REAL(KIND=JWRB) :: PHIWA
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
! ----------------------------------------------------------------------------
!

IF (LHOOK) CALL DR_HOOK('CALCPHIWA',0,ZHOOK_HANDLE)

!/ 0) --- split integral into low/high frequency contributions ------------- /
!
!
!     Th=2pi,f=inf                 Th=2pi,f=FR(NFRE)       Th=2pi,f=inf  
!          / /                       / /                    / /
!          | | S(f,Th) df dTh  =     | | S(f,Th) df dTh +   | | S(f,Th) df dTh
!         / /                       / /                    / /
!     Th=0,f=0                   Th=0,f=0               Th=0,f=FR(NFRE)
!
!
!                              =     LF_contribution    +  HF_contribution
!
!
!/ 1) --- low frequency contributions to the integral ---------------------- /
!         -- Direct summation over available freq. bins up to FR(NFRE)
!
!
SPOS_LF    = SUM(SUM(SPOS,1) * DFIM)
SNEG_LF    = SUM(SUM(SNEG,1) * DFIM)

PHIWA_LF   = G * ROWATER * ( SPOS_LF +  SNEG_LF )

!/ 2) --- high frequency contributions to the integral --------------------- /
!         -- Assume spectral slope for S_IN(F) is proportional to F**(-2), then 
!            integral collapses into easy analytic solution
!
!          
!   Th=2pi,f=inf  
!     / /
!     | | S(f,Th) df dTh = FR(NFRE) * DELTH * SUM(S(:,NFRE))
!    / /
!  Th=0,f=FR(NFRE)
!
!
! Determine value of spectrum at NFRE (i.e. at highest frequency). 
!  - Note, direction dimension must remain
ZA_SPOS  = SPOS(:,NFRE)
ZA_SNEG  = SNEG(:,NFRE)

SPOS_HF  = FR(NFRE) * DELTH * SUM(ZA_SPOS)
SNEG_HF  = FR(NFRE) * DELTH * SUM(ZA_SNEG)

PHIWA_HF = G * ROWATER * ( SPOS_HF +  SNEG_HF )

!/ 3) --- summate low + high frequency contributions to the integral ------- /
PHIWA = PHIWA_LF + PHIWA_HF

IF (LHOOK) CALL DR_HOOK('CALCPHIWA',1,ZHOOK_HANDLE)

END FUNCTION CALCPHIWA
