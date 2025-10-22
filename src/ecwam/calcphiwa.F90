FUNCTION CALCPHIWA(SPOS,SNEG,DF) RESULT(PHIWA)

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
      USE YOWFRED  , ONLY : FR       ,DELTH
      USE YOWPARAM , ONLY : NANG     ,NFRE    
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

!----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN)  :: SPOS ! POS Sin(sigma) in [m2/Hz]
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN)  :: SNEG ! NEG Sin(sigma) in [m2/Hz]
      REAL(KIND=JWRB), DIMENSION(NFRE),      INTENT(IN)  :: DF   ! freq. bandwidths in [Hz]
      
      REAL(KIND=JWRB), DIMENSION(NANG)      :: ZA_SPOS, ZA_SNEG
      REAL(KIND=JWRB), DIMENSION(NANG)      :: SPOS_HF, SNEG_HF, SPOS_LF, SNEG_LF
      REAL(KIND=JWRB), DIMENSION(NANG)      :: SPOS_TOTAL, SNEG_TOTAL

      REAL(KIND=JWRB) :: PHIWA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
! ----------------------------------------------------------------------------
!

      IF (LHOOK) CALL DR_HOOK('CALCPHIWA',0,ZHOOK_HANDLE)

      !/ 0) --- split integral into low/high frequency contributions ------------- /
      !
      !
      !     Th=2pi,f=inf              Th=2pi,f=FR(NFRE)    Th=2pi,f=inf  
      !          / /                    / /                 / /
      !          | | S(f) df dTh  =     | | S(f) df dTh +   | | S(f) df dTh
      !         / /                    / /                 / /
      !     Th=0,f=0                Th=0,f=0            Th=0,f=FR(NFRE)
      !
      !
      !                           =     LF_contribution +  HF_contribution
      !
      !
      !/ 1) --- low frequency contributions to the integral ---------------------- /
      !         -- Direct summation over available freq. bins up to FR(NFRE)
      !
      SPOS_LF    = SUM(SPOS,1) * DELTH * FR
      SNEG_LF    = SUM(SNEG,1) * DELTH * FR

      !/ 2) --- high frequency contributions to the integral --------------------- /
      !         -- Assume spectral slope for S_IN(F) is proportional to F**(-2), then 
      !            integral collapses into easy analytic solution
      !
      !          
      !   Th=2pi,f=inf  
      !     / /
      !     | | S(f) df dTh = ZPI * S(NFRE) / FR(NFRE)
      !    / /
      !  Th=0,f=FR(NFRE)
      !
      !
      ! Determine value of spectrum at NFRE (i.e. at highest frequency). 
      !  - Note, direction dimension must remain
      ZA_SPOS = SPOS(:,NFRE)
      ZA_SNEG = SNEG(:,NFRE)

      SPOS_HF    = (ZPI*ZA_SPOS)/FR(NFRE)
      SNEG_HF    = (ZPI*ZA_SNEG)/FR(NFRE)

      !/ 3) --- summate low + high frequency contributions to the integral ------- /
      SPOS_TOTAL = SPOS_LF + SPOS_HF
      SNEG_TOTAL = SNEG_LF + SNEG_HF

      !/ 4) --- compute the flux ------------------------------------------------- /      
      PHIWA      = G * ROWATER * ( SUM(SPOS_TOTAL) +  SUM(SNEG_TOTAL) )

      IF (LHOOK) CALL DR_HOOK('CALCPHIWA',1,ZHOOK_HANDLE)

      END FUNCTION CALCPHIWA
