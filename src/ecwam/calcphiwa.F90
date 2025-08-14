FUNCTION CALCPHIWA(SPOS,SNEG,DSII) RESULT(PHIWA)

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
        USE YOWPCONS , ONLY : G        ,ROWATER
        USE YOWFRED  , ONLY : DELTH
        USE YOWPARAM , ONLY : NANG     ,NFRE    
        USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
  
  !----------------------------------------------------------------------
  
        IMPLICIT NONE
  
        REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN)  :: SPOS ! POS Sin(sigma) in [m2/rad-Hz]
        REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN)  :: SNEG ! NEG Sin(sigma) in [m2/rad-Hz]
        REAL(KIND=JWRB), DIMENSION(NFRE),      INTENT(IN)  :: DSII ! freq. bandwidths in [radians]
  
        REAL(KIND=JWRB), DIMENSION(NFRE) :: SPOSDENSIG, SNEGDENSIG
  
        REAL(KIND=JWRB) :: PHIWA
        REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
  ! ----------------------------------------------------------------------------
  !
  
        IF (LHOOK) CALL DR_HOOK('CALCPHIWA',0,ZHOOK_HANDLE)
  
        SPOSDENSIG = SUM(SPOS,1) * DELTH
        SNEGDENSIG = SUM(SNEG,1) * DELTH
        PHIWA      = G * ROWATER * ( SUM(SPOSDENSIG*DSII) + SUM(SNEGDENSIG*DSII) )
  
        IF (LHOOK) CALL DR_HOOK('CALCPHIWA',1,ZHOOK_HANDLE)
  
        END FUNCTION CALCPHIWA
  