MODULE YOWWVDISP

  !-----------------------------------------------------------------
  !  MODULE: YOWWVDISP
  !
  !  PURPOSE:
  !    PROVIDES DISPERSION RELATION FUNCTIONS FOR WATER WAVES
  !    INCLUDING SURFACE TENSION EFFECTS.
  !
  !  AVAILABLE FUNCTIONS:
  !    FOMEG_GC : INTRINSIC (ANGULAR) FREQUENCY  (RAD/S)
  !    FVG_GC   : GROUP VELOCITY                (M/S)
  !    FC_GC    : PHASE VELOCITY                (M/S)
  !
  !  NOTE:
  !    - PHYSICAL CONSTANTS (G, SURFT) ARE PASSED AS ARGUMENTS
  !    - DESIGNED FOR USE IN SPECTRAL WAVE MODELS (ECWAM STYLE)
  !-----------------------------------------------------------------

  USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

  IMPLICIT NONE

  PRIVATE

  !-----------------------------------------------------------------
  !  PUBLIC INTERFACE
  !-----------------------------------------------------------------
  PUBLIC :: FOMEG_GC
  PUBLIC :: FVG_GC
  PUBLIC :: FC_GC

CONTAINS

  !-----------------------------------------------------------------
  !  FUNCTION: FOMEG_GC
  !
  !  PURPOSE:
  !    COMPUTES THE INTRINSIC WAVE FREQUENCY GIVEN WAVENUMBER
  !
  !  FORMULATION:
  !    OMEGA = SQRT( G*K + SURFT*K^3 )
  !
  !  INPUT:
  !    XWNB   : WAVENUMBER (1/M)
  !    G      : GRAVITY (M/S^2)
  !    SURFT  : SURFACE TENSION (N/M)
  !
  !  OUTPUT:
  !    OMEGA  : ANGULAR FREQUENCY (RAD/S)
  !-----------------------------------------------------------------
  FUNCTION FOMEG_GC(XWNB, G, SURFT) RESULT(OMEGA)

    REAL(KIND=JWRB), INTENT(IN) :: XWNB
    REAL(KIND=JWRB), INTENT(IN) :: G
    REAL(KIND=JWRB), INTENT(IN) :: SURFT

    REAL(KIND=JWRB)             :: OMEGA

    OMEGA = SQRT(G * XWNB + SURFT * XWNB**3)

  END FUNCTION FOMEG_GC


  !-----------------------------------------------------------------
  !  FUNCTION: FVG_GC
  !
  !  PURPOSE:
  !    COMPUTES THE GROUP VELOCITY FROM THE DISPERSION RELATION
  !
  !  FORMULATION:
  !    CG = (1 / (2 * OMEGA)) * ( G + 3 * SURFT * K^2 )
  !
  !  INPUT:
  !    XWNB   : WAVENUMBER (1/M)
  !    G      : GRAVITY (M/S^2)
  !    SURFT  : SURFACE TENSION (N/M)
  !
  !  OUTPUT:
  !    CG     : GROUP VELOCITY (M/S)
  !-----------------------------------------------------------------
  FUNCTION FVG_GC(XWNB, G, SURFT) RESULT(CG)

    REAL(KIND=JWRB), INTENT(IN) :: XWNB
    REAL(KIND=JWRB), INTENT(IN) :: G
    REAL(KIND=JWRB), INTENT(IN) :: SURFT

    REAL(KIND=JWRB)             :: CG

    CG = 0.5_JWRB / FOMEG_GC(XWNB, G, SURFT) * (G + 3.0_JWRB * SURFT * XWNB**2)

  END FUNCTION FVG_GC


  !-----------------------------------------------------------------
  !  FUNCTION: FC_GC
  !
  !  PURPOSE:
  !    COMPUTES THE PHASE VELOCITY
  !
  !  FORMULATION:
  !    CP = OMEGA / K
  !
  !  INPUT:
  !    XWNB   : WAVENUMBER (1/M)
  !    G      : GRAVITY (M/S^2)
  !    SURFT  : SURFACE TENSION (N/M)
  !
  !  OUTPUT:
  !    CP     : PHASE VELOCITY (M/S)
  !
  !  NOTE:
  !    - ASSUMES XWNB > 0 (VALID FOR SPECTRAL GRIDS)
  !    - NO EXPLICIT ZERO-DIVISION CHECK FOR PERFORMANCE REASONS
  !-----------------------------------------------------------------
  FUNCTION FC_GC(XWNB, G, SURFT) RESULT(CP)

    REAL(KIND=JWRB), INTENT(IN) :: XWNB
    REAL(KIND=JWRB), INTENT(IN) :: G
    REAL(KIND=JWRB), INTENT(IN) :: SURFT

    REAL(KIND=JWRB)             :: CP

    CP = FOMEG_GC(XWNB, G, SURFT) / XWNB

  END FUNCTION FC_GC

END MODULE YOWWVDISP
