! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE Z0WAVE_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE Z0WAVE_CUF_PARAMETRISE (KIJS, KIJL, US, TAUW, UTOP, Z0, Z0B, CHRNCK, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *Z0WAVE* - DETERMINE THE SEA STATE DEPENDENT ROUGHNESS LENGTH.
    
    !*    PURPOSE.
    !     --------
    
    !       COMPUTE ROUGHNESS LENGTH.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *Z0WAVE (KIJS, KIJL, US, TAUW, UTOP, Z0, Z0B, CHRNCK)
    !          *KIJS* - INDEX OF FIRST GRIDPOINT.
    !          *KIJL* - INDEX OF LAST GRIDPOINT.
    !          *US*   - OUTPUT BLOCK OF SURFACE STRESSES.
    !          *TAUW* - INPUT BLOCK OF WAVE STRESSES.
    !          *UTOP* - WIND SPEED.
    !          *Z0*   - OUTPUT BLOCK OF ROUGHNESS LENGTH.
    !          *Z0B*  - BACKGROUND ROUGHNESS LENGTH.
    !          *CHRNCK- CHARNOCK COEFFICIENT
    
    !     METHOD.
    !     -------
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ---------
    
    !       NONE.
    
    ! ----------------------------------------------------------------------
    
    USE CHNKMIN_CUF_PARAMETRISE_MOD, ONLY: CHNKMIN_CUF_PARAMETRISE
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWCOUP, ONLY: LLCAPCHNK_D
    USE YOWPCONS, ONLY: G_D, GM1_D
    USE YOWPHYS, ONLY: ALPHA_D
    USE YOWTABL, ONLY: EPS1
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: US, TAUW, UTOP
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: Z0, Z0B, CHRNCK
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    REAL(KIND=JWRB) :: UST2, UST3, ARG
    REAL(KIND=JWRB) :: ALPHAOG
    
    ! ----------------------------------------------------------------------
    
    
    IF (LLCAPCHNK_D) THEN
      ALPHAOG = CHNKMIN_CUF_PARAMETRISE(UTOP(IJ))*GM1_D
    ELSE
      ALPHAOG = ALPHA_D*GM1_D
    END IF
    
    UST2 = US(IJ)**2
    UST3 = US(IJ)**3
    ARG = MAX(UST2 - TAUW(IJ), EPS1)
    Z0(IJ) = ALPHAOG*UST3 / SQRT(ARG)
    Z0B(IJ) = ALPHAOG*UST2
    CHRNCK(IJ) = G_D*Z0(IJ) / UST2
    
    
  END SUBROUTINE Z0WAVE_CUF_PARAMETRISE
END MODULE Z0WAVE_CUF_PARAMETRISE_MOD
