! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE FRCUTINDEX_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE FRCUTINDEX_CUF_PARAMETRISE (KIJS, KIJL, FM, FMWS, UFRIC, CICOVER, MIJ, RHOWGDFTH, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *FRCUTINDEX* - RETURNS THE LAST FREQUENCY INDEX OF
    !                    PROGNOSTIC PART OF SPECTRUM.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *FRCUTINDEX (KIJS, KIJL, FM, FMWS, CICOVER, MIJ, RHOWGDFTH)
    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *FM*     - MEAN FREQUENCY
    !          *FMWS*   - MEAN FREQUENCY OF WINDSEA
    !          *UFRIC*  - FRICTION VELOCITY IN M/S
    !          *CICOVER*- CICOVER
    !          *MIJ*    - LAST FREQUENCY INDEX for imposing high frequency tail
    !          *RHOWGDFTH - WATER DENSITY * G * DF * DTHETA
    !                       FOR TRAPEZOIDAL INTEGRATION BETWEEN FR(1) and FR(MIJ)
    !                       !!!!!!!!  RHOWGDFTH=0 FOR FR > FR(MIJ)
    
    
    !     METHOD.
    !     -------
    
    !*    COMPUTES LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
    !*    FREQUENCIES LE 2.5*MAX(FMWS,FM).
    
    
    !!! be aware that if this is NOT used, for iphys=1, the cumulative dissipation has to be
    !!! re-activated (see module yowphys) !!!
    
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWFRED, ONLY: FR_D, DFIM_D, FRATIO, FLOGSPRDM1_D, ZPIFR_D, DELTH_D, RHOWG_DFIM_D, FRIC
    USE YOWICE, ONLY: CITHRSH_TAIL_D
    USE YOWPARAM, ONLY: NFRE_D
    USE YOWPCONS, ONLY: G_D, EPSMIN
    USE YOWPHYS, ONLY: TAILFACTOR_D, TAILFACTOR_PM_D
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    INTEGER(KIND=JWIM), INTENT(OUT), DEVICE :: MIJ(KIJL)
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: FM, FMWS
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: UFRIC, CICOVER
    REAL(KIND=JWRB), INTENT(OUT), DEVICE, DIMENSION(NFRE_LOKI_PARAM) :: RHOWGDFTH
    
    
    INTEGER(KIND=JWIM) :: M
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    REAL(KIND=JWRB) :: FPMH, FPPM, FM2, FPM, FPM4
    
    ! ----------------------------------------------------------------------
    
    
    !*    COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
    !*    FREQUENCIES LE MAX(TAILFACTOR*MAX(FMNWS,FM),TAILFACTOR_PM*FPM),
    !*    WHERE FPM IS THE PIERSON-MOSKOWITZ FREQUENCY BASED ON FRICTION
    !*    VELOCITY. (FPM=G/(FRIC*ZPI*USTAR))
    !     ------------------------------------------------------------
    
    FPMH = TAILFACTOR_D / FR_D(1)
    FPPM = TAILFACTOR_PM_D*G_D / (FRIC*ZPIFR_D(1))
    
    IF (CICOVER(IJ) <= CITHRSH_TAIL_D) THEN
      FM2 = MAX(FMWS, FM)*FPMH
      FPM = FPPM / MAX(UFRIC(IJ), EPSMIN)
      FPM4 = MAX(FM2, FPM)
      MIJ(IJ) = NINT(LOG10(FPM4)*FLOGSPRDM1_D) + 1
      MIJ(IJ) = MIN(MAX(1, MIJ(IJ)), NFRE_D)
    ELSE
      MIJ(IJ) = NFRE_D
    END IF
    
    !     SET RHOWGDFTH
    DO M=1,MIJ(IJ)
      RHOWGDFTH(M) = RHOWG_DFIM_D(M)
    END DO
    IF (MIJ(IJ) /= NFRE_D) RHOWGDFTH(MIJ(IJ)) = 0.5_JWRB*RHOWGDFTH(MIJ(IJ))
    DO M=MIJ(IJ) + 1,NFRE_D
      RHOWGDFTH(M) = 0.0_JWRB
    END DO
    
    
  END SUBROUTINE FRCUTINDEX_CUF_PARAMETRISE
END MODULE FRCUTINDEX_CUF_PARAMETRISE_MOD
