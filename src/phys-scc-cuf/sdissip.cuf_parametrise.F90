! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDISSIP_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE SDISSIP_CUF_PARAMETRISE (KIJS, KIJL, FL1, FLD, SL, INDEP, WAVNUM, XK2CG, EMEAN, F1MEAN, XKMEAN,  &
  & UFRIC, COSWDIF, RAORW, IJ)
    ! ----------------------------------------------------------------------
    
    !**** *SDISSIP* - COMPUTATION OF DEEP WATER DISSIPATION SOURCE FUNCTION.
    
    
    !*    PURPOSE.
    !     --------
    !       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
    !       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
    !       OF DISSIPATION SOURCE FUNCTION.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *SDISSIP (KIJS, KIJL, FL1, FLD, SL, *
    !                        INDEP, WAVNUM, XK2CG,
    !                        EMEAN, F1MEAN, XKMEAN,*
    !                        UFRIC, COSWDIF, RAORW)*
    !         *KIJS* - INDEX OF FIRST GRIDPOINT
    !         *KIJL* - INDEX OF LAST GRIDPOINT
    !         *FL1*  - SPECTRUM.
    !         *FLD*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*  - TOTAL SOURCE FUNCTION ARRAY
    !       *INDEP*  - DEPTH INDEX
    !       *WAVNUM* - WAVE NUMBER
    !       *XK2CG*  - (WAVNUM)**2 * GROUP SPEED
    !        *EMEAN* - MEAN ENERGY DENSITY
    !       *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
    !       *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.
    !       *UFRIC*  - FRICTION VELOCITY IN M/S.
    !       *RAORW*  - RATIO AIR DENSITY TO WATER DENSITY
    !       *COSWDIF*-  COS(TH(K)-WDWAVE(IJ))
    
    ! ----------------------------------------------------------------------
    USE SDISSIP_JAN_CUF_PARAMETRISE_MOD, ONLY: SDISSIP_JAN_CUF_PARAMETRISE
    USE SDISSIP_ARD_CUF_PARAMETRISE_MOD, ONLY: SDISSIP_ARD_CUF_PARAMETRISE
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWPARAM, ONLY: NANG_D, NFRE_D
    USE YOWSTAT, ONLY: IPHYS_D
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    REAL(KIND=JWRB), INTENT(INOUT), DEVICE, DIMENSION(NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FLD, SL
    INTEGER(KIND=JWIM), INTENT(IN), DEVICE, DIMENSION(KIJL) :: INDEP
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NFRE_LOKI_PARAM) :: WAVNUM, XK2CG
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: EMEAN, F1MEAN, XKMEAN
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: UFRIC
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: RAORW
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(NANG_LOKI_PARAM) :: COSWDIF
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    
    ! ----------------------------------------------------------------------
    
    
    SELECT CASE (IPHYS_D)
    CASE (0)
      CALL SDISSIP_JAN_CUF_PARAMETRISE(KIJS, KIJL, FL1, FLD, SL, WAVNUM, EMEAN, F1MEAN, XKMEAN, IJ)
      
    CASE (1)
      CALL SDISSIP_ARD_CUF_PARAMETRISE(KIJS, KIJL, FL1, FLD, SL, INDEP, WAVNUM, XK2CG, UFRIC, COSWDIF, RAORW, IJ)
    END SELECT
    
    
  END SUBROUTINE SDISSIP_CUF_PARAMETRISE
END MODULE SDISSIP_CUF_PARAMETRISE_MOD
