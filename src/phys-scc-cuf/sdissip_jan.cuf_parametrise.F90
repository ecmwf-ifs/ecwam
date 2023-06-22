! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDISSIP_JAN_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE SDISSIP_JAN_CUF_PARAMETRISE (KIJS, KIJL, FL1, FLD, SL, WAVNUM, EMEAN, F1MEAN, XKMEAN, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *SDISSIP_JAN* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.
    
    !     S.D.HASSELMANN.
    !     MODIFIED TO SHALLOW WATER : G. KOMEN , P. JANSSEN
    !     OPTIMIZATION : L. ZAMBRESKY
    !     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
    !     J. BIDLOT   ECMWF  NOVEMBER 2004  REFORMULATION BASED ON XKMEAN
    !                                       AND F1MEAN.
    !                        AUGUST 2020 Added small viscous dissipation term
    
    !*    PURPOSE.
    !     --------
    !       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
    !       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
    !       OF DISSIPATION SOURCE FUNCTION.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *SDISSIP_JAN (KIJS, KIJ, FL1, FLD, SL,
    !                            WAVNUM,
    !                            EMEAN,F1MEAN, XKMEAN,)*
    !          *FL1*    - SPECTRUM.
    !          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*     - TOTAL SOURCE FUNCTION ARRAY
    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *WAVNUM* - WAVE NUMBER
    !          *EMEAN*  - MEAN ENERGY DENSITY
    !          *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
    !          *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.
    
    
    !     METHOD.
    !     -------
    
    !       SEE REFERENCES.
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    !       G.KOMEN, S. HASSELMANN AND K. HASSELMANN, ON THE EXISTENCE
    !          OF A FULLY DEVELOPED WINDSEA SPECTRUM, JGR, 1984.
    
    ! ---------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWFRED, ONLY: FR_D, DELTH_D, DFIM_D, FRATIO
    USE YOWPARAM, ONLY: NANG_D, NFRE_D
    USE YOWPCONS, ONLY: G_D, ZPI_D, ZPI4GM2_D
    USE YOWPHYS, ONLY: CDIS_D, DELTA_SDIS_D, RNU_D, CDISVIS_D
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 24
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    REAL(KIND=JWRB), INTENT(INOUT), DEVICE, DIMENSION(NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FLD, SL
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NFRE_LOKI_PARAM) :: WAVNUM
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: EMEAN, F1MEAN, XKMEAN
    
    
    INTEGER(KIND=JWIM) :: K, M
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    REAL(KIND=JWRB) :: SCDFM, CONSD, CONSS, DELTA_SDISM1, CVIS
    REAL(KIND=JWRB) :: TEMP1, SDS, X
    REAL(KIND=JWRB) :: XK2
    
    ! ----------------------------------------------------------------------
    
    
    !*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
    !*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
    !        --------------------------------------------------------------
    
    DELTA_SDISM1 = 1.0_JWRB - DELTA_SDIS_D
    
    CONSS = CDIS_D*ZPI_D
    SDS = CONSS*F1MEAN*EMEAN**2*XKMEAN**4
    
    DO M=1,NFRE_D
      X = WAVNUM(IJ, M) / XKMEAN
      XK2 = WAVNUM(IJ, M)**2
      
      CVIS = RNU_D*CDISVIS_D
      TEMP1 = SDS*X*(DELTA_SDISM1 + DELTA_SDIS_D*X) + CVIS*XK2
      
      DO K=1,NANG_D
        FLD(K, M) = FLD(K, M) + TEMP1
        SL(K, M) = SL(K, M) + TEMP1*FL1(IJ, K, M)
      END DO
      
    END DO
    
    
  END SUBROUTINE SDISSIP_JAN_CUF_PARAMETRISE
END MODULE SDISSIP_JAN_CUF_PARAMETRISE_MOD
