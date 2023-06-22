! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE CIMSSTRN_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE CIMSSTRN_CUF_PARAMETRISE (KIJS, KIJL, FL1, WAVNUM, DEPTH, CITHICK, STRN, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *CIMSSTRN* - COMPUTATION OF THE MEAN SQUARE WAVE STRAIN IN SEA ICE.
    
    !     J. BIDLOT  ECMWF  JANUARY 2013.
    
    !*    PURPOSE.
    !     --------
    
    !       COMPUTES MEAN SQUARE WAVE STRAIN AT EACH GRID POINT.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *CIMSSTRN (KIJS, KIJL, FL1, WAVNUM, DEPTH, CITHICK, STRN)*
    !              *KIJS*    - INDEX OF FIRST GRIDPOINT
    !              *KIJL*    - INDEX OF LAST GRIDPOINT
    !              *FL1*     - SPECTRUM.
    !              *WAVNUM*  - OPEN WATER WAVE NUMBER
    !              *DEPTH*   - WATER DEPTH
    !              *CITHICK* - SEA ICE THICKNESS
    !              *STRN*    - MEAN SQUARE WAVE STRAIN IN ICE (OUTPUT).
    
    !     METHOD.
    !     -------
    
    !      !!! IT ASSUMES SO DEFAULT SETTING FOR THE MECHANICAL PROPERTIES OF
    !          THE SEA ICE (SEE AKI_ICE) !!!!!!!
    
    !       NONE.
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    !       NONE.
    
    ! ----------------------------------------------------------------------
    
    USE AKI_ICE_CUF_PARAMETRISE_MOD, ONLY: AKI_ICE_CUF_PARAMETRISE
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWICE, ONLY: FLMIN
    USE YOWFRED, ONLY: FR_D, DFIM_D, DELTH_D
    USE YOWPARAM, ONLY: NANG_D, NFRE_D
    USE YOWPCONS, ONLY: G_D, ZPI_D, ROWATER
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NFRE_LOKI_PARAM) :: WAVNUM
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: DEPTH
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: CITHICK
    REAL(KIND=JWRB), INTENT(OUT), DEVICE, DIMENSION(KIJL) :: STRN
    
    
    INTEGER(KIND=JWIM) :: M, K
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    REAL(KIND=JWRB) :: F1LIM
    REAL(KIND=JWRB) :: XKI, E, SUME
    
    ! ----------------------------------------------------------------------
    
    
    !*    1. INITIALISE
    !        ----------
    
    F1LIM = FLMIN / DELTH_D
    
    STRN(IJ) = 0.0_JWRB
    
    ! ----------------------------------------------------------------------
    
    !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
    !        ------------------------------------------
    
    DO M=1,NFRE_D
      XKI = AKI_ICE_CUF_PARAMETRISE(G_D, WAVNUM(IJ, M), DEPTH(IJ), ROWATER, CITHICK(IJ))
      E = 0.5_JWRB*CITHICK(IJ)*XKI**3 / WAVNUM(IJ, M)
      
      SUME = 0.0_JWRB
      DO K=1,NANG_D
        SUME = SUME + FL1(IJ, K, M)
      END DO
      
      IF (SUME > F1LIM) THEN
        STRN(IJ) = STRN(IJ) + E**2*SUME*DFIM_D(M)
      END IF
      
    END DO
    
    
  END SUBROUTINE CIMSSTRN_CUF_PARAMETRISE
END MODULE CIMSSTRN_CUF_PARAMETRISE_MOD
