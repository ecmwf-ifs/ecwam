! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE STOKESTRN_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE STOKESTRN_CUF_PARAMETRISE (KIJS, KIJL, FL1, WAVNUM, STOKFAC, DEPTH, WSWAVE, WDWAVE, CICOVER,  &
  & CITHICK, USTOKES, VSTOKES, STRNMS, NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *STOKESTRN* - WRAPPER TO CALL STOKESDRIFT and CIMSSTRN
    
    !*    PURPOSE.
    !     --------
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *STOKESTRN (KIJS, KIJL, FL1, WAVNUM, STOKFAC, DEPTH, FF_NOW, INTFLDS, WAM2NEMO)
    
    !          *KIJS*    - INDEX OF FIRST GRIDPOINT.
    !          *KIJL*    - INDEX OF LAST GRIDPOINT.
    !          *FL1*     - SPECTRUM(INPUT).
    !          *WAVNUM*  - WAVE NUMBER.
    !          *STOKFAC* - STOKES DRIFT FACTOR.
    !          *DEPTH*   - WATER DEPTH.
    !          *FF_NOW*  - FORCING FIELDS AT CURRENT TIME.
    !          *INTFLDS* - INTEGRATED/DERIVED PARAMETERS
    !          *WAM2NEMO*- WAVE FIELDS PASSED TO NEMO
    
    ! ----------------------------------------------------------------------
    
    USE STOKESDRIFT_CUF_PARAMETRISE_MOD, ONLY: STOKESDRIFT_CUF_PARAMETRISE
    USE CIMSSTRN_CUF_PARAMETRISE_MOD, ONLY: CIMSSTRN_CUF_PARAMETRISE
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU, JWRO
    USE YOWDRVTYPE, ONLY: FORCING_FIELDS, INTGT_PARAM_FIELDS, WAVE2OCEAN
    
    USE YOWCOUP, ONLY: LWCOU_D, LWNEMOCOU_D, LWNEMOCOUSEND_D, LWNEMOCOUSTK_D, LWNEMOCOUSTRN_D
    USE YOWPARAM, ONLY: NANG_D, NFRE_D
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NFRE_LOKI_PARAM) :: WAVNUM, STOKFAC
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: DEPTH
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: WSWAVE, WDWAVE, CICOVER, CITHICK
    REAL(KIND=JWRB), INTENT(INOUT), DEVICE, DIMENSION(KIJL) :: USTOKES, VSTOKES, STRNMS
    REAL(KIND=JWRO), INTENT(INOUT), DEVICE, DIMENSION(KIJL) :: NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN
    
    
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    ! ----------------------------------------------------------------------
    
    
    CALL STOKESDRIFT_CUF_PARAMETRISE(KIJS, KIJL, FL1, STOKFAC, WSWAVE, WDWAVE, CICOVER, USTOKES, VSTOKES, IJ)
    
    IF (LWNEMOCOUSTRN_D) CALL CIMSSTRN_CUF_PARAMETRISE(KIJS, KIJL, FL1, WAVNUM, DEPTH, CITHICK, STRNMS, IJ)
    
    
    IF (LWNEMOCOU_D .and. (LWNEMOCOUSEND_D .and. LWCOU_D .or. .not.LWCOU_D)) THEN
      IF (LWNEMOCOUSTK_D) THEN
        NEMOUSTOKES(IJ) = USTOKES(IJ)
        NEMOVSTOKES(IJ) = VSTOKES(IJ)
      ELSE
        NEMOUSTOKES(IJ) = 0.0_JWRO
        NEMOVSTOKES(IJ) = 0.0_JWRO
      END IF
      
      IF (LWNEMOCOUSTRN_D) NEMOSTRN(IJ) = STRNMS(IJ)
    END IF
    
    
    ! ----------------------------------------------------------------------
    
  END SUBROUTINE STOKESTRN_CUF_PARAMETRISE
END MODULE STOKESTRN_CUF_PARAMETRISE_MOD
