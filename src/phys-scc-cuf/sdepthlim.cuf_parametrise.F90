! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDEPTHLIM_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE SDEPTHLIM_CUF_PARAMETRISE (KIJS, KIJL, EMAXDPT, FL1, IJ)
    ! ----------------------------------------------------------------------
    !     J. BIDLOT    ECMWF   NOVEMBER 2017
    
    !*    PURPOSE.
    !     --------
    !     LIMITS THE SPECTRAL VARIANCE SUCH THAT THE TOTAL VARIANCE
    !     DOES NOT EXCEED THE MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
    
    !**   INTERFACE.
    !     ----------
    !     *CALL* *SDEPTHLIM((KIJS, KIJL, EMAXDPT, FL1)
    !          *KIJS*   - LOCAL INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - LOCAL  INDEX OF LAST GRIDPOIN
    !          *EMAXDPT - MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
    !          *FL1*    - SPECTRUM.
    
    
    !     METHOD.
    !     -------
    
    !     EXTERNALS.
    !     ----------
    
    !     REFERENCE.
    !     ----------
    !     NONE
    
    ! ----------------------------------------------------------------------
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWPARAM, ONLY: NANG_D, NFRE_D
    USE YOWPCONS, ONLY: EPSMIN
    USE YOWFRED, ONLY: FR_D, DFIM_D, DELTH_D, WETAIL
    
    
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 24
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: EMAXDPT
    REAL(KIND=JWRB), INTENT(INOUT), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    
    INTEGER(KIND=JWIM) :: K, M
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    REAL(KIND=JWRB) :: DELT25
    REAL(KIND=JWRB) :: EM
    REAL(KIND=JWRB) :: TEMP
    LOGICAL :: LLEPSMIN
    
    ! ----------------------------------------------------------------------
    
    
    EM = EPSMIN
    DO M=1,NFRE_D
      K = 1
      TEMP = FL1(IJ, K, M)
      DO K=2,NANG_D
        TEMP = TEMP + FL1(IJ, K, M)
      END DO
      EM = EM + DFIM_D(M)*TEMP
    END DO
    ! ----------------------------------------------------------------------
    
    !*    3. ADD TAIL ENERGY.
    !        ----------------
    
    DELT25 = WETAIL*FR_D(NFRE_D)*DELTH_D
    EM = EM + DELT25*TEMP
    
    EM = MIN(EMAXDPT(IJ) / EM, 1.0_JWRB)
    
    DO M=1,NFRE_D
      DO K=1,NANG_D
        FL1(IJ, K, M) = MAX(FL1(IJ, K, M)*EM, EPSMIN)
      END DO
    END DO
    
    
  END SUBROUTINE SDEPTHLIM_CUF_PARAMETRISE
END MODULE SDEPTHLIM_CUF_PARAMETRISE_MOD
