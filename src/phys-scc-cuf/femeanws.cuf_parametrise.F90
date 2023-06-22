! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE FEMEANWS_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE FEMEANWS_CUF_PARAMETRISE (KIJS, KIJL, FL1, XLLWS, FM, IJ, EM)
    
    ! ----------------------------------------------------------------------
    
    !**** *FEMEANWS* - COMPUTATION OF MEAN ENERGY, MEAN FREQUENCY
    !                  FOR WINDSEA PART OF THE SPECTRUM AS DETERMINED
    !                  BY XLLWS
    
    !*    PURPOSE.
    !     --------
    
    !       COMPUTE MEAN FREQUENCY AT EACH GRID POINT FOR PART OF THE
    !       SPECTRUM WHERE XLLWS IS NON ZERO.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *FEMEANWS (KIJS, KIJL, FL1, XLLWS, EM, FM)*
    !              *KIJS*   - INDEX OF FIRST GRIDPOINT
    !              *KIJL*   - INDEX OF LAST GRIDPOINT
    !              *FL1*    - SPECTRUM.
    !              *XLLWS* - TOTAL WINDSEA MASK FROM INPUT SOURCE TERM
    !              *EM*     - MEAN WAVE ENERGY (OUTPUT)
    !              *FM*     - MEAN WAVE FREQUENCY (OUTPUT)
    
    !     METHOD.
    !     -------
    
    !       NONE.
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    !       NONE.
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWFRED, ONLY: FR_D, DFIM_D, DFIMOFR_D, DELTH_D, WETAIL, FRTAIL
    USE YOWPARAM, ONLY: NANG_D, NFRE_D
    USE YOWPCONS, ONLY: EPSMIN
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 24
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1, XLLWS
    REAL(KIND=JWRB), INTENT(OUT), DEVICE :: FM
    REAL(KIND=JWRB), OPTIONAL, INTENT(OUT), DEVICE :: EM
    
    
    INTEGER(KIND=JWIM) :: M, K
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    REAL(KIND=JWRB) :: DELT25, DELT2
    REAL(KIND=JWRB) :: TEMP2, EM_LOC
    
    ! ----------------------------------------------------------------------
    
    
    !*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
    !        ------------------------------------------------
    
    EM_LOC = EPSMIN
    FM = EPSMIN
    
    DELT25 = WETAIL*FR_D(NFRE_D)*DELTH_D
    DELT2 = FRTAIL*DELTH_D
    
    
    !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
    !        ------------------------------------------
    
    DO M=1,NFRE_D
      TEMP2 = 0.0_JWRB
      DO K=1,NANG_D
        TEMP2 = TEMP2 + XLLWS(IJ, K, M)*FL1(IJ, K, M)
      END DO
      EM_LOC = EM_LOC + DFIM_D(M)*TEMP2
      FM = FM + DFIMOFR_D(M)*TEMP2
    END DO
    
    !*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
    !*       NORMALIZE WITH TOTAL ENERGY.
    !        ------------------------------------------
    
    EM_LOC = EM_LOC + DELT25*TEMP2
    FM = FM + DELT2*TEMP2
    FM = EM_LOC / FM
    
    IF (PRESENT(EM)) THEN
      EM = EM_LOC
    END IF
    
    
  END SUBROUTINE FEMEANWS_CUF_PARAMETRISE
END MODULE FEMEANWS_CUF_PARAMETRISE_MOD
