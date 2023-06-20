! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE HALPHAP_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE HALPHAP_CUF_PARAMETRISE (KIJS, KIJL, WAVNUM, COSWDIF, FL1, HALP, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *HALPHAP* - COMPUTATION OF 1/2 PHILLIPS PARAMETER
    
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *HALPHAP(KIJS, KIJL, WAVNUM, UDIR, FL1, HALP)
    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *WAVNUM* - WAVE NUMBER
    !          *COSWDIF*- COSINE ( WIND SPEED DIRECTION - WAVE DIRECTIONS)
    !          *FL1*    - SPECTRA
    !          *HALP*   - 1/2 PHILLIPS PARAMETER
    
    !     METHOD.
    !     -------
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWFRED, ONLY: FR_D, TH_D, FR5_D, DELTH_D, WETAIL, FRTAIL, DFIM_D, DFIMOFR_D
    USE YOWPARAM, ONLY: NANG_D, NFRE_D, NANG_PARAM
    USE YOWPCONS, ONLY: G_D, ZPI_D, ZPI4GM2_D, EPSMIN
    USE YOWPHYS, ONLY: ALPHAPMAX_D
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE_LOKI_PARAM) :: WAVNUM
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(NANG_LOKI_PARAM) :: COSWDIF
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    REAL(KIND=JWRB), INTENT(OUT), DEVICE :: HALP
    
    INTEGER(KIND=JWIM) :: K, M
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    REAL(KIND=JWRB) :: ZLNFRNFRE
    REAL(KIND=JWRB) :: DELT25, DELT2, DEL2
    REAL(KIND=JWRB) :: TEMP1, TEMP2
    REAL(KIND=JWRB) :: ALPHAP
    REAL(KIND=JWRB) :: XMSS, EM, FM, F1D
    REAL(KIND=JWRB), DIMENSION(NANG_PARAM) :: FLWD
    
    ! ----------------------------------------------------------------------
    
    
    ZLNFRNFRE = LOG(FR_D(NFRE_D))
    
    DELT25 = WETAIL*FR_D(NFRE_D)*DELTH_D
    DELT2 = FRTAIL*DELTH_D
    
    ! Find spectrum in wind direction
    DO M=1,NFRE_D
      DO K=1,NANG_D
        FLWD(K) = FL1(IJ, K, M)*0.5_JWRB + 0.5_JWRB*SIGN(1.0_JWRB, COSWDIF(K))
      END DO
      
      XMSS = 0._JWRB
      TEMP1 = DFIM_D(M)*WAVNUM(IJ, M)**2
      TEMP2 = 0.0_JWRB
      DO K=1,NANG_D
        TEMP2 = TEMP2 + FLWD(K)
      END DO
      XMSS = XMSS + TEMP1*TEMP2
      
      K = 1
      EM = 0._JWRB
      FM = 0._JWRB
      TEMP2 = MAX(FLWD(K), EPSMIN)
      DO K=2,NANG_D
        TEMP2 = TEMP2 + MAX(FLWD(K), EPSMIN)
      END DO
      EM = EM + TEMP2*DFIM_D(M)
      FM = FM + DFIMOFR_D(M)*TEMP2
    END DO
    
    DO K=1,NANG_D
      FLWD(K) = FL1(IJ, K, NFRE_D)*0.5_JWRB + 0.5_JWRB*SIGN(1.0_JWRB, COSWDIF(K))
    END DO
    
    EM = EM + DELT25*TEMP2
    FM = FM + DELT2*TEMP2
    FM = EM / FM
    FM = MAX(FM, FR_D(1))
    
    IF (EM > 0.0_JWRB .and. FM < FR_D(NFRE_D - 2)) THEN
      ALPHAP = XMSS / (ZLNFRNFRE - LOG(FM))
      IF (ALPHAP > ALPHAPMAX_D) THEN
        ! some odd cases, revert to tail value
        F1D = 0.0_JWRB
        DO K=1,NANG_D
          F1D = F1D + FLWD(K)*DELTH_D
        END DO
        ALPHAP = ZPI4GM2_D*FR5_D(NFRE_D)*F1D
      END IF
    ELSE
      F1D = 0.0_JWRB
      DO K=1,NANG_D
        F1D = F1D + FLWD(K)*DELTH_D
      END DO
      ALPHAP = ZPI4GM2_D*FR5_D(NFRE_D)*F1D
    END IF
    
    !     1/2 ALPHAP:
    HALP = 0.5_JWRB*MIN(ALPHAP, ALPHAPMAX_D)
    
    
  END SUBROUTINE HALPHAP_CUF_PARAMETRISE
END MODULE HALPHAP_CUF_PARAMETRISE_MOD
