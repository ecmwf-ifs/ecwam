! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE FKMEAN_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE FKMEAN_CUF_PARAMETRISE (KIJS, KIJL, FL1, WAVNUM, EM, FM1, F1, AK, XK, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *FKMEAN* - COMPUTATION OF MEAN FREQUENCIES AT EACH GRID POINT
    !                AND MEAN WAVE NUMBER (based in  sqrt(k)*F moment) .
    !                COMPUTATION OF THE MEAN WAVE ENERGY WAS ALSO
    !                ADDED SUCH THAT A CALL TO FKMEAN DOES NOT NEED
    
    
    !*    PURPOSE.
    !     --------
    
    !       COMPUTE MEAN FREQUENCIES AND WAVE NUMBER AT EACH GRID POINT.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *FKMEAN (KIJS, KIJL, FL1, WAVNUM, EM, FM1, F1, AK, XK)*
    !              *KIJS*    - LOCAL INDEX OF FIRST GRIDPOINT
    !              *KIJL*    - LOCAL INDEX OF LAST GRIDPOINT
    !              *FL1*     - SPECTRUM.
    !              *WAVNUM*  - WAVE NUMBER.
    !              *EM*      - MEAN WAVE ENERGY
    !              *FM1*     - MEAN WAVE FREQUENCY BASED ON (1/f)*FL1 INTEGRATION
    !              *F1*      - MEAN WAVE FREQUENCY BASED ON f*FL1 INTEGRATION
    !              *AK*      - MEAN WAVE NUMBER  BASED ON sqrt(1/k)*FL1 INTGRATION
    !                          ONLY FOR SHALLOW WATER RUNS.
    !!!                        AK IS STILL NEEDED IN SNONLIN !!!!
    !!!                        IF THE OLD FORMULATION IS USED.
    !              *XK*      - MEAN WAVE NUMBER  BASED ON sqrt(k)*FL1 INTEGRATION
    !                          ONLY FOR SHALLOW WATER RUNS.
    
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
    
    USE YOWFRED, ONLY: FR_D, DFIM_D, DFIMOFR_D, DFIMFR_D, DELTH_D, WETAIL, FRTAIL, WP1TAIL
    USE YOWPARAM, ONLY: NANG_D, NFRE_D
    USE YOWPCONS, ONLY: G_D, ZPI_D, EPSMIN
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 24
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NFRE_LOKI_PARAM) :: WAVNUM
    
    REAL(KIND=JWRB), INTENT(OUT), DEVICE :: EM, FM1, F1, AK, XK
    
    
    INTEGER(KIND=JWIM) :: M, K
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    REAL(KIND=JWRB) :: DELT25, COEFM1, COEF1, COEFA, COEFX, SQRTK
    REAL(KIND=JWRB) :: TEMPA, TEMPX, TEMP2
    
    ! ----------------------------------------------------------------------
    
    
    
    !*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
    !        ------------------------------------------------
    
    EM = EPSMIN
    FM1 = EPSMIN
    F1 = EPSMIN
    AK = EPSMIN
    XK = EPSMIN
    
    DELT25 = WETAIL*FR_D(NFRE_D)*DELTH_D
    COEFM1 = FRTAIL*DELTH_D
    COEF1 = WP1TAIL*DELTH_D*FR_D(NFRE_D)**2
    COEFA = COEFM1*SQRT(G_D) / ZPI_D
    COEFX = COEF1*(ZPI_D / SQRT(G_D))
    
    !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
    !        ------------------------------------------
    
    !*    2.2 SHALLOW WATER INTEGRATION.
    !         --------------------------
    
    DO M=1,NFRE_D
      SQRTK = SQRT(WAVNUM(IJ, M))
      TEMPA = DFIM_D(M) / SQRTK
      TEMPX = SQRTK*DFIM_D(M)
      K = 1
      TEMP2 = FL1(IJ, K, M)
      DO K=2,NANG_D
        TEMP2 = TEMP2 + FL1(IJ, K, M)
      END DO
      EM = EM + DFIM_D(M)*TEMP2
      FM1 = FM1 + DFIMOFR_D(M)*TEMP2
      F1 = F1 + DFIMFR_D(M)*TEMP2
      AK = AK + TEMPA*TEMP2
      XK = XK + TEMPX*TEMP2
    END DO
    
    !*      ADD TAIL CORRECTION TO MEAN FREQUENCY AND
    !*      NORMALIZE WITH TOTAL ENERGY.
    EM = EM + DELT25*TEMP2
    FM1 = FM1 + COEFM1*TEMP2
    FM1 = EM / FM1
    F1 = F1 + COEF1*TEMP2
    F1 = F1 / EM
    AK = AK + COEFA*TEMP2
    AK = (EM / AK)**2
    XK = XK + COEFX*TEMP2
    XK = (XK / EM)**2
    
    
  END SUBROUTINE FKMEAN_CUF_PARAMETRISE
END MODULE FKMEAN_CUF_PARAMETRISE_MOD
