! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDISSIP_ARD_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE SDISSIP_ARD_CUF_PARAMETRISE (KIJS, KIJL, FL1, FLD, SL, INDEP, WAVNUM, XK2CG, UFRIC, COSWDIF,  &
  & RAORW, IJ)
    ! ----------------------------------------------------------------------
    
    !**** *SDISSIP_ARD* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.
    
    !     LOTFI AOUF       METEO FRANCE 2013
    !     FABRICE ARDHUIN  IFREMER  2013
    
    
    !*    PURPOSE.
    !     --------
    !       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
    !       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
    !       OF DISSIPATION SOURCE FUNCTION.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *SDISSIP_ARD (KIJS, KIJL, FL1, FLD,SL,*
    !                            INDEP, WAVNUM, XK2CG,
    !                            UFRIC, COSWDIF, RAORW)*
    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *FL1*    - SPECTRUM.
    !          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*     - TOTAL SOURCE FUNCTION ARRAY
    !          *INDEP*  - DEPTH INDEX
    !          *WAVNUM* - WAVE NUMBER
    !          *XK2CG*  - (WAVE NUMBER)**2 * GROUP SPEED
    !          *UFRIC*  - FRICTION VELOCITY IN M/S.
    !          *RAORW*  - RATIO AIR DENSITY TO WATER DENSITY
    !          *COSWDIF*-  COS(TH(K)-WDWAVE(IJ))
    
    
    !     METHOD.
    !     -------
    
    !       SEE REFERENCES.
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    !       ARDHUIN et AL. JPO DOI:10.1175/20110JPO4324.1
    
    
    ! ----------------------------------------------------------------------
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWFRED, ONLY: FR_D, TH_D, ZPIFR_D
    USE YOWPCONS, ONLY: G_D, ZPI_D
    USE YOWPARAM, ONLY: NANG_D, NFRE_D, NANG_PARAM
    USE YOWPHYS, ONLY: SDSBR, ISDSDTH, ISB, IPSAT, SSDSC2, SSDSC4, SSDSC6, MICHE, SSDSC3, SSDSBRF1, BRKPBCOEF, SSDSC5,  &
    & NSDSNTH_D, NDIKCUMUL_D, INDICESSAT_D, SATWEIGHTS_D, CUMULW_D
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 24
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    REAL(KIND=JWRB), INTENT(INOUT), DEVICE, DIMENSION(NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FLD, SL
    INTEGER(KIND=JWIM), INTENT(IN), DEVICE, DIMENSION(KIJL) :: INDEP
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NFRE_LOKI_PARAM) :: WAVNUM, XK2CG
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: UFRIC
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: RAORW
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(NANG_LOKI_PARAM) :: COSWDIF
    
    
    INTEGER(KIND=JWIM) :: K, M, I, J, M2, K2, KK
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    REAL(KIND=JWRB) :: TPIINV, TPIINVH, TMP01, TMP03
    REAL(KIND=JWRB) :: EPSR, SSDSC6M1, ZCOEF, ZCOEFM1
    
    
    REAL(KIND=JWRB) :: SSDSC2_SIG
    REAL(KIND=JWRB) :: FACTURB, BTH, BTH0
    REAL(KIND=JWRB), DIMENSION(NANG_PARAM) :: SCUMUL, D
    
    REAL(KIND=JWRB) :: RENEWALFREQ
    
    ! ----------------------------------------------------------------------
    
    
    ! INITIALISATION
    
    EPSR = SQRT(SDSBR)
    
    TPIINV = 1.0_JWRB / ZPI_D
    TPIINVH = 0.5_JWRB*TPIINV
    TMP03 = 1.0_JWRB / (SDSBR*MICHE)
    SSDSC6M1 = 1._JWRB - SSDSC6
    
    DO M=1,NFRE_D
      
      ! SATURATION TERM
      SSDSC2_SIG = SSDSC2*ZPIFR_D(M)
      ZCOEF = SSDSC2_SIG*SSDSC6
      ZCOEFM1 = SSDSC2_SIG*SSDSC6M1
      
      ! COMPUTE SATURATION SPECTRUM
      BTH0 = 0.0_JWRB
      
      DO K=1,NANG_D
        BTH = 0.0_JWRB
        ! integrates in directional sector
        DO K2=1,NSDSNTH_D*2 + 1
          KK = INDICESSAT_D(K, K2)
          BTH = BTH + SATWEIGHTS_D(K, K2)*FL1(IJ, KK, M)
        END DO
        BTH = BTH*WAVNUM(IJ, M)*TPIINV*XK2CG(IJ, M)
        BTH0 = MAX(BTH0, BTH)
        
        D(K) = ZCOEFM1*MAX(0._JWRB, BTH*TMP03 - SSDSC4)**IPSAT
        
        SCUMUL(K) = MAX(SQRT(ABS(BTH)) - EPSR, 0._JWRB)**2
      END DO
      
      DO K=1,NANG_D
        ! cumulative term
        D(K) = D(K) + ZCOEF*MAX(0._JWRB, BTH0*TMP03 - SSDSC4)**IPSAT
        IF (BTH0 <= SDSBR) THEN
          SCUMUL(K) = 0._JWRB
        END IF
        
      END DO
      
      IF (M > NDIKCUMUL_D) THEN
        ! CUMULATIVE TERM
        IF (SSDSC3 /= 0.0_JWRB) THEN
          
          DO K=1,NANG_D
            ! Correction of saturation level for shallow-water kinematics
            ! Cumulative effect based on lambda   (breaking probability is
            ! the expected rate of sweeping by larger breaking waves)
            
            RENEWALFREQ = 0.0_JWRB
            
            DO M2=1,M - NDIKCUMUL_D
              DO K2=1,NANG_D
                KK = ABS(K2 - K)
                IF (KK > NANG_D / 2) KK = KK - NANG_D / 2
                ! Integrates over frequencies M2 and directions K2 to
                ! Integration is performed from M2=1 to a frequency lower than M: IK-NDIKCUMUL
                RENEWALFREQ = RENEWALFREQ + CUMULW_D(INDEP(IJ), KK, M2, M)*SCUMUL(K2)
              END DO
            END DO
            
            D(K) = D(K) + RENEWALFREQ
          END DO
        END IF
      END IF
      
      !       WAVE-TURBULENCE INTERACTION TERM
      IF (SSDSC5 /= 0.0_JWRB) THEN
        TMP01 = 2._JWRB*SSDSC5 / G_D
        FACTURB = TMP01*RAORW*UFRIC(IJ)*UFRIC(IJ)
        DO K=1,NANG_D
          D(K) = D(K) - ZPIFR_D(M)*WAVNUM(IJ, M)*FACTURB*COSWDIF(K)
        END DO
      END IF
      
      
      ! ADD ALL CONTRIBUTIONS TO SOURCE TERM
      DO K=1,NANG_D
        SL(K, M) = SL(K, M) + D(K)*FL1(IJ, K, M)
        FLD(K, M) = FLD(K, M) + D(K)
      END DO
    END DO
    
    
  END SUBROUTINE SDISSIP_ARD_CUF_PARAMETRISE
END MODULE SDISSIP_ARD_CUF_PARAMETRISE_MOD
