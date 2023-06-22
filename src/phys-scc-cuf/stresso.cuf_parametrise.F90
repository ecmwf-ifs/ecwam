! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE STRESSO_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE STRESSO_CUF_PARAMETRISE (KIJS, KIJL, MIJ, RHOWGDFTH, FL1, SL, SPOS, CINV, WDWAVE, UFRIC, Z0M,  &
  & AIRD, RNFAC, COSWDIF, SINWDIF2, TAUW, TAUWDIR, PHIWA, LLPHIWA, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *STRESSO* - COMPUTATION OF WAVE STRESS.
    
    !     H. GUNTHER      GKSS/ECMWF NOVEMBER   1989 CODE MOVED FROM SINPUT.
    !     P.A.E.M. JANSSEN     KNMI  AUGUST     1990
    !     J. BIDLOT            ECMWF FEBRUARY   1996-97
    !     S. ABDALLA           ECMWF OCTOBER    2001 INTRODUCTION OF VARIABLE
    !                                                AIR DENSITY
    !     P.A.E.M. JANSSEN     ECMWF            2011  ADD FLUX CALULATIONS
    
    !*    PURPOSE.
    !     --------
    
    !       COMPUTE NORMALIZED WAVE STRESS FROM INPUT SOURCE FUNCTION
    
    !**   INTERFACE.
    !     ----------
    
    !        *CALL* *STRESSO (KIJS, KIJL, MIJ, RHOWGDFTH,
    !                         FL1, SL, SPOS,
    !    &                    CINV,
    !    &                    WDWAVE, UFRIC, Z0M, AIRD, RNFAC,
    !    &                    COSWDIF, SINWDIF2,
    !    &                    TAUW, TAUWDIR, PHIWA)*
    !         *KIJS*        - INDEX OF FIRST GRIDPOINT.
    !         *KIJL*        - INDEX OF LAST GRIDPOINT.
    !         *MIJ*         - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
    !         *RHOWGDFTH    - WATER DENSITY * G * DF * DTHETA
    !         *FL1*         - WAVE SPECTRUM.
    !         *SL*          - WIND INPUT SOURCE FUNCTION ARRAY (positive and negative contributions).
    !         *SPOS*        - POSITIVE WIND INPUT SOURCE FUNCTION ARRAY.
    !         *CINV*        - INVERSE PHASE VELOCITY.
    !         *WDWAVE*      - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
    !                         NOTATION (POINTING ANGLE OF WIND VECTOR,
    !                         CLOCKWISE FROM NORTH).
    !         *UFRIC*       - FRICTION VELOCITY IN M/S.
    !         *Z0M*         - ROUGHNESS LENGTH IN M.
    !         *AIRD*        - AIR DENSITY IN KG/M**3.
    !         *RNFAC*       - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
    !         *COSWDIF*     - COS(TH(K)-WDWAVE(IJ))
    !         *SINWDIF2*    - SIN(TH(K)-WDWAVE(IJ))**2
    !         *TAUW*        - KINEMATIC WAVE STRESS IN (M/S)**2
    !         *TAUWDIR*     - KINEMATIC WAVE STRESS DIRECTION
    !         *PHIWA*       - ENERGY FLUX FROM WIND INTO WAVES INTEGRATED
    !                         OVER THE FULL FREQUENCY RANGE.
    !         *LLPHIWA*     - TRUE IF PHIWA NEEDS TO BE COMPUTED
    
    !     METHOD.
    !     -------
    
    !       THE INPUT SOURCE FUNCTION IS INTEGRATED OVER FREQUENCY
    !       AND DIRECTIONS.
    !       BECAUSE ARRAY *SPOS* IS USED, ONLY THE INPUT SOURCE
    !       HAS TO BE STORED IN *SPOS* (CALL FIRST SINPUT, THEN
    !       STRESSO, AND THEN THE REST OF THE SOURCE FUNCTIONS)
    
    !     REFERENCE.
    !     ----------
    !       P. JANSSEN,
    
    ! ----------------------------------------------------------------------
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWCOUP, ONLY: LLGCBZ0_D
    USE YOWFRED, ONLY: FR_D, RHOWG_DFIM_D, DELTH_D, TH_D, COSTH_D, SINTH_D, FR5_D
    USE YOWPARAM, ONLY: NANG_D, NFRE_D
    USE YOWPHYS, ONLY: TAUWSHELTER_D
    USE YOWTABL, ONLY: EPS1
    USE YOWSTAT, ONLY: IPHYS_D
    
    USE TAU_PHI_HF_CUF_PARAMETRISE_MOD, ONLY: TAU_PHI_HF_CUF_PARAMETRISE
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    INTEGER(KIND=JWIM), INTENT(IN), DEVICE :: MIJ(KIJL)
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(NFRE_LOKI_PARAM) :: RHOWGDFTH
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: SL, SPOS
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NFRE_LOKI_PARAM) :: CINV
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: WDWAVE, UFRIC, Z0M, AIRD
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: RNFAC
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(NANG_LOKI_PARAM) :: COSWDIF, SINWDIF2
    REAL(KIND=JWRB), INTENT(OUT), DEVICE, DIMENSION(KIJL) :: TAUW, TAUWDIR
    REAL(KIND=JWRB), INTENT(OUT), DEVICE :: PHIWA
    LOGICAL, VALUE, INTENT(IN) :: LLPHIWA
    
    
    INTEGER(KIND=JWIM) :: M, K, I, J, II
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    REAL(KIND=JWRB) :: TAUTOUS2
    REAL(KIND=JWRB) :: COSW, FCOSW2
    REAL(KIND=JWRB) :: XSTRESS, YSTRESS
    REAL(KIND=JWRB), DEVICE :: TAUHF, PHIHF
    REAL(KIND=JWRB) :: USDIRP
    REAL(KIND=JWRB), DEVICE :: UST
    
    REAL(KIND=JWRB) :: CMRHOWGDFTH
    REAL(KIND=JWRB) :: TAUX, TAUY, TAUPX, TAUPY
    REAL(KIND=JWRB) :: SUMT, SUMX, SUMY
    
    LOGICAL :: LTAUWSHELTER
    
    ! ----------------------------------------------------------------------
    
    
    PHIWA = 0.0_JWRB
    XSTRESS = 0.0_JWRB
    YSTRESS = 0.0_JWRB
    
    !*    CONTRIBUTION TO THE WAVE STRESS FROM THE NEGATIVE PART OF THE WIND INPUT
    !     ------------------------------------------------------------------------
    
    IF (LLPHIWA) THEN
      !     full energy flux due to negative Sinput (SL-SPOS)
      !     we assume that above NFRE, the contibutions can be neglected
      DO M=1,NFRE_D
        DO K=1,NANG_D
          PHIWA = PHIWA + (SL(K, M) - SPOS(K, M))*RHOWG_DFIM_D(M)
        END DO
      END DO
    END IF
    
    !*    CALCULATE LOW-FREQUENCY CONTRIBUTION TO STRESS AND ENERGY FLUX (positive sinput).
    !     ---------------------------------------------------------------------------------
    DO M=1,NFRE_D
      !     THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
      K = 1
      SUMX = SPOS(K, M)*SINTH_D(K)
      SUMY = SPOS(K, M)*COSTH_D(K)
      DO K=2,NANG_D
        SUMX = SUMX + SPOS(K, M)*SINTH_D(K)
        SUMY = SUMY + SPOS(K, M)*COSTH_D(K)
      END DO
      CMRHOWGDFTH = RHOWGDFTH(M)*CINV(IJ, M)
      XSTRESS = XSTRESS + CMRHOWGDFTH*SUMX
      YSTRESS = YSTRESS + CMRHOWGDFTH*SUMY
    END DO
    
    !     TAUW is the kinematic wave stress !
    XSTRESS = XSTRESS / MAX(AIRD(IJ), 1.0_JWRB)
    YSTRESS = YSTRESS / MAX(AIRD(IJ), 1.0_JWRB)
    
    IF (LLPHIWA) THEN
      DO M=1,NFRE_D
        !       THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
        K = 1
        SUMT = SPOS(K, M)
        DO K=2,NANG_D
          SUMT = SUMT + SPOS(K, M)
        END DO
        PHIWA = PHIWA + RHOWGDFTH(M)*SUMT
      END DO
    END IF
    
    !*    CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS and energy flux (positive sinput).
    !     ----------------------------------------------------------------------------------
    
    IF (IPHYS_D == 0 .or. TAUWSHELTER_D == 0.0_JWRB) THEN
      LTAUWSHELTER = .false.
      USDIRP = WDWAVE(IJ)
      UST = UFRIC(IJ)
    ELSE
      LTAUWSHELTER = .true.
      TAUX = UFRIC(IJ)**2*SIN(WDWAVE(IJ))
      TAUY = UFRIC(IJ)**2*COS(WDWAVE(IJ))
      TAUPX = TAUX - TAUWSHELTER_D*XSTRESS
      TAUPY = TAUY - TAUWSHELTER_D*YSTRESS
      USDIRP = ATAN2(TAUPX, TAUPY)
      UST = (TAUPX**2 + TAUPY**2)**0.25_JWRB
    END IF
    
    CALL TAU_PHI_HF_CUF_PARAMETRISE(KIJS, KIJL, MIJ, LTAUWSHELTER, UFRIC, Z0M, FL1, AIRD, RNFAC, COSWDIF, SINWDIF2, UST, TAUHF,  &
    & PHIHF, LLPHIWA, IJ)
    
    XSTRESS = XSTRESS + TAUHF*SIN(USDIRP)
    YSTRESS = YSTRESS + TAUHF*COS(USDIRP)
    TAUW(IJ) = SQRT(XSTRESS**2 + YSTRESS**2)
    TAUW(IJ) = MAX(TAUW(IJ), 0.0_JWRB)
    TAUWDIR(IJ) = ATAN2(XSTRESS, YSTRESS)
    
    IF (.not.LLGCBZ0_D) THEN
      TAUTOUS2 = 1.0_JWRB / (1.0_JWRB + EPS1)
      TAUW(IJ) = MIN(TAUW(IJ), UFRIC(IJ)**2*TAUTOUS2)
    END IF
    
    IF (LLPHIWA) THEN
      PHIWA = PHIWA + PHIHF
    END IF
    
    
  END SUBROUTINE STRESSO_CUF_PARAMETRISE
END MODULE STRESSO_CUF_PARAMETRISE_MOD
