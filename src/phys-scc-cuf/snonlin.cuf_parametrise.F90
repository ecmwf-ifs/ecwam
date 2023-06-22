! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SNONLIN_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE SNONLIN_CUF_PARAMETRISE (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN, ENH, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *SNONLIN* - COMPUTATION OF NONLINEAR TRANSFER RATE AND ITS
    !****             FUNCTIONAL DERIVATIVE (DIAGONAL TERMS ONLY) AND
    !****             ADDITION TO CORRESPONDING NET EXPRESSIONS.
    
    !     S.D. HASSELMANN.  MPI
    
    !     G. KOMEN, P. JANSSEN   KNMI             MODIFIED TO SHALLOW WATER
    !     H. GUENTHER, L. ZAMBRESKY               OPTIMIZED
    !     H. GUENTHER       GKSS/ECMWF  JUNE 1991 INTERACTIONS BETWEEN DIAG-
    !                                             AND PROGNOSTIC PART.
    !     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
    !     P. JANSSEN  ECMWF  JUNE 2005       IMPROVED SCALING IN SHALLOW
    !                                        WATER
    !     J. BIDLOT   ECMWF  AUGUST 2006     KEEP THE OLD FORMULATION
    !                                        UNDER A SWITCH (ISNONLIN = 0 for OLD
    !                                                                 = 1 for NEW
    !                                        BE AWARE THAT THE OLD FORMULATION
    !                                        REQUIRES THE MEAN WAVE NUMBER AKMEAN.
    !     J. BIDLOT   ECMWF  JANUARY 2012    ADD EXTENSION TO LOW FREQUENCIES
    !                                        OPTIMISATION FOR IBM.
    
    !*    PURPOSE.
    !     --------
    
    !       SEE ABOVE.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *SNONLIN (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN)*
    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *FL1*    - SPECTRUM.
    !          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*     - TOTAL SOURCE FUNCTION ARRAY.
    !          *WAVNUM* - WAVE NUMBER.
    !          *DEPTH*  - WATER DEPTH.
    !          *AKMEAN* - MEAN WAVE NUMBER  BASED ON sqrt(1/k)*F INTGRATION
    
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
    USE TRANSF_SNL_CUF_PARAMETRISE_MOD, ONLY: TRANSF_SNL_CUF_PARAMETRISE
    USE TRANSF_CUF_PARAMETRISE_MOD, ONLY: TRANSF_CUF_PARAMETRISE
    USE PEAK_ANG_CUF_PARAMETRISE_MOD, ONLY: PEAK_ANG_CUF_PARAMETRISE
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWFRED, ONLY: FR_D, FRATIO, ZPIFR_D
    USE YOWPARAM, ONLY: NANG_D, NFRE_D
    USE YOWINDN, ONLY: IKP_D, IKP1_D, IKM_D, IKM1_D, K1W_D, K2W_D, K11W_D, K21W_D, AF11_D, FKLAP_D, FKLAP1_D, FKLAM_D, FKLAM1_D,  &
    & DAL1_D, DAL2_D, MFRSTLW_D, MLSTHG_D, KFRH_D, INLCOEF_D, RNLCOEF_D
    USE YOWSTAT, ONLY: ISNONLIN_D
    USE YOWPCONS, ONLY: GM1_D
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 24
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    REAL(KIND=JWRB), INTENT(INOUT), DEVICE, DIMENSION(NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FLD, SL
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NFRE_LOKI_PARAM) :: WAVNUM
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: DEPTH
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: AKMEAN
    REAL(KIND=JWRB), INTENT(INOUT), DEVICE, DIMENSION(KIJL, MLSTHG_D) :: ENH
    
    INTEGER(KIND=JWIM) :: K, M, MC, KH, K1, K2, K11, K21
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: MP, MP1, MM, MM1, IC, IP, IP1, IM, IM1
    INTEGER(KIND=JWIM) :: MFR1STFR, MFRLSTFR
    
    REAL(KIND=JWRB), PARAMETER :: ENH_MAX = 10.0_JWRB
    REAL(KIND=JWRB), PARAMETER :: ENH_MIN = 0.1_JWRB    ! to prevent ENH to become too small
    REAL(KIND=JWRB) :: XK
    
    REAL(KIND=JWRB) :: FTAIL, FKLAMP, GW1, GW2, GW3, GW4
    REAL(KIND=JWRB) :: FKLAMPA, FKLAMPB, FKLAMP2, FKLAMP1, FKLAPA2
    REAL(KIND=JWRB) :: FKLAPB2, FKLAP12, FKLAP22, FKLAMM, FKLAMM1
    REAL(KIND=JWRB) :: GW5, GW6, GW7, GW8, FKLAMMA, FKLAMMB, FKLAMM2
    REAL(KIND=JWRB) :: FKLAMA2, FKLAMB2, FKLAM12, FKLAM22
    REAL(KIND=JWRB) :: SAP, SAM, FIJ, FAD1, FAD2, FCEN
    
    REAL(KIND=JWRB), DEVICE :: XNU, SIG_TH
    REAL(KIND=JWRB) :: FTEMP, AD, DELAD, DELAP, DELAM, ENHFR
    
    ! ----------------------------------------------------------------------
    
    
    
    !*    1. SHALLOW WATER SCALING
    !        ---------------------
    
    SELECT CASE (ISNONLIN_D)
    CASE (0)
      ENHFR = MAX(0.75_JWRB*DEPTH(IJ)*AKMEAN, 0.5_JWRB)
      ENHFR = 1.0_JWRB + (5.5_JWRB / ENHFR)*(1.0_JWRB - .833_JWRB*ENHFR)*EXP(-1.25_JWRB*ENHFR)
      DO MC=1,MLSTHG_D
        ENH(IJ, MC) = ENHFR
      END DO
      
    CASE (1)
      DO MC=1,NFRE_D
        ENH(IJ, MC) = MAX(MIN(ENH_MAX, TRANSF_CUF_PARAMETRISE(WAVNUM(IJ, MC), DEPTH(IJ))), ENH_MIN)
      END DO
      DO MC=NFRE_D + 1,MLSTHG_D
        XK = GM1_D*(ZPIFR_D(NFRE_D)*FRATIO**(MC - NFRE_D))**2
        ENH(IJ, MC) = MAX(MIN(ENH_MAX, TRANSF_CUF_PARAMETRISE(XK, DEPTH(IJ))), ENH_MIN)
      END DO
      
    CASE (2)
      CALL PEAK_ANG_CUF_PARAMETRISE(KIJS, KIJL, FL1, XNU, SIG_TH, IJ)
      DO MC=1,NFRE_D
        ENH(IJ, MC) = TRANSF_SNL_CUF_PARAMETRISE(WAVNUM(IJ, MC), DEPTH(IJ), XNU, SIG_TH)
      END DO
      DO MC=NFRE_D + 1,MLSTHG_D
        XK = GM1_D*(ZPIFR_D(NFRE_D)*FRATIO**(MC - NFRE_D))**2
        ENH(IJ, MC) = TRANSF_SNL_CUF_PARAMETRISE(XK, DEPTH(IJ), XNU, SIG_TH)
      END DO
    END SELECT
    
    
    !*    2. FREQUENCY LOOP.
    !        ---------------
    
    MFR1STFR = -MFRSTLW_D + 1
    MFRLSTFR = NFRE_D - KFRH_D + MFR1STFR
    
    DO MC=1,MLSTHG_D
      MP = IKP_D(MC)
      MP1 = IKP1_D(MC)
      MM = IKM_D(MC)
      MM1 = IKM1_D(MC)
      IC = INLCOEF_D(1, MC)
      IP = INLCOEF_D(2, MC)
      IP1 = INLCOEF_D(3, MC)
      IM = INLCOEF_D(4, MC)
      IM1 = INLCOEF_D(5, MC)
      
      FTAIL = RNLCOEF_D(1, MC)
      
      FKLAMP = FKLAP_D(MC)
      FKLAMP1 = FKLAP1_D(MC)
      GW1 = RNLCOEF_D(2, MC)
      GW2 = RNLCOEF_D(3, MC)
      GW3 = RNLCOEF_D(4, MC)
      GW4 = RNLCOEF_D(5, MC)
      FKLAMPA = RNLCOEF_D(6, MC)
      FKLAMPB = RNLCOEF_D(7, MC)
      FKLAMP2 = RNLCOEF_D(8, MC)
      FKLAMP1 = RNLCOEF_D(9, MC)
      FKLAPA2 = RNLCOEF_D(10, MC)
      FKLAPB2 = RNLCOEF_D(11, MC)
      FKLAP12 = RNLCOEF_D(12, MC)
      FKLAP22 = RNLCOEF_D(13, MC)
      
      FKLAMM = FKLAM_D(MC)
      FKLAMM1 = FKLAM1_D(MC)
      GW5 = RNLCOEF_D(14, MC)
      GW6 = RNLCOEF_D(15, MC)
      GW7 = RNLCOEF_D(16, MC)
      GW8 = RNLCOEF_D(17, MC)
      FKLAMMA = RNLCOEF_D(18, MC)
      FKLAMMB = RNLCOEF_D(19, MC)
      FKLAMM2 = RNLCOEF_D(20, MC)
      FKLAMM1 = RNLCOEF_D(21, MC)
      FKLAMA2 = RNLCOEF_D(22, MC)
      FKLAMB2 = RNLCOEF_D(23, MC)
      FKLAM12 = RNLCOEF_D(24, MC)
      FKLAM22 = RNLCOEF_D(25, MC)
      
      FTEMP = AF11_D(MC)*ENH(IJ, MC)
      
      
      IF (MC > MFR1STFR .and. MC < MFRLSTFR) THEN
        !       the interactions for MC are all within the fully resolved spectral domain
        
        DO KH=1,2
          DO K=1,NANG_D
            K1 = K1W_D(K, KH)
            K2 = K2W_D(K, KH)
            K11 = K11W_D(K, KH)
            K21 = K21W_D(K, KH)
            
            !*    2.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND
            !*            DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
            !             ----------------------------------------------
            SAP = GW1*FL1(IJ, K1, IP) + GW2*FL1(IJ, K11, IP) + GW3*FL1(IJ, K1, IP1) + GW4*FL1(IJ, K11, IP1)
            SAM = GW5*FL1(IJ, K2, IM) + GW6*FL1(IJ, K21, IM) + GW7*FL1(IJ, K2, IM1) + GW8*FL1(IJ, K21, IM1)
            !!!! not needed ftail always=1.                FIJ = FL1(IJ,K  ,IC )*FTAIL
            FIJ = FL1(IJ, K, IC)
            FAD1 = FIJ*(SAP + SAM)
            FAD2 = FAD1 - 2.0_JWRB*SAP*SAM
            FAD1 = FAD1 + FAD2
            FCEN = FTEMP*FIJ
            AD = FAD2*FCEN
            DELAD = FAD1*FTEMP
            DELAP = (FIJ - 2.0_JWRB*SAM)*DAL1_D*FCEN
            DELAM = (FIJ - 2.0_JWRB*SAP)*DAL2_D*FCEN
            
            SL(K, MC) = SL(K, MC) - 2.0_JWRB*AD
            FLD(K, MC) = FLD(K, MC) - 2.0_JWRB*DELAD
            SL(K2, MM) = SL(K2, MM) + AD*FKLAMM1
            FLD(K2, MM) = FLD(K2, MM) + DELAM*FKLAM12
            SL(K21, MM) = SL(K21, MM) + AD*FKLAMM2
            FLD(K21, MM) = FLD(K21, MM) + DELAM*FKLAM22
            SL(K2, MM1) = SL(K2, MM1) + AD*FKLAMMA
            FLD(K2, MM1) = FLD(K2, MM1) + DELAM*FKLAMA2
            SL(K21, MM1) = SL(K21, MM1) + AD*FKLAMMB
            FLD(K21, MM1) = FLD(K21, MM1) + DELAM*FKLAMB2
            SL(K1, MP) = SL(K1, MP) + AD*FKLAMP1
            FLD(K1, MP) = FLD(K1, MP) + DELAP*FKLAP12
            SL(K11, MP) = SL(K11, MP) + AD*FKLAMP2
            FLD(K11, MP) = FLD(K11, MP) + DELAP*FKLAP22
            SL(K1, MP1) = SL(K1, MP1) + AD*FKLAMPA
            FLD(K1, MP1) = FLD(K1, MP1) + DELAP*FKLAPA2
            SL(K11, MP1) = SL(K11, MP1) + AD*FKLAMPB
            FLD(K11, MP1) = FLD(K11, MP1) + DELAP*FKLAPB2
          END DO
        END DO
        
      ELSE IF (MC >= MFRLSTFR) THEN
        DO KH=1,2
          DO K=1,NANG_D
            K1 = K1W_D(K, KH)
            K2 = K2W_D(K, KH)
            K11 = K11W_D(K, KH)
            K21 = K21W_D(K, KH)
            
            SAP = GW1*FL1(IJ, K1, IP) + GW2*FL1(IJ, K11, IP) + GW3*FL1(IJ, K1, IP1) + GW4*FL1(IJ, K11, IP1)
            SAM = GW5*FL1(IJ, K2, IM) + GW6*FL1(IJ, K21, IM) + GW7*FL1(IJ, K2, IM1) + GW8*FL1(IJ, K21, IM1)
            FIJ = FL1(IJ, K, IC)*FTAIL
            FAD1 = FIJ*(SAP + SAM)
            FAD2 = FAD1 - 2.0_JWRB*SAP*SAM
            FAD1 = FAD1 + FAD2
            FCEN = FTEMP*FIJ
            AD = FAD2*FCEN
            DELAD = FAD1*FTEMP
            DELAP = (FIJ - 2.0_JWRB*SAM)*DAL1_D*FCEN
            DELAM = (FIJ - 2.0_JWRB*SAP)*DAL2_D*FCEN
            
            SL(K2, MM) = SL(K2, MM) + AD*FKLAMM1
            FLD(K2, MM) = FLD(K2, MM) + DELAM*FKLAM12
            SL(K21, MM) = SL(K21, MM) + AD*FKLAMM2
            FLD(K21, MM) = FLD(K21, MM) + DELAM*FKLAM22
            
            IF (MM1 <= NFRE_D) THEN
              SL(K2, MM1) = SL(K2, MM1) + AD*FKLAMMA
              FLD(K2, MM1) = FLD(K2, MM1) + DELAM*FKLAMA2
              SL(K21, MM1) = SL(K21, MM1) + AD*FKLAMMB
              FLD(K21, MM1) = FLD(K21, MM1) + DELAM*FKLAMB2
              
              IF (MC <= NFRE_D) THEN
                SL(K, MC) = SL(K, MC) - 2.0_JWRB*AD
                FLD(K, MC) = FLD(K, MC) - 2.0_JWRB*DELAD
                
                IF (MP <= NFRE_D) THEN
                  SL(K1, MP) = SL(K1, MP) + AD*FKLAMP1
                  FLD(K1, MP) = FLD(K1, MP) + DELAP*FKLAP12
                  SL(K11, MP) = SL(K11, MP) + AD*FKLAMP2
                  FLD(K11, MP) = FLD(K11, MP) + DELAP*FKLAP22
                  
                  IF (MP1 <= NFRE_D) THEN
                    SL(K1, MP1) = SL(K1, MP1) + AD*FKLAMPA
                    FLD(K1, MP1) = FLD(K1, MP1) + DELAP*FKLAPA2
                    SL(K11, MP1) = SL(K11, MP1) + AD*FKLAMPB
                    FLD(K11, MP1) = FLD(K11, MP1) + DELAP*FKLAPB2
                  END IF
                END IF
              END IF
            END IF
          END DO
        END DO
        
      ELSE
        
        DO KH=1,2
          DO K=1,NANG_D
            K1 = K1W_D(K, KH)
            K2 = K2W_D(K, KH)
            K11 = K11W_D(K, KH)
            K21 = K21W_D(K, KH)
            
            SAP = GW1*FL1(IJ, K1, IP) + GW2*FL1(IJ, K11, IP) + GW3*FL1(IJ, K1, IP1) + GW4*FL1(IJ, K11, IP1)
            SAM = GW5*FL1(IJ, K2, IM) + GW6*FL1(IJ, K21, IM) + GW7*FL1(IJ, K2, IM1) + GW8*FL1(IJ, K21, IM1)
            FIJ = FL1(IJ, K, IC)*FTAIL
            FAD1 = FIJ*(SAP + SAM)
            FAD2 = FAD1 - 2.0_JWRB*SAP*SAM
            FAD1 = FAD1 + FAD2
            FCEN = FTEMP*FIJ
            AD = FAD2*FCEN
            DELAD = FAD1*FTEMP
            DELAP = (FIJ - 2.0_JWRB*SAM)*DAL1_D*FCEN
            DELAM = (FIJ - 2.0_JWRB*SAP)*DAL2_D*FCEN
            
            IF (MM1 >= 1) THEN
              SL(K2, MM1) = SL(K2, MM1) + AD*FKLAMMA
              FLD(K2, MM1) = FLD(K2, MM1) + DELAM*FKLAMA2
              SL(K21, MM1) = SL(K21, MM1) + AD*FKLAMMB
              FLD(K21, MM1) = FLD(K21, MM1) + DELAM*FKLAMB2
            END IF
            
            SL(K, MC) = SL(K, MC) - 2.0_JWRB*AD
            FLD(K, MC) = FLD(K, MC) - 2.0_JWRB*DELAD
            SL(K1, MP) = SL(K1, MP) + AD*FKLAMP1
            FLD(K1, MP) = FLD(K1, MP) + DELAP*FKLAP12
            SL(K11, MP) = SL(K11, MP) + AD*FKLAMP2
            FLD(K11, MP) = FLD(K11, MP) + DELAP*FKLAP22
            SL(K1, MP1) = SL(K1, MP1) + AD*FKLAMPA
            FLD(K1, MP1) = FLD(K1, MP1) + DELAP*FKLAPA2
            SL(K11, MP1) = SL(K11, MP1) + AD*FKLAMPB
            FLD(K11, MP1) = FLD(K11, MP1) + DELAP*FKLAPB2
          END DO
        END DO
        
      END IF
      
      !*    BRANCH BACK TO 2. FOR NEXT FREQUENCY.
      
    END DO
    
    
  END SUBROUTINE SNONLIN_CUF_PARAMETRISE
END MODULE SNONLIN_CUF_PARAMETRISE_MOD
