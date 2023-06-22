! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE STRESS_GC_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) FUNCTION STRESS_GC_CUF_PARAMETRISE (ANG_GC, USTAR, Z0, Z0MIN, HALP, RNFAC)
    
    !***  DETERMINE WAVE INDUCED STRESS FOR GRAV-CAP WAVES
    
    !     AUTHOR: PETER JANSSEN
    !     ------
    
    !     REFERENCES:
    !     ----------
    
    !     VIERS PAPER EQ.(29)
    !     FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
    
    !----------------------------------------------------------------------
    
    USE NS_GC_CUF_PARAMETRISE_MOD, ONLY: NS_GC_CUF_PARAMETRISE
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWCOUP, ONLY: LLNORMAGAM_D
    USE YOWFRED, ONLY: NWAV_GC_D, OMEGA_GC_D, XK_GC_D, OMXKM3_GC_D, CM_GC_D, C2OSQRTVG_GC_D, XKMSQRTVGOC2_GC_D, OM3GMKM_GC_D,  &
    & DELKCC_GC_NS_D, DELKCC_OMXKM3_GC_D
    USE YOWPCONS, ONLY: G_D, EPSUS, SURFT
    USE YOWPHYS, ONLY: XKAPPA, ZALP_D, BETAMAXOXKAPPA2_D, BMAXOKAP_D, RN1_RN_D
    
    
    
    !----------------------------------------------------------------------
    
    IMPLICIT NONE
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 24
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    REAL(KIND=JWRB) :: STRESS_GC_CUF_PARAMETRISE
    REAL(KIND=JWRB), INTENT(IN) :: ANG_GC    ! factor to account for angular spreading of the input.
    REAL(KIND=JWRB), INTENT(IN) :: USTAR    ! friction velocity
    REAL(KIND=JWRB), INTENT(IN) :: Z0    !  surface roughness
    REAL(KIND=JWRB), INTENT(IN) :: Z0MIN    ! minimum surface roughness
    REAL(KIND=JWRB), INTENT(IN) :: HALP    ! 1/2 Phillips parameter
    REAL(KIND=JWRB), INTENT(IN) :: RNFAC    ! wind dependent factor used in the growth renormalisation
    
    
    INTEGER(KIND=JWIM) :: NS
    INTEGER(KIND=JWIM) :: I
    
    REAL(KIND=JWRB) :: XLAMBDA    ! Correction factor in the wave growth for gravity-capillary waves
    ! XLAMBDA = 1.0_JWRB + XLAMA * TANH(XLAMB * USTAR**NLAM)
    REAL(KIND=JWRB), PARAMETER :: XLAMA = 0.25_JWRB
    REAL(KIND=JWRB), PARAMETER :: XLAMB = 4.0_JWRB
    INTEGER(KIND=JWIM), PARAMETER :: NLAM = 4
    
!$loki routine seq
    REAL(KIND=JWRB) :: TAUWCG_MIN
    REAL(KIND=JWRB) :: TAUWCG
    REAL(KIND=JWRB) :: ZABHRC
    REAL(KIND=JWRB) :: X, XLOG, ZLOG, ZLOG2X
    REAL(KIND=JWRB) :: CONST, ZN
    REAL(KIND=JWRB) :: GAMNORMA    ! RENORMALISATION FACTOR OF THE GROWTH RATE
    REAL(KIND=JWRB) :: GAM_W
    
    !     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
    
    ! ----------------------------------------------------------------------
    
    
    !*    1.0  DETERMINE GRAV_CAP SPECTRUM, TAUWHF.
    !          ------------------------------------
    
    !     FIND NS:
    NS = NS_GC_CUF_PARAMETRISE(USTAR)
    
    TAUWCG_MIN = (USTAR*(Z0MIN / Z0))**2
    
    XLAMBDA = 1.0_JWRB + XLAMA*TANH(XLAMB*USTAR**NLAM)
    
    ZABHRC = ANG_GC*BETAMAXOXKAPPA2_D*HALP*C2OSQRTVG_GC_D(NS)
    IF (LLNORMAGAM_D) THEN
      CONST = RNFAC*BMAXOKAP_D*HALP*C2OSQRTVG_GC_D(NS) / MAX(USTAR, EPSUS)
    ELSE
      CONST = 0.0_JWRB
    END IF
    
    DO I=NS,NWAV_GC_D
      !       GROWTHRATE BY WIND WITHOUT the multiplicative factor representing the ratio of air density to water density (eps)
      !       and BETAMAXOXKAPPA2
      X = USTAR*CM_GC_D(I)
      XLOG = LOG(XK_GC_D(I)*Z0) + XKAPPA / (X + ZALP_D)
      ZLOG = XLOG - LOG(XLAMBDA)
      ZLOG = MIN(ZLOG, 0.0_JWRB)
      ZLOG2X = ZLOG*ZLOG*X
    END DO
    
    GAM_W = ZLOG2X*ZLOG2X*EXP(XLOG)*OM3GMKM_GC_D(NS)
    ZN = CONST*XKMSQRTVGOC2_GC_D(NS)*GAM_W
    GAMNORMA = (1.0_JWRB + RN1_RN_D*ZN) / (1.0_JWRB + ZN)
    TAUWCG = GAM_W*DELKCC_GC_NS_D(NS)*OMXKM3_GC_D(NS)*GAMNORMA
    DO I=NS + 1,NWAV_GC_D
      !       ANALYTICAL FORM INERTIAL SUB RANGE F(k) = k**(-4)*BB
      !       BB = HALP * C2OSQRTVG_GC(NS)*SQRT(VG_GC(I))/C_GC(I)**2
      !       Tauwcg : (rhow * g /rhoa) * integral of (1/c) * gammma * F(k)  k dk
      !       with omega=g*k and omega=k*c,  then
      !       Tauwcg : (rhow /rhoa) * integral of omega * gammma * F(k)  k dk
      !       but gamma is computed wihtout the rhoa/rhow factor so
      !       Tauwcg : integral of omega * gammma_wam * F(k)  k dk
      !       It should be done in vector form with actual directional spreading information
      !       It simplified here by using the ANG_GC factor.
      GAM_W = ZLOG2X*ZLOG2X*EXP(XLOG)*OM3GMKM_GC_D(I)
      ZN = CONST*XKMSQRTVGOC2_GC_D(I)*GAM_W
      GAMNORMA = (1.0_JWRB + RN1_RN_D*ZN) / (1.0_JWRB + ZN)
      TAUWCG = TAUWCG + GAM_W*DELKCC_OMXKM3_GC_D(I)*GAMNORMA
    END DO
    STRESS_GC_CUF_PARAMETRISE = MAX(ZABHRC*TAUWCG, TAUWCG_MIN)
    
    
  END FUNCTION STRESS_GC_CUF_PARAMETRISE
END MODULE STRESS_GC_CUF_PARAMETRISE_MOD
