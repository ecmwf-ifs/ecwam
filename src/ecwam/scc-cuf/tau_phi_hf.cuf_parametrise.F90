! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE TAU_PHI_HF_CUF_PARAMETRISE_MOD
  !CONTAINED SUBROUTINES:
  ! - OMEGAGC
  ! - TAU_PHI_HF
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE OMEGAGC_CUF_PARAMETRISE (UST, NS, XKS, OMS)
    
    !***  DETERMINE THE CUT-OFF ANGULAR FREQUENCY FOR THE GRAV-CAPILLARY WAVES
    !     !!!! rounded to the closest index of XK_GC  !!!!!
    
    !     AUTHOR: PETER JANSSEN
    !     ------
    
    !     REFERENCES:
    !     ----------
    
    !     VIERS PAPER EQ.(29)
    
    !----------------------------------------------------------------------
    
    USE NS_GC_CUF_PARAMETRISE_MOD, ONLY: NS_GC_CUF_PARAMETRISE
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWFRED, ONLY: OMEGA_GC_D, XK_GC_D
    
    
    !----------------------------------------------------------------------
    
    IMPLICIT NONE
!$loki routine seq
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    REAL(KIND=JWRB), INTENT(IN) :: UST
    INTEGER(KIND=JWIM), INTENT(OUT) :: NS    ! index in array XK_GC corresponding to XKS and OMS
    REAL(KIND=JWRB), INTENT(OUT) :: XKS    ! cut-off wave number
    REAL(KIND=JWRB), INTENT(OUT) :: OMS    ! cut-off angular frequency
    
    
    
    ! ----------------------------------------------------------------------
    
    
    NS = NS_GC_CUF_PARAMETRISE(UST)
    XKS = XK_GC_D(NS)
    OMS = OMEGA_GC_D(NS)
    
    
  END SUBROUTINE OMEGAGC_CUF_PARAMETRISE
  
  ATTRIBUTES(DEVICE) SUBROUTINE TAU_PHI_HF_CUF_PARAMETRISE (KIJS, KIJL, MIJ, LTAUWSHELTER, UFRIC, Z0M, FL1, AIRD, RNFAC,  &
  & COSWDIF, SINWDIF2, UST, TAUHF, PHIHF, LLPHIHF, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *TAU_PHI_HF* - COMPUTATION OF HIGH-FREQUENCY STRESS.
    !                                   HIGH-FREQUENCY ENERGY FLUX.
    
    !     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90
    !     JEAN BIDLOT  ECMWF  JANUARY 2017
    
    !*    PURPOSE.
    !     ---------
    
    !       COMPUTE HIGH-FREQUENCY WAVE STRESS AND ENERGY FLUX
    
    !**   INTERFACE.
    !     ---------
    
    !       *CALL* *TAU_PHI_HF(KIJS, KIJL, MIJ, LTAUWSHELTER, UFRIC, UST, Z0M,
    !                          FL1, AIRD, RNFAC,
    !                          COSWDIF, SINWDIF2,
    !                          UST, TAUHF, PHIHF, LLPHIHF)
    !          *KIJS*         - INDEX OF FIRST GRIDPOINT
    !          *KIJL*         - INDEX OF LAST GRIDPOINT
    !          *MIJ*          - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
    !          *LTAUWSHELTER* - if true then TAUWSHELTER
    !          *FL1*          - WAVE SPECTRUM.
    !          *AIRD*         - AIR DENSITY IN KG/M**3.
    !          *RNFAC*        - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
    !          *UFRIC*        - FRICTION VELOCITY
    !          *COSWDIF*      - COS(TH(K)-WDWAVE(IJ))
    !          *SINWDIF2*     - SIN(TH(K)-WDWAVE(IJ))**2
    !          *UST*          - REDUCED FRICTION VELOCITY DUE TO SHELTERING
    !          *Z0M*          - ROUGHNESS LENGTH
    !          *TAUHF*        - HIGH-FREQUENCY STRESS
    !          *PHIHF*        - HIGH-FREQUENCY ENERGY FLUX INTO OCEAN
    !          *LLPHIHF*      - TRUE IF PHIHF NEEDS TO COMPUTED
    
    
    !     METHOD.
    !     -------
    
    !       IT NEEDS A CALL TO INIT_X0TAUHF TO INITIALISE
    !       SEE REFERENCE FOR WAVE STRESS CALCULATION.
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    !       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWCOUP, ONLY: X0TAUHF_D, JTOT_TAUHF, WTAUHF_D, LLGCBZ0_D, LLNORMAGAM_D
    USE YOWFRED, ONLY: ZPIFR_D, FR5_D, TH_D, DELTH_D
    USE YOWPARAM, ONLY: NANG_D, NFRE_D
    USE YOWPCONS, ONLY: G_D, GM1_D, ZPI_D, ZPI4GM1_D, ZPI4GM2_D
    USE YOWPHYS, ONLY: ZALP_D, XKAPPA, TAUWSHELTER_D, GAMNCONST_D
    USE YOWTEST, ONLY: IU06
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    INTEGER(KIND=JWIM), INTENT(IN) :: MIJ(KIJL)
    LOGICAL, VALUE, INTENT(IN) :: LTAUWSHELTER
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: UFRIC, Z0M
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: AIRD
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: RNFAC
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(NANG_LOKI_PARAM) :: COSWDIF, SINWDIF2
    REAL(KIND=JWRB), INTENT(INOUT), DEVICE :: UST
    REAL(KIND=JWRB), INTENT(OUT), DEVICE :: TAUHF, PHIHF
    LOGICAL, VALUE, INTENT(IN) :: LLPHIHF
    
    
    INTEGER(KIND=JWIM) :: J, K
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: NS
    
    REAL(KIND=JWRB), PARAMETER :: ZSUPMAX = 0.0_JWRB    !  LOG(1.)
    REAL(KIND=JWRB) :: OMEGA, OMEGACC
    REAL(KIND=JWRB) :: X0G
    REAL(KIND=JWRB) :: YC, Y, CM1, ZX, ZARG, ZLOG, ZBETA
    REAL(KIND=JWRB) :: FNC, FNC2
    REAL(KIND=JWRB) :: GAMNORMA    ! RENORMALISATION FACTOR OF THE GROWTH RATE
    REAL(KIND=JWRB) :: ZNZ, CONFG
    REAL(KIND=JWRB) :: COSW, FCOSW2
    
    REAL(KIND=JWRB) :: XKS, OMS
    REAL(KIND=JWRB) :: SQRTZ0OG, ZSUP, ZINF, DELZ
    REAL(KIND=JWRB) :: TAUL, XLOGGZ0, SQRTGZ0
    REAL(KIND=JWRB) :: USTPH
    REAL(KIND=JWRB) :: CONST1, CONST2, CONSTTAU, CONSTPHI
    REAL(KIND=JWRB) :: F1DCOS2, F1DCOS3
    REAL(KIND=JWRB) :: F1D, F1DSIN2
    
    ! ----------------------------------------------------------------------
    
    
    IF (LLGCBZ0_D) THEN
!$loki inline
      CALL OMEGAGC_CUF_PARAMETRISE(UFRIC(IJ), NS, XKS, OMS)
    END IF
    
    !     See INIT_X0TAUHF
    X0G = X0TAUHF_D*G_D
    
    IF (LLPHIHF) USTPH = UST
    
    !*    COMPUTE THE INTEGRALS
    !     ---------------------
    
    XLOGGZ0 = LOG(G_D*Z0M(IJ))
    OMEGACC = MAX(ZPIFR_D(MIJ(IJ)), X0G / UST)
    SQRTZ0OG = SQRT(Z0M(IJ)*GM1_D)
    SQRTGZ0 = 1.0_JWRB / SQRTZ0OG
    YC = OMEGACC*SQRTZ0OG
    ZINF = LOG(YC)
    
    CONSTTAU = ZPI4GM2_D*FR5_D(MIJ(IJ))
    
    K = 1
    COSW = MAX(COSWDIF(K), 0.0_JWRB)
    FCOSW2 = FL1(IJ, K, MIJ(IJ))*COSW**2
    F1DCOS3 = FCOSW2*COSW
    F1DCOS2 = FCOSW2
    F1DSIN2 = FL1(IJ, K, MIJ(IJ))*SINWDIF2(K)
    F1D = FL1(IJ, K, MIJ(IJ))
    DO K=2,NANG_D
      COSW = MAX(COSWDIF(K), 0.0_JWRB)
      FCOSW2 = FL1(IJ, K, MIJ(IJ))*COSW**2
      F1DCOS3 = F1DCOS3 + FCOSW2*COSW
      F1DCOS2 = F1DCOS2 + FCOSW2
      F1DSIN2 = F1DSIN2 + FL1(IJ, K, MIJ(IJ))*SINWDIF2(K)
      F1D = F1D + FL1(IJ, K, MIJ(IJ))
    END DO
    F1DCOS3 = DELTH_D*F1DCOS3
    F1DCOS2 = DELTH_D*F1DCOS2
    F1DSIN2 = DELTH_D*F1DSIN2
    F1D = DELTH_D*F1D
    
    IF (LLNORMAGAM_D) THEN
      CONFG = GAMNCONST_D*FR5_D(MIJ(IJ))*RNFAC*SQRTGZ0
      CONST1 = CONFG*F1DSIN2
      CONST2 = CONFG*F1D
    ELSE
      CONST1 = 0.0_JWRB
      CONST2 = 0.0_JWRB
    END IF
    
    
    !     TAUHF :
    IF (LLGCBZ0_D) THEN
      ZSUP = MIN(LOG(OMS*SQRTZ0OG), ZSUPMAX)
    ELSE
      ZSUP = ZSUPMAX
    END IF
    
    TAUL = UST**2
    DELZ = MAX((ZSUP - ZINF) / REAL(JTOT_TAUHF - 1, kind=JWRB), 0.0_JWRB)
    TAUHF = 0.0_JWRB
    
    ! Intergrals are integrated following a change of variable : Z=LOG(Y)
    IF (LTAUWSHELTER) THEN
      DO J=1,JTOT_TAUHF
        Y = EXP(ZINF + REAL(J - 1, kind=JWRB)*DELZ)
        OMEGA = Y*SQRTGZ0
        CM1 = OMEGA*GM1_D
        ZX = UST*CM1 + ZALP_D
        ZARG = XKAPPA / ZX
        ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
        ZLOG = MIN(ZLOG, 0.0_JWRB)
        ZBETA = ZLOG**4*EXP(ZLOG)
        ZNZ = ZBETA*UST*Y
        GAMNORMA = (1.0_JWRB + CONST1*ZNZ) / (1.0_JWRB + CONST2*ZNZ)
        FNC2 = F1DCOS3*CONSTTAU*ZBETA*TAUL*WTAUHF_D(J)*DELZ*GAMNORMA
        TAUL = MAX(TAUL - TAUWSHELTER_D*FNC2, 0.0_JWRB)
        
        UST = SQRT(TAUL)
        TAUHF = TAUHF + FNC2
      END DO
    ELSE
      DO J=1,JTOT_TAUHF
        Y = EXP(ZINF + REAL(J - 1, kind=JWRB)*DELZ)
        OMEGA = Y*SQRTGZ0
        CM1 = OMEGA*GM1_D
        ZX = UST*CM1 + ZALP_D
        ZARG = XKAPPA / ZX
        ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
        ZLOG = MIN(ZLOG, 0.0_JWRB)
        ZBETA = ZLOG**4*EXP(ZLOG)
        FNC2 = ZBETA*WTAUHF_D(J)
        ZNZ = ZBETA*UST*Y
        GAMNORMA = (1.0_JWRB + CONST1*ZNZ) / (1.0_JWRB + CONST2*ZNZ)
        TAUHF = TAUHF + FNC2*GAMNORMA
      END DO
      TAUHF = F1DCOS3*CONSTTAU*TAUL*TAUHF*DELZ
    END IF
    
    
    PHIHF = 0.0_JWRB
    IF (LLPHIHF) THEN
      !       PHIHF:
      !       We are neglecting the gravity-capillary contribution
      !       Recompute DELZ over the full interval
      TAUL = USTPH**2
      ZSUP = ZSUPMAX
      DELZ = MAX((ZSUP - ZINF) / REAL(JTOT_TAUHF - 1, kind=JWRB), 0.0_JWRB)
      
      CONSTPHI = AIRD(IJ)*ZPI4GM1_D*FR5_D(MIJ(IJ))
      
      ! Intergrals are integrated following a change of variable : Z=LOG(Y)
      IF (LTAUWSHELTER) THEN
        DO J=1,JTOT_TAUHF
          Y = EXP(ZINF + REAL(J - 1, kind=JWRB)*DELZ)
          OMEGA = Y*SQRTGZ0
          CM1 = OMEGA*GM1_D
          ZX = USTPH*CM1 + ZALP_D
          ZARG = XKAPPA / ZX
          ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
          ZLOG = MIN(ZLOG, 0.0_JWRB)
          ZBETA = ZLOG**4*EXP(ZLOG)
          ZNZ = ZBETA*UST*Y
          GAMNORMA = (1.0_JWRB + CONST1*ZNZ) / (1.0_JWRB + CONST2*ZNZ)
          FNC2 = ZBETA*TAUL*WTAUHF_D(J)*DELZ*GAMNORMA
          TAUL = MAX(TAUL - TAUWSHELTER_D*F1DCOS3*CONSTTAU*FNC2, 0.0_JWRB)
          USTPH = SQRT(TAUL)
          PHIHF = PHIHF + FNC2 / Y
        END DO
        PHIHF = F1DCOS2*CONSTPHI*SQRTZ0OG*PHIHF
      ELSE
        DO J=1,JTOT_TAUHF
          Y = EXP(ZINF + REAL(J - 1, kind=JWRB)*DELZ)
          OMEGA = Y*SQRTGZ0
          CM1 = OMEGA*GM1_D
          ZX = USTPH*CM1 + ZALP_D
          ZARG = XKAPPA / ZX
          ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
          ZLOG = MIN(ZLOG, 0.0_JWRB)
          ZBETA = ZLOG**4*EXP(ZLOG)
          ZNZ = ZBETA*UST*Y
          GAMNORMA = (1.0_JWRB + CONST1*ZNZ) / (1.0_JWRB + CONST2*ZNZ)
          FNC2 = ZBETA*WTAUHF_D(J)*GAMNORMA
          PHIHF = PHIHF + FNC2 / Y
        END DO
        PHIHF = F1DCOS2*CONSTPHI*SQRTZ0OG*TAUL*PHIHF*DELZ
      END IF
    END IF
    
    
  END SUBROUTINE TAU_PHI_HF_CUF_PARAMETRISE
END MODULE TAU_PHI_HF_CUF_PARAMETRISE_MOD
MODULE MEANSQS_GC_CUF_PARAMETRISE_MOD
  !CONTAINED SUBROUTINES:
  ! - MEANSQS_GC
  CONTAINS
  SUBROUTINE MEANSQS_GC (XKMSS, KIJS, KIJL, HALP, USTAR, XMSSCG, FRGC)
    
    !***  DETERMINE MSS FOR GRAV-CAP WAVES UP TO WAVE NUMBER XKMSS
    
    !     AUTHOR: PETER JANSSEN
    !     ------
    
    !     REFERENCES:
    !     ----------
    
    !     VIERS PAPER EQ.(29)
    
    !----------------------------------------------------------------------
    
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWFRED, ONLY: NWAV_GC, XLOGKRATIOM1_GC, XKM_GC, VG_GC, C2OSQRTVG_GC, DELKCC_GC, DELKCC_GC_NS
    USE YOWPCONS, ONLY: G, ZPI, SURFT
    
    USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
    USE TAU_PHI_HF_MOD, ONLY: OMEGAGC
    
    !----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    
    REAL(KIND=JWRB), INTENT(IN) :: XKMSS    ! WAVE NUMBER CUT-OFF
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: HALP    ! 1/2 Phillips parameter
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: USTAR    ! friction velocity
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: XMSSCG    ! mean square slope for gravity-capillary waves
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: FRGC    ! Frequency from which the gravity-capillary spectrum is approximated
    
    
    INTEGER(KIND=JWIM) :: IJ, I, NE
    INTEGER(KIND=JWIM), DIMENSION(KIJL) :: NS
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    REAL(KIND=JWRB), DIMENSION(KIJL) :: XKS, OMS, COEF
    
    !     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
#include "gc_dispersion.h"
    
    ! ----------------------------------------------------------------------
    
    IF (LHOOK) CALL DR_HOOK('MEANSQS_GC', 0, ZHOOK_HANDLE)
    
    NE = MIN(MAX(NINT(LOG(XKMSS*XKM_GC(1))*XLOGKRATIOM1_GC), 1), NWAV_GC)
    
    DO IJ=KIJS,KIJL
      CALL OMEGAGC(USTAR(IJ), NS(IJ), XKS(IJ), OMS(IJ))
      FRGC(IJ) = OMS(IJ) / ZPI
      IF (XKS(IJ) > XKMSS) THEN
        NS(IJ) = NE
        XMSSCG(IJ) = 0.0_JWRB
      ELSE
        XMSSCG(IJ) = DELKCC_GC_NS(NS(IJ))*XKM_GC(NS(IJ))
      END IF
    END DO
    
    DO IJ=KIJS,KIJL
      DO I=NS(IJ) + 1,NE
        !         ANALYTICAL FORM INERTIAL SUB RANGE F(k) = k**(-4)*BB
        !         BB = COEF(IJ)*SQRT(VG_GC(I))/C_GC(I)**2
        !         mss :  integral of k**2 F(k)  k dk
        XMSSCG(IJ) = XMSSCG(IJ) + DELKCC_GC(I)*XKM_GC(I)
      END DO
      COEF(IJ) = C2OSQRTVG_GC(NS(IJ))*HALP(IJ)
      XMSSCG(IJ) = XMSSCG(IJ)*COEF(IJ)
    END DO
    
    IF (LHOOK) CALL DR_HOOK('MEANSQS_GC', 1, ZHOOK_HANDLE)
    
  END SUBROUTINE MEANSQS_GC
END MODULE MEANSQS_GC_CUF_PARAMETRISE_MOD
