! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE TRANSF_SNL_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) FUNCTION TRANSF_SNL_CUF_PARAMETRISE (XK0, D, XNU, SIG_TH)
    
    !***  DETERMINE NARROW BAND LIMIT NONLINEAR TRANSFER FUNCTION
    
    !     VARIABLE       TYPE         PURPOSE
    !     --------       ----         -------
    
    !     XK0            REAL         WAVE NUMBER
    !     D              REAL         DEPTH
    !     XNU            REAL         RELATIVE SPECTRAL WIDTH
    !     SIG_TH         REAL         RELATIVE WIDTH in DIRECTION
    
    !----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWPCONS, ONLY: G_D, DKMAX
    USE YOWSHAL, ONLY: BATHYMAX, XKDMIN
    
    
    !----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    REAL(KIND=JWRB) :: TRANSF_SNL_CUF_PARAMETRISE
    REAL(KIND=JWRB), INTENT(IN) :: XK0, D, XNU, SIG_TH
    
    REAL(KIND=JWRB), PARAMETER :: EPS = 0.0001_JWRB
    REAL(KIND=JWRB), PARAMETER :: TRANSF_SNL_MIN = 0.1_JWRB
    REAL(KIND=JWRB), PARAMETER :: TRANSF_SNL_MAX = 10._JWRB
    
!$loki routine seq
    REAL(KIND=JWRB) :: X, XK, T_0, T_0_SQ, OM, C_0, V_G, V_G_SQ, DV_G
    REAL(KIND=JWRB) :: XNL_1, XNL_2, XNL_3, XNL_4, XNL
    REAL(KIND=JWRB) :: C_S_SQ, ALP, ZFAC
    
    !----------------------------------------------------------------------
    
    
    !*    1. DETERMINE TRANSFER FUNCTION.
    !     ------------------------------
    
    IF (D < BATHYMAX .and. D > 0._JWRB) THEN
      X = XK0*D
      IF (X > DKMAX) THEN
        TRANSF_SNL_CUF_PARAMETRISE = 1._JWRB
      ELSE
        XK = MAX(XK0, XKDMIN / D)
        X = XK*D
        T_0 = TANH(X)
        T_0_SQ = T_0**2
        OM = SQRT(G_D*XK*T_0)
        C_0 = OM / XK
        C_S_SQ = G_D*D
        IF (X < EPS) THEN
          V_G = C_0
        ELSE
          V_G = 0.5_JWRB*C_0*(1._JWRB + 2._JWRB*X / SINH(2._JWRB*X))
        END IF
        V_G_SQ = V_G**2
        DV_G = (T_0 - X*(1. - T_0_SQ))**2 + 4._JWRB*X**2*T_0_SQ*(1._JWRB - T_0_SQ)
        
        XNL_1 = (9._JWRB*T_0_SQ**2 - 10._JWRB*T_0_SQ + 9._JWRB) / (8._JWRB*T_0_SQ*T_0)
        XNL_2 = ((2._JWRB*V_G - 0.5_JWRB*C_0)**2 / (G_D*D - V_G_SQ) + 1._JWRB) / X
        XNL_4 = 1. / (4._JWRB*T_0)*(2._JWRB*C_0 + V_G*(1._JWRB - T_0_SQ))**2 / (C_S_SQ - V_G_SQ)
        ALP = (1. - V_G_SQ / C_S_SQ)*C_0**2 / V_G_SQ
        ZFAC = SIG_TH**2 / (SIG_TH**2 + ALP*XNU**2)
        XNL_3 = ZFAC*XNL_4
        
        XNL = XNL_1 - XNL_2 + XNL_3
        TRANSF_SNL_CUF_PARAMETRISE = XNL**2 / (DV_G*T_0_SQ**4)
        TRANSF_SNL_CUF_PARAMETRISE = MAX(MIN(TRANSF_SNL_MAX, TRANSF_SNL_CUF_PARAMETRISE), TRANSF_SNL_MIN)
      END IF
    ELSE
      TRANSF_SNL_CUF_PARAMETRISE = 1._JWRB
    END IF
    
    
    
  END FUNCTION TRANSF_SNL_CUF_PARAMETRISE
END MODULE TRANSF_SNL_CUF_PARAMETRISE_MOD
