! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

FUNCTION STRESS_GC(ANG_GC, USTAR, Z0, Z0MIN, HALP, RNFAC) RESULT(STRESS_GC_RES)

!***  DETERMINE WAVE INDUCED STRESS FOR GRAV-CAP WAVES

!     AUTHOR: PETER JANSSEN
!     ------

!     REFERENCES:
!     ----------

!     VIERS PAPER EQ.(29)
!     FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

!----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LLNORMAGAM
      USE YOWFRED  , ONLY : NWAV_GC, OMEGA_GC, XK_GC, &
     &                      OMXKM3_GC, CM_GC, C2OSQRTVG_GC, XKMSQRTVGOC2_GC, &
     &                      OM3GMKM_GC, DELKCC_GC_NS, DELKCC_OMXKM3_GC
      USE YOWPCONS , ONLY : G, EPSUS, SURFT
      USE YOWPHYS  , ONLY : XKAPPA, ZALP, BETAMAXOXKAPPA2, BMAXOKAP, RN1_RN

      USE YOMHOOK  ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN) :: ANG_GC  ! factor to account for angular spreading of the input.
      REAL(KIND=JWRB), INTENT(IN) :: USTAR ! friction velocity
      REAL(KIND=JWRB), INTENT(IN) :: Z0 !  surface roughness
      REAL(KIND=JWRB), INTENT(IN) :: Z0MIN ! minimum surface roughness
      REAL(KIND=JWRB), INTENT(IN) :: HALP  ! 1/2 Phillips parameter
      REAL(KIND=JWRB), INTENT(IN) :: RNFAC  ! wind dependent factor used in the growth renormalisation


      INTEGER(KIND=JWIM) :: NS
      INTEGER(KIND=JWIM) :: I

      REAL(KIND=JWRB) :: XLAMBDA  ! Correction factor in the wave growth for gravity-capillary waves
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
      REAL(KIND=JWRB) :: GAMNORMA ! RENORMALISATION FACTOR OF THE GROWTH RATE
      REAL(KIND=JWRB) :: GAM_W
      REAL(KIND=JWRB) :: STRESS_GC_RES
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
   
!     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION RELATIONS
#include "gc_dispersion.h"
#include "ns_gc.intfb.h"

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('STRESS_GC',0,ZHOOK_HANDLE)

!*    1.0  DETERMINE GRAV_CAP SPECTRUM, TAUWHF.
!          ------------------------------------

!     FIND NS:
      NS = NS_GC(USTAR)

      TAUWCG_MIN = (USTAR*(Z0MIN/Z0))**2

      XLAMBDA = 1.0_JWRB + XLAMA * TANH(XLAMB * USTAR**NLAM)

      ZABHRC = ANG_GC * BETAMAXOXKAPPA2 * HALP * C2OSQRTVG_GC(NS)
      IF (LLNORMAGAM) THEN
        CONST = RNFAC * BMAXOKAP * HALP * C2OSQRTVG_GC(NS) /MAX(USTAR,EPSUS)
      ELSE
        CONST = 0.0_JWRB
      ENDIF

!     GAM_W GROWTHRATE BY WIND WITHOUT the multiplicative factor representing the ratio of air density to water density (eps)
!     and BETAMAXOXKAPPA2
      X       = USTAR*CM_GC(NS)
      XLOG    = LOG(XK_GC(NS)*Z0) + XKAPPA/(X + ZALP)
      ZLOG    = XLOG - LOG(XLAMBDA) 
      ZLOG    = MIN(ZLOG, 0.0_JWRB)
      ZLOG2X  = ZLOG*ZLOG*X
      GAM_W   = ZLOG2X*ZLOG2X*EXP(XLOG)*OM3GMKM_GC(NS)

      ZN = CONST*XKMSQRTVGOC2_GC(NS)*GAM_W
      GAMNORMA = (1.0_JWRB + RN1_RN*ZN)/(1.0_JWRB + ZN)

      TAUWCG = GAM_W * DELKCC_GC_NS(NS) * OMXKM3_GC(NS) * GAMNORMA

      DO I = NS+1, NWAV_GC
!       ANALYTICAL FORM INERTIAL SUB RANGE F(k) = k**(-4)*BB
!       BB = HALP * C2OSQRTVG_GC(NS)*SQRT(VG_GC(I))/C_GC(I)**2
!       Tauwcg : (rhow * g /rhoa) * integral of (1/c) * gammma * F(k)  k dk 
!       with omega=g*k and omega=k*c,  then
!       Tauwcg : (rhow /rhoa) * integral of omega * gammma * F(k)  k dk
!       but gamma is computed wihtout the rhoa/rhow factor so
!       Tauwcg : integral of omega * gammma_wam * F(k)  k dk
!       It should be done in vector form with actual directional spreading information
!       It simplified here by using the ANG_GC factor.

        X       = USTAR*CM_GC(I)
        XLOG    = LOG(XK_GC(I)*Z0) + XKAPPA/(X + ZALP)
        ZLOG    = XLOG - LOG(XLAMBDA) 
        ZLOG    = MIN(ZLOG, 0.0_JWRB)
        ZLOG2X  = ZLOG*ZLOG*X
        GAM_W   = ZLOG2X*ZLOG2X*EXP(XLOG)*OM3GMKM_GC(I)

        ZN = CONST*XKMSQRTVGOC2_GC(I)*GAM_W
        GAMNORMA = (1.0_JWRB + RN1_RN*ZN)/(1.0_JWRB + ZN)

        TAUWCG = TAUWCG + GAM_W * DELKCC_OMXKM3_GC(I) * GAMNORMA
      ENDDO
      STRESS_GC_RES = MAX(ZABHRC * TAUWCG, TAUWCG_MIN)

IF (LHOOK) CALL DR_HOOK('STRESS_GC',1,ZHOOK_HANDLE)
 
END FUNCTION STRESS_GC
