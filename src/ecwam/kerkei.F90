! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE KERKEI (X, KER, KEI)
!**********************************************************************
! Computes the values of the zeroth order Kelvin function Ker and Kei
! These functions are used to determine the friction factor fw as a
! function of the bottom roughness length assuming a linear profile
! of eddy viscosity (See Grant and Madsen, 1979)
!**********************************************************************

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE
#include "kzeone.intfb.h"

      REAL(KIND=JWRU) :: ZR, ZI, CYR, CYI, CYR1, CYI1

      INTEGER(KIND=JWIM) :: NZ, IERR
      REAL(KIND=JWRB) :: X, KER, KEI

      ZR = X*0.50_JWRU*SQRT(2.0_JWRU)
      ZI = ZR
      CALL KZEONE (ZR, ZI, CYR, CYI, CYR1, CYI1)
      KER = CYR/EXP(ZR)
      KEI = CYI/EXP(ZR)
      END SUBROUTINE KERKEI
