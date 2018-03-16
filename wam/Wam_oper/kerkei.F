! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE KERKEI(X,KER,KEI)
!**********************************************************************
! Computes the values of the zeroth order Kelvin function Ker and Kei
! These functions are used to determine the friction factor fw as a 
! function of the bottom roughness length assuming a linear profile
! of eddy viscosity (See Grant and Madsen, 1979)
!**********************************************************************
      IMPLICIT NONE
#include "kzeone.intfb.h"
   
      DOUBLE PRECISION ZR,ZI,CYR,CYI,CYR1,CYI1
      INTEGER NZ,IERR
      REAL X,KER,KEI
   
      ZR=X*.50D0*SQRT(2.0D0)
      ZI=ZR
      CALL KZEONE(ZR, ZI, CYR, CYI,CYR1,CYI1)
      KER=CYR/EXP(ZR)
      KEI=CYI/EXP(ZR)
      END SUBROUTINE KERKEI
