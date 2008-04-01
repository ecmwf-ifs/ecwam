      SUBROUTINE H_MAX(C4,NW,IJS,IJL,HMAX_N)

!***  DETERMINE EXPECTED MAXIMUM WAVE HEIGHT.

!     PURPOSE:
!     -------

!     DETERMINE EXPECTED MAXIMUM WAVE HEIGHT ACCORDING TO
!                 /
!         <H_MAX>=| dh h p(h)
!                 /
!     WHERE p(h) IS DERIVED IN MORI AND JANSSEN(2006)

!     VARIABLE       TYPE         PURPOSE
!     --------       ----         -------

!     C_4            REAL         KURTOSIS(<ETA>^4/(3<ETA^2>^2)-1)
!     N_W            INTEGER      NUMBER OF WAVES
!     IJS            INTEGER      FIRST INDEX
!     IJL            INTEGER      LAST INDEX

!     OUTPUT:
!     ------

!     HMAX_N         REAL         MAXIMUM WAVE HEIGHT NORMALIZED WITH
!                                 SIGNIFICANT WAVE HEIGHT H_S. 

!     AUTHOR:
!     ------

!     P.A.E.M. JANSSEN, JANUARY 2008

!----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER :: IJ, IJS, IJL
      REAL :: GAMMA_E, PI, CONST, C4M, XN, Z0, B_Z, ARG, Z 
      INTEGER, DIMENSION(IJS:IJL):: NW
      REAL, DIMENSION(IJS:IJL):: C4,HMAX_N

!*    1. DETERMINE EXPECTED MAXIMUM WAVE HEIGHT.
!     -----------------------------------------
    
      GAMMA_E = 0.57721566  ! EULER CONSTANT
      PI = 4.*ATAN(1.)

      CONST = 0.5*(GAMMA_E**2+PI**2/6.)

      DO IJ=IJS,IJL
         C4M = MIN(1.,C4(IJ))
         C4M = MAX(-0.33,C4M)
         XN = FLOAT(NW(IJ))

         Z0 = 0.5*LOG(XN)
         B_Z = 2.*Z0*(Z0-1.)
         ARG = 1.+C4M*(B_Z-GAMMA_E*(1.-2.*Z0)-CONST)
         ARG = MAX(ARG,0.1)
         Z = Z0+0.5*(GAMMA_E+LOG(ARG))

         HMAX_N(IJ) = SQRT(Z)
      ENDDO

      RETURN

      END SUBROUTINE H_MAX
