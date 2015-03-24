      SUBROUTINE H_MAX(C3,C4,NW,IJS,IJL,HMAX_N)

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

!     C_3            REAL         SKEWNESS
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

!     FURTHER UPDATES TO THE FREAK WAVE WARNING SYSTEM. PETER A.E.M. 
!     JANSSEN AND J.-R. BIDLOT, TO BE PUBLISHED AS ECMWF TECH MEMO
!
!----------------------------------------------------------------------

      USE YOWPCONS , ONLY : PI     , C4MAX    ,C4MIN 
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

!----------------------------------------------------------------------
 
      IMPLICIT NONE
      INTEGER :: IJ, IJS, IJL
      INTEGER, DIMENSION(IJS:IJL):: NW
      REAL :: GAMMA_E, ZETA_3, C3M, C4M, XN, Z0, C_Z, D_Z, E_Z, ARG, Z 
      REAL :: ALPHA, BETA, G1, G2, G3
      REAL :: ZHOOK_HANDLE

      REAL, DIMENSION(IJS:IJL):: C3,C4,HMAX_N

!----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('H_MAX',0,ZHOOK_HANDLE)
#endif

!*    1. DETERMINE EXPECTED MAXIMUM WAVE HEIGHT.
!     -----------------------------------------
    
      GAMMA_E = 0.57721566  ! EULER CONSTANT
      ZETA_3  = 1.20206     ! RIEMANN ZETA FUNCTION
!
!     G1, G2, G3 ARE FIRST, SECOND AND THIRD DERIVATIVE OF
!     GAMMA FUNCTION AT Z=1.
!
      G1 = -GAMMA_E
      G2 = (GAMMA_E**2+PI**2/6.)
      G3 = -2.*ZETA_3-GAMMA_E*PI**2/2.-GAMMA_E**3

      DO IJ=IJS,IJL

        IF (NW(IJ).GT.1) THEN
          C3M = C3(IJ)
          C4M = MIN(C4MAX,C4(IJ))
          C4M = MAX(C4MIN,C4M)
          XN = REAL(NW(IJ))

          Z0 = 0.5*LOG(XN)
          ALPHA = 2.*Z0*(Z0-1.)+(1.-2.*Z0)*G1+0.5*G2

          C_Z = Z0*(4.*Z0**2-12.*Z0+6.)-(6.*Z0**2-12.*Z0+3.)*G1
          D_Z = 3.*(Z0-1.)*G2 
          E_Z = -0.5*G3
          BETA = C_Z+D_Z+E_Z

          ARG = 1.+C4M*ALPHA+C3M**2*BETA       
          ARG = MAX(ARG,0.1)

          HMAX_N(IJ) = SQRT(Z0)*(1.+(GAMMA_E+LOG(ARG))/(4.*Z0))
        ELSE
          HMAX_N(IJ) = 0.
        ENDIF

      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('H_MAX',1,ZHOOK_HANDLE)
#endif
      RETURN

      END SUBROUTINE H_MAX
