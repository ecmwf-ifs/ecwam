      SUBROUTINE WSIGSTAR (IJS, IJL, USNEW, Z0NEW, WSTAR, SIG_N)
! ----------------------------------------------------------------------

!**** *WSIGSTAR* - COMPUTATION OF THE RELATIVE STANDARD DEVIATION OF USTAR.

!*    PURPOSE.
!     ---------

!     COMPUTES THE STANDARD DEVIATION OF USTAR DUE TO SMALL SCALE GUSTINESS
!     RELATIVE TO USTAR

!**   INTERFACE.
!     ----------

!     *CALL* *WSIGSTAR (IJS, IJL, USNEW, Z0NEW, WSTAR, SIG_N)
!             *IJS*   - INDEX OF FIRST GRIDPOINT.
!             *IJL*   - INDEX OF LAST GRIDPOINT.
!             *USNEW* - NEW FRICTION VELOCITY IN M/S.
!             *Z0NEW* - ROUGHNESS LENGTH IN M.
!             *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
!             *SIG_N* - ESTINATED RELATIVE STANDARD DEVIATION OF USTAR.


!     METHOD.
!     -------

!     USE PANOFSKY (1991) TO EXPRESS THE STANDARD DEVIATION OF U10 IN TERMS
!     USTAR AND (Zi/L) THE MONIN-OBOKHOV LENGTH (Zi THE INVERSION HEIGHT).
!     (but with the background gustiness set to 0.)
!     and USTAR=SQRT(Cd)*U10 to DERIVE THE STANDARD DEVIATION OF USTAR.
!     WITH CD=A+B*U10 (see below).

!     EXTERNALS.
!     ----------

!       NONE.

!     MODIFICATIONS
!     -------------

!     REFERENCE.
!     ----------

!     SEE SECTION 3.2.1 OF THE WAM DOCUMENTATION.
!     USE HERSBACH 2011 FOR CD(U10) (SEE ALSO EDSON et al. 2013)

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : XKAPPA
      USE YOWWIND  , ONLY : WSPMIN
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      REAL, PARAMETER :: BG_GUST = 0.0 ! NO BACKGROUND GUSTINESS (S0 12. IS NOT USED)
      REAL, PARAMETER :: ONETHIRD = 1./3.
      REAL, PARAMETER :: LOG10 = LOG(10.)
      REAL, PARAMETER :: C1 = 1.03E-3
      REAL, PARAMETER :: C2 = 0.04E-3
      REAL, PARAMETER :: P1 = 1.48
      REAL, PARAMETER :: P2 = -0.21

      INTEGER :: IJ,IJS,IJL

      REAL :: U10, C_D, DC_DDU, SIG_CONV
      REAL :: XKAPPAD
      REAL :: U10M1, C2U10P1, U10P2
      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(IJS:IJL) :: USNEW, Z0NEW, WSTAR
      REAL, DIMENSION(IJS:IJL) :: SIG_N

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WSIGSTAR',0,ZHOOK_HANDLE)
#endif

! ----------------------------------------------------------------------


      XKAPPAD=1.0/XKAPPA
      DO IJ=IJS,IJL
!
!       IN THE FOLLOWING U10 IS ESTIMATED ASSUMING EVERYTHING IS
!       BASED ON U*
!
        U10 = USNEW(IJ)*XKAPPAD*(LOG10-LOG(Z0NEW(IJ)))
        U10 = MAX(U10,WSPMIN)
        U10M1=1./U10
        C2U10P1=C2*U10**P1
        U10P2=U10**P2
        C_D = (C1 + C2U10P1)*U10P2
        DC_DDU = (P2*C1+(P1+P2)*C2U10P1)*U10P2*U10M1
        SIG_CONV = 1. + 0.5*U10/C_D*DC_DDU
        SIG_N(IJ) = MIN(0.1, SIG_CONV * U10M1 * 
     &                        (BG_GUST*USNEW(IJ)**3+
     &                         0.5*XKAPPA*WSTAR(IJ)**3)**ONETHIRD
     &                 )
      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WSIGSTAR',1,ZHOOK_HANDLE)
#endif

      RETURN
      END SUBROUTINE WSIGSTAR
