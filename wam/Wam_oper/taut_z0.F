      SUBROUTINE TAUT_Z0(IJS, IJL, XLEV, UTOP, TAUW, USTAR, Z0)

! ----------------------------------------------------------------------

!**** *TAUT_Z0* - COMPUTATION OF TOTAL STRESS AND ROUGHNESS LENGTH SCALE.


!**   INTERFACE.
!     ----------

!       *CALL* *TAUT_Z0(IJS, IJL, XLEV, UTOP, TAUW, XLEV, USTAR, Z0)
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *XLEV* - WIND SPEED REFERENCE LEVEL (m)
!          *UTOP* - WIND SPEED AT REFERENCE LEVEL XLEV
!          *TAUW* - WAVE STRESS.
!          *USTAR*- FRICTION VELOCITY
!          *Z0*   - ROUGHNESS LENGTH 



!     METHOD.
!     -------

!       A STEADY STATE WIND PROFILE IS ASSUMED.
!       THE WIND STRESS IS COMPUTED USING THE ROUGHNESSLENGTH

!                  Z1=Z0/SQRT(1-TAUW/TAU)

!       WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-
!       INDUCED STRESS AND TAU IS THE TOTAL STRESS.
!       WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH TAUW/TAU < 1.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : ALPHA    ,XKAPPA  ,XNUM
      USE YOWPCONS , ONLY : G        ,EPSUS   ,ACD      ,BCD
      USE YOWTABL  , ONLY : JPLEVT   ,EPS1 
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL
      REAL, INTENT(IN) :: XLEV
      REAL,DIMENSION(IJS:IJL), INTENT(IN) :: UTOP, TAUW
      REAL,DIMENSION(IJS:IJL), INTENT(OUT) :: USTAR, Z0

      INTEGER, PARAMETER :: NITER=15

      REAL, PARAMETER :: TWOXMP1=3.

      INTEGER :: IJ, ITER

      REAL :: XLOGXL, ALPHAOG, XKUTOP, XOLOGZ0
      REAL :: USTOLD, TAUOLD, TAUNEW, X, F, DELF
      REAL :: USTM1, Z0TOT, Z0CH, Z0VIS, ZZ
      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('TAUT_Z0',0,ZHOOK_HANDLE)
#endif

      XLOGXL=LOG(XLEV)
      ALPHAOG=ALPHA/G

      DO IJ=IJS,IJL
        XKUTOP = XKAPPA*UTOP(IJ)
        USTOLD = UTOP(IJ)*SQRT(ACD+BCD*UTOP(IJ))
        TAUOLD = MAX(USTOLD**2,TAUW(IJ)+EPS1)
        USTAR(IJ) = SQRT(TAUOLD)
        USTM1 = 1./MAX(USTAR(IJ),EPSUS) 

        DO ITER=1,NITER
          X = TAUW(IJ)/TAUOLD
          Z0CH = ALPHAOG*TAUOLD/SQRT(1.-X)
          Z0VIS = XNUM*USTM1
          Z0TOT = Z0CH+Z0VIS

          XOLOGZ0= 1./(XLOGXL-LOG(Z0TOT))
          F = USTAR(IJ)-XKUTOP*XOLOGZ0
          ZZ = USTM1*(Z0CH*(2.-TWOXMP1*X)/(1.-X)-Z0VIS)/Z0TOT
          DELF= 1.-XKUTOP*XOLOGZ0**2*ZZ

          USTAR(IJ) = USTAR(IJ)-F/DELF
          TAUNEW = MAX(USTAR(IJ)**2,TAUW(IJ)+EPS1)
          USTAR(IJ) = SQRT(TAUNEW)
          IF(TAUNEW.EQ.TAUOLD) EXIT
          USTM1 = 1./MAX(USTAR(IJ),EPSUS) 
          TAUOLD = TAUNEW
        ENDDO

        Z0(IJ)=Z0CH

      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('TAUT_Z0',1,ZHOOK_HANDLE)
#endif
      END SUBROUTINE TAUT_Z0
