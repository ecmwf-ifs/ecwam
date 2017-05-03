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

      USE YOWCOUP  , ONLY : ALPHA    ,XKAPPA
      USE YOWPCONS , ONLY : G        ,ACD      ,BCD
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

        DO ITER=1,NITER
          X = TAUW(IJ)/TAUOLD
          Z0(IJ) = ALPHAOG*TAUOLD/SQRT(1.-X)
          XOLOGZ0= 1./(XLOGXL-LOG(Z0(IJ)))
          F = USTAR(IJ)-XKUTOP*XOLOGZ0
          DELF = 1.-XKUTOP*XOLOGZ0**2*(2.-TWOXMP1*X)/(USTAR(IJ)*(1.-X))
          USTAR(IJ) = USTAR(IJ)-F/DELF
          TAUNEW = MAX(USTAR(IJ)**2,TAUW(IJ)+EPS1)
          USTAR(IJ) = SQRT(TAUNEW)
          IF(TAUNEW.EQ.TAUOLD) EXIT
          TAUOLD = TAUNEW
        ENDDO

      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('TAUT_Z0',1,ZHOOK_HANDLE)
#endif
      END SUBROUTINE TAUT_Z0
