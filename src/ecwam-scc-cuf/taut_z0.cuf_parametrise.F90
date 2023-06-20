! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE TAUT_Z0_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE TAUT_Z0_CUF_PARAMETRISE (KIJS, KIJL, IUSFG, HALP, UTOP, UDIR, TAUW, TAUWDIR, RNFAC, USTAR, Z0,  &
  & Z0B, CHRNCK, IJ)
    
    ! ----------------------------------------------------------------------
    
    !**** *TAUT_Z0* - COMPUTATION OF TOTAL STRESS AND ROUGHNESS LENGTH SCALE.
    
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *TAUT_Z0(KIJS, KIJL, IUSFG, FL1, WAVNUM,
    !                       UTOP, UDIR, TAUW, TAUWDIR, RNFAC,
    !                       USTAR, Z0, Z0B, CHRNCK)
    !          *KIJS*    - INDEX OF FIRST GRIDPOINT
    !          *KIJL*    - INDEX OF LAST GRIDPOINT
    !          *IUSFG*   - IF = 1 THEN USE THE FRICTION VELOCITY (US) AS FIRST GUESS in TAUT_Z0
    !                           0 DO NOT USE THE FIELD USTAR
    !          *FL1*     - 2D-SPECTRA
    !          *WAVNUM*  - WAVE NUMBER
    !          *HALP*    - 1/2 PHILLIPS PARAMETER
    !          *UTOP*    - WIND SPEED AT REFERENCE LEVEL XNLEV
    !          *UDIR*    - WIND SPEED DIRECTION AT REFERENCE LEVEL XNLEV
    !          *TAUW*    - WAVE STRESS.
    !          *TAUWDIR* - WAVE STRESS DIRECTION.
    !          *RNFAC*   - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
    !          *USTAR*   - FRICTION VELOCITY
    !          *Z0*      - ROUGHNESS LENGTH
    !          *Z0B*     - BACKGROUND ROUGHNESS LENGTH
    !          *CHRNCK*  - CHARNOCK COEFFICIENT
    
    !     METHOD.
    !     -------
    
    !       A STEADY STATE WIND PROFILE IS ASSUMED.
    !       THE WIND STRESS IS COMPUTED USING THE ROUGHNESS LENGTH
    
    !                  Z1=Z0/SQRT(1-TAUW/TAU)
    
    !       WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-
    !       INDUCED STRESS AND TAU IS THE TOTAL STRESS.
    !       WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH TAUW/TAU < 1.
    
    !       IT WAS EXTENDED TO INCLUDE THE GRAVITY-CAPILLARY MODEL FOR THE CALCULATION
    !       OF THE BACKGROUND ROUGHNESS.
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    !       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
    
    ! ----------------------------------------------------------------------
    
    USE STRESS_GC_CUF_PARAMETRISE_MOD, ONLY: STRESS_GC_CUF_PARAMETRISE
    USE CHNKMIN_CUF_PARAMETRISE_MOD, ONLY: CHNKMIN_CUF_PARAMETRISE
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWCOUP, ONLY: LLCAPCHNK_D, LLGCBZ0_D
    USE YOWPARAM, ONLY: NANG_D, NFRE_D
    USE YOWPCONS, ONLY: G_D, GM1_D, EPSUS, EPSMIN, ACD, BCD, CDMAX
    USE YOWPHYS, ONLY: XKAPPA, XNLEV, RNU_D, RNUM_D, ALPHA_D, ALPHAMIN_D, ALPHAMAX, ANG_GC_A_D, ANG_GC_B_D, ANG_GC_C_D
    USE YOWTABL, ONLY: EPS1
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL, IUSFG
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: UTOP, UDIR, TAUW, TAUWDIR
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: HALP, RNFAC
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: USTAR
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: Z0, Z0B, CHRNCK
    
    
    INTEGER(KIND=JWIM), PARAMETER :: NITER = 17
    
    REAL(KIND=JWRB), PARAMETER :: TWOXMP1 = 3.0_JWRB
    
    INTEGER(KIND=JWIM) :: ITER
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: IFRPH
    
    ! Cd and Z0 from Hersbach 2010, ECMWF Tech Memo (without the viscous part)
    !     CD = ACDLIN + BCDLIN*SQRT(PCHAR) * U10
    REAL(KIND=JWRB), PARAMETER :: ACDLIN = 0.0008_JWRB
    REAL(KIND=JWRB), PARAMETER :: BCDLIN = 0.00047_JWRB
    REAL(KIND=JWRB) :: ALPHAGM1
    
    REAL(KIND=JWRB), PARAMETER :: Z0MIN = 0.000001_JWRB
    REAL(KIND=JWRB) :: PCE_GC
    REAL(KIND=JWRB) :: Z0MINRST
    REAL(KIND=JWRB) :: CHARNOCK_MIN
    REAL(KIND=JWRB) :: COSDIFF
    REAL(KIND=JWRB) :: ZCHAR
    REAL(KIND=JWRB) :: US2TOTAUW, USMAX
    REAL(KIND=JWRB) :: XLOGXL, XKUTOP, XOLOGZ0
    REAL(KIND=JWRB) :: USTOLD, USTNEW, TAUOLD, TAUNEW, X, F, DELF, CDFG
    REAL(KIND=JWRB) :: USTM1, Z0TOT, Z0CH, Z0VIS, HZ0VISO1MX, ZZ
    REAL(KIND=JWRB) :: CONST, TAUV, DEL
    REAL(KIND=JWRB) :: RNUEFF, RNUKAPPAM1
    REAL(KIND=JWRB) :: ALPHAOG, XMIN
    REAL(KIND=JWRB) :: W1
    REAL(KIND=JWRB) :: TAUWACT, TAUWEFF
    REAL(KIND=JWRB) :: ANG_GC, TAUUNR
    
    LOGICAL :: LLCOSDIFF
    
    ! ----------------------------------------------------------------------
    
    
    XLOGXL = LOG(XNLEV)
    US2TOTAUW = 1.0_JWRB + EPS1
    
    !     ONLY take the contribution of TAUW that is in the wind direction
    COSDIFF = COS(UDIR(IJ) - TAUWDIR(IJ))
    TAUWACT = MAX(TAUW(IJ)*COSDIFF, EPSMIN)
    LLCOSDIFF = COSDIFF > 0.9_JWRB
    
    !  USING THE CG MODEL:
    IF (LLGCBZ0_D) THEN
      
      IF (LLCAPCHNK_D) THEN
        CHARNOCK_MIN = CHNKMIN_CUF_PARAMETRISE(UTOP(IJ))
        ALPHAOG = CHARNOCK_MIN*GM1_D
      ELSE
        ALPHAOG = 0.0_JWRB
      END IF
      
      USMAX = MAX(-0.21339_JWRB + 0.093698_JWRB*UTOP(IJ) - 0.0020944_JWRB*UTOP(IJ)**2 + 5.5091E-5_JWRB*UTOP(IJ)**3, 0.03_JWRB)
      TAUWEFF = MIN(TAUWACT*US2TOTAUW, USMAX**2)
      
      RNUEFF = 0.04_JWRB*RNU_D
      
      RNUKAPPAM1 = RNUEFF / XKAPPA
      
      PCE_GC = 0.001_JWRB*IUSFG + (1 - IUSFG)*0.005_JWRB
      
      IF (IUSFG == 0) THEN
        ALPHAGM1 = ALPHA_D*GM1_D
        IF (UTOP(IJ) < 1.0_JWRB) THEN
          CDFG = 0.002_JWRB
        ELSE IF (LLCOSDIFF) THEN
          X = MIN(TAUWACT / MAX(USTAR(IJ), EPSUS)**2, 0.99_JWRB)
          ZCHAR = MIN(ALPHAGM1*USTAR(IJ)**2 / SQRT(1.0_JWRB - X), 0.05_JWRB*EXP(-0.05_JWRB*(UTOP(IJ) - 35._JWRB)))
          ZCHAR = MIN(ZCHAR, ALPHAMAX)
          CDFG = ACDLIN + BCDLIN*SQRT(ZCHAR)*UTOP(IJ)
        ELSE
          CDFG = MAX(MIN(0.0006_JWRB + 0.00008_JWRB*UTOP(IJ), 0.001_JWRB + 0.0018_JWRB*EXP(-0.05_JWRB*(UTOP(IJ) - 33._JWRB))),  &
          & 0.001_JWRB)
        END IF
        USTAR(IJ) = UTOP(IJ)*SQRT(CDFG)
      END IF
      
      W1 = 0.85_JWRB - 0.05_JWRB*(TANH(10.0_JWRB*(UTOP(IJ) - 5.0_JWRB)) + 1.0_JWRB)
      
      XKUTOP = XKAPPA*UTOP(IJ)
      
      USTOLD = USTAR(IJ)
      TAUOLD = USTOLD**2
      
      DO ITER=1,NITER
        !         Z0 IS DERIVED FROM THE NEUTRAL LOG PROFILE: UTOP = (USTAR/XKAPPA)*LOG((XNLEV+Z0)/Z0)
        Z0(IJ) = MAX(XNLEV / (EXP(MIN(XKUTOP / USTOLD, 50.0_JWRB)) - 1.0_JWRB), Z0MIN)
        ! Viscous kinematic stress nu_air * dU/dz at z=0 of the neutral log profile reduced by factor 25 (0.04)
        TAUV = RNUKAPPAM1*USTOLD / Z0(IJ)
        
        ANG_GC = ANG_GC_A_D + ANG_GC_B_D*TANH(ANG_GC_C_D*TAUOLD)
        
        TAUUNR = STRESS_GC_CUF_PARAMETRISE(ANG_GC, USTAR(IJ), Z0(IJ), Z0MIN, HALP, RNFAC)
        
        !         TOTAL kinematic STRESS:
        TAUNEW = TAUWEFF + TAUV + TAUUNR
        USTNEW = SQRT(TAUNEW)
        USTAR(IJ) = W1*USTOLD + (1.0_JWRB - W1)*USTNEW
        
        !         CONVERGENCE ?
        DEL = USTAR(IJ) - USTOLD
        IF (ABS(DEL) < PCE_GC*USTAR(IJ)) EXIT
        TAUOLD = USTAR(IJ)**2
        USTOLD = USTAR(IJ)
      END DO
      ! protection just in case there is no convergence
      IF (ITER > NITER) THEN
        CDFG = MAX(MIN(0.0006_JWRB + 0.00008_JWRB*UTOP(IJ), 0.001_JWRB + 0.0018_JWRB*EXP(-0.05_JWRB*(UTOP(IJ) - 33._JWRB))),  &
        & 0.001_JWRB)
        USTAR(IJ) = UTOP(IJ)*SQRT(CDFG)
        Z0MINRST = USTAR(IJ)**2*ALPHA_D*GM1_D
        Z0(IJ) = MAX(XNLEV / (EXP(XKUTOP / USTAR(IJ)) - 1.0_JWRB), Z0MINRST)
        Z0B(IJ) = Z0MINRST
      ELSE
        Z0(IJ) = MAX(XNLEV / (EXP(XKUTOP / USTAR(IJ)) - 1.0_JWRB), Z0MIN)
        Z0B(IJ) = Z0(IJ)*SQRT(TAUUNR / TAUOLD)
      END IF
      
      !       Refine solution
      X = TAUWEFF / TAUOLD
      
      IF (X < 0.99_JWRB) THEN
        USTOLD = USTAR(IJ)
        TAUOLD = MAX(USTOLD**2, TAUWEFF)
        
        DO ITER=1,NITER
          X = MIN(TAUWEFF / TAUOLD, 0.99_JWRB)
          USTM1 = 1.0_JWRB / MAX(USTOLD, EPSUS)
          !!!! Limit how small z0 could become
          !!!! This is a bit of a compromise to limit very low Charnock for intermediate high winds (15 -25 m/s)
          !!!! It is not ideal !!!
          Z0(IJ) = MAX(XNLEV / (EXP(MIN(XKUTOP / USTOLD, 50.0_JWRB)) - 1.0_JWRB), Z0MIN)
          
          TAUUNR = STRESS_GC_CUF_PARAMETRISE(ANG_GC, USTOLD, Z0(IJ), Z0MIN, HALP, RNFAC)
          
          Z0B(IJ) = MAX(Z0(IJ)*SQRT(TAUUNR / TAUOLD), ALPHAOG*TAUOLD)
          Z0VIS = RNUM_D*USTM1
          HZ0VISO1MX = 0.5_JWRB*Z0VIS / (1.0_JWRB - X)
          Z0(IJ) = HZ0VISO1MX + SQRT(HZ0VISO1MX**2 + Z0B(IJ)**2 / (1.0_JWRB - X))
          
          XOLOGZ0 = 1.0_JWRB / (XLOGXL - LOG(Z0(IJ)))
          F = USTOLD - XKUTOP*XOLOGZ0
          ZZ = 2.0_JWRB*USTM1*(3.0_JWRB*Z0B(IJ)**2 + 0.5_JWRB*Z0VIS*Z0(IJ) - Z0(IJ)**2) / (2.0_JWRB*Z0(IJ)**2*(1.0_JWRB - X) -  &
          & Z0VIS*Z0(IJ))
          
          DELF = 1.0_JWRB - XKUTOP*XOLOGZ0**2*ZZ
          IF (DELF /= 0.0_JWRB) USTAR(IJ) = USTOLD - F / DELF
          
          !           CONVERGENCE ?
          DEL = USTAR(IJ) - USTOLD
          
          IF (ABS(DEL) < PCE_GC*USTAR(IJ)) EXIT
          USTOLD = USTAR(IJ)
          TAUOLD = MAX(USTOLD**2, TAUWEFF)
        END DO
        ! protection just in case there is no convergence
        IF (ITER > NITER) THEN
          CDFG = MAX(MIN(0.0006_JWRB + 0.00008_JWRB*UTOP(IJ), 0.001_JWRB + 0.0018_JWRB*EXP(-0.05_JWRB*(UTOP(IJ) - 33._JWRB))),  &
          & 0.001_JWRB)
          USTAR(IJ) = UTOP(IJ)*SQRT(CDFG)
          Z0MINRST = USTAR(IJ)**2*ALPHA_D*GM1_D
          Z0(IJ) = MAX(XNLEV / (EXP(XKUTOP / USTAR(IJ)) - 1.0_JWRB), Z0MINRST)
          Z0B(IJ) = Z0MINRST
          CHRNCK(IJ) = MAX(G_D*Z0(IJ) / USTAR(IJ)**2, ALPHAMIN_D)
        ELSE
          CHRNCK(IJ) = MAX(G_D*(Z0B(IJ) / SQRT(1.0_JWRB - X)) / MAX(USTAR(IJ), EPSUS)**2, ALPHAMIN_D)
        END IF
        
      ELSE
        USTM1 = 1.0_JWRB / MAX(USTAR(IJ), EPSUS)
        Z0VIS = RNUM_D*USTM1
        CHRNCK(IJ) = MAX(G_D*(Z0(IJ) - Z0VIS)*USTM1**2, ALPHAMIN_D)
      END IF
      
      
      
    ELSE
      
      TAUWEFF = TAUWACT*US2TOTAUW
      
      IF (LLCAPCHNK_D) THEN
        CHARNOCK_MIN = CHNKMIN_CUF_PARAMETRISE(UTOP(IJ))
        XMIN = 0.15_JWRB*(ALPHA_D - CHARNOCK_MIN)
        ALPHAOG = CHARNOCK_MIN*GM1_D
      ELSE
        XMIN = 0.0_JWRB
        ALPHAOG = ALPHA_D*GM1_D
      END IF
      
      XKUTOP = XKAPPA*UTOP(IJ)
      
      USTOLD = (1 - IUSFG)*UTOP(IJ)*SQRT(MIN(ACD + BCD*UTOP(IJ), CDMAX)) + IUSFG*USTAR(IJ)
      TAUOLD = MAX(USTOLD**2, TAUWEFF)
      USTAR(IJ) = SQRT(TAUOLD)
      USTM1 = 1.0_JWRB / MAX(USTAR(IJ), EPSUS)
      
      DO ITER=1,NITER
        X = MAX(TAUWACT / TAUOLD, XMIN)
        Z0CH = ALPHAOG*TAUOLD / SQRT(1.0_JWRB - X)
        Z0VIS = RNUM_D*USTM1
        Z0TOT = Z0CH + Z0VIS
        
        XOLOGZ0 = 1.0_JWRB / (XLOGXL - LOG(Z0TOT))
        F = USTAR(IJ) - XKUTOP*XOLOGZ0
        ZZ = USTM1*(Z0CH*(2.0_JWRB - TWOXMP1*X) / (1.0_JWRB - X) - Z0VIS) / Z0TOT
        DELF = 1.0_JWRB - XKUTOP*XOLOGZ0**2*ZZ
        
        IF (DELF /= 0.0_JWRB) USTAR(IJ) = USTAR(IJ) - F / DELF
        TAUNEW = MAX(USTAR(IJ)**2, TAUWEFF)
        USTAR(IJ) = SQRT(TAUNEW)
        IF (TAUNEW == TAUOLD) EXIT
        USTM1 = 1.0_JWRB / MAX(USTAR(IJ), EPSUS)
        TAUOLD = TAUNEW
      END DO
      
      Z0(IJ) = Z0CH
      Z0B(IJ) = ALPHAOG*TAUOLD
      CHRNCK(IJ) = MAX(G_D*Z0(IJ)*USTM1**2, ALPHAMIN_D)
      
      
    END IF
    
  END SUBROUTINE TAUT_Z0_CUF_PARAMETRISE
END MODULE TAUT_Z0_CUF_PARAMETRISE_MOD
