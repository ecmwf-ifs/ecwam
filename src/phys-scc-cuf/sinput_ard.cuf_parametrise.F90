! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SINPUT_ARD_CUF_PARAMETRISE_MOD
  !CONTAINED SUBROUTINES:
  ! - WSIGSTAR
  ! - SINPUT_ARD
  ! - SINPUT_JAN
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE WSIGSTAR_CUF_PARAMETRISE (WSWAVE, UFRIC, Z0M, WSTAR, SIG_N)
    ! ----------------------------------------------------------------------
    
    !**** *WSIGSTAR* - COMPUTATION OF THE RELATIVE STANDARD DEVIATION OF USTAR.
    
    !*    PURPOSE.
    !     ---------
    
    !     COMPUTES THE STANDARD DEVIATION OF USTAR DUE TO SMALL SCALE GUSTINESS
    !     RELATIVE TO USTAR
    
    !**   INTERFACE.
    !     ----------
    
    !     *CALL* *WSIGSTAR (KIJS, KIJL, WSWAVE, UFRIC, Z0M, WSTAR, SIG_N)
    !             *KIJS*   - INDEX OF FIRST GRIDPOINT.
    !             *KIJL*   - INDEX OF LAST GRIDPOINT.
    !             *WSWAVE* - 10M WIND SPEED (m/s).
    !             *UFRIC*  - NEW FRICTION VELOCITY IN M/S.
    !             *Z0M*    - ROUGHNESS LENGTH IN M.
    !             *WSTAR*  - FREE CONVECTION VELOCITY SCALE (M/S).
    !             *SIG_N*  - ESTINATED RELATIVE STANDARD DEVIATION OF USTAR.
    
    !     METHOD.
    !     -------
    
    !     USE PANOFSKY (1991) TO EXPRESS THE STANDARD DEVIATION OF U10 IN TERMS
    !     USTAR AND  w* THE CONVECTIVE VELOCITY SCALE.
    !     (but with the background gustiness set to 0.)
    !     and USTAR=SQRT(Cd)*U10 to DERIVE THE STANDARD DEVIATION OF USTAR.
    !     WITH CD=A+B*U10 (see below).
    
    !     REFERENCE.
    !     ----------
    
    !     SEE SECTION 3.2.1 OF THE WAM DOCUMENTATION.
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWCOUP, ONLY: LLGCBZ0_D
    USE YOWPCONS, ONLY: G_D, EPSUS, ACDLIN, BCDLIN
    USE YOWPHYS, ONLY: XKAPPA, RNUM_D, ALPHAMIN_D, ALPHAMAX
    USE YOWWIND, ONLY: WSPMIN_D
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    REAL(KIND=JWRB), INTENT(IN) :: WSWAVE, UFRIC, Z0M, WSTAR
    REAL(KIND=JWRB), INTENT(OUT) :: SIG_N
    
    REAL(KIND=JWRB), PARAMETER :: BG_GUST = 0.0_JWRB    ! NO BACKGROUND GUSTINESS (S0 12. IS NOT USED)
    REAL(KIND=JWRB), PARAMETER :: ONETHIRD = 1.0_JWRB / 3.0_JWRB
    REAL(KIND=JWRB), PARAMETER :: SIG_NMAX = 0.9_JWRB    ! MAX OF RELATIVE STANDARD DEVIATION OF USTAR
    
    REAL(KIND=JWRB), PARAMETER :: LOG10 = LOG(10.0_JWRB)
    REAL(KIND=JWRB), PARAMETER :: C1 = 1.03E-3_JWRB
    REAL(KIND=JWRB), PARAMETER :: C2 = 0.04E-3_JWRB
    REAL(KIND=JWRB), PARAMETER :: P1 = 1.48_JWRB
    REAL(KIND=JWRB), PARAMETER :: P2 = -0.21_JWRB
    
!$loki routine seq
    REAL(KIND=JWRB) :: ZCHAR, C_D, DC_DDU, SIG_CONV
    REAL(KIND=JWRB) :: XKAPPAD, U10, C2U10P1, U10P2
    REAL(KIND=JWRB) :: BCD, U10M1, ZN, Z0VIS
    
    
    ! ----------------------------------------------------------------------
    
    
    IF (LLGCBZ0_D) THEN
      ZN = RNUM_D
      
      U10M1 = 1.0_JWRB / MAX(WSWAVE, WSPMIN_D)
      ! CHARNOCK:
      Z0VIS = ZN / MAX(UFRIC, EPSUS)
      ZCHAR = G_D*(Z0M - Z0VIS) / MAX(UFRIC**2, EPSUS)
      ZCHAR = MAX(MIN(ZCHAR, ALPHAMAX), ALPHAMIN_D)
      
      BCD = BCDLIN*SQRT(ZCHAR)
      C_D = ACDLIN + BCD*WSWAVE
      DC_DDU = BCD
      SIG_CONV = 1.0_JWRB + 0.5_JWRB*WSWAVE / C_D*DC_DDU
      SIG_N = MIN(SIG_NMAX, SIG_CONV*U10M1*(BG_GUST*UFRIC**3 + 0.5_JWRB*XKAPPA*WSTAR**3)**ONETHIRD)
    ELSE
      ZN = 0.0_JWRB
      
      !!! for consistency I have kept the old method, even though the new method above could be used,
      !!! but until LLGCBZ0 is the default, keep the old scheme whe it is not...
      !
      !       IN THE FOLLOWING U10 IS ESTIMATED ASSUMING EVERYTHING IS
      !       BASED ON U*
      !
      XKAPPAD = 1.0_JWRB / XKAPPA
      U10 = UFRIC*XKAPPAD*(LOG10 - LOG(Z0M))
      U10 = MAX(U10, WSPMIN_D)
      U10M1 = 1.0_JWRB / U10
      C2U10P1 = C2*U10**P1
      U10P2 = U10**P2
      C_D = (C1 + C2U10P1)*U10P2
      DC_DDU = (P2*C1 + (P1 + P2)*C2U10P1)*U10P2*U10M1
      SIG_CONV = 1.0_JWRB + 0.5_JWRB*U10 / C_D*DC_DDU
      SIG_N = MIN(SIG_NMAX, SIG_CONV*U10M1*(BG_GUST*UFRIC**3 + 0.5_JWRB*XKAPPA*WSTAR**3)**ONETHIRD)
    END IF
    
    
  END SUBROUTINE WSIGSTAR_CUF_PARAMETRISE
  ATTRIBUTES(DEVICE) SUBROUTINE SINPUT_ARD_CUF_PARAMETRISE (NGST, LLSNEG, KIJS, KIJL, FL1, WAVNUM, CINV, XK2CG, WDWAVE, WSWAVE,  &
  & UFRIC, Z0M, COSWDIF, SINWDIF2, RAORW, WSTAR, RNFAC, FLD, SL, SPOS, XLLWS, IJ)
    ! ----------------------------------------------------------------------
    
    !**** *SINPUT_ARD* - COMPUTATION OF INPUT SOURCE FUNCTION.
    
    
    !*    PURPOSE.
    !     ---------
    
    !       COMPUTE THE WIND INPUT SOURCE TRERM BASED ON ARDHUIN ET AL. 2010.
    
    !       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
    !       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
    !       INPUT SOURCE FUNCTION.
    !
    !       GUSTINESS IS INTRODUCED FOLL0WING THE APPROACH OF JANSSEN(1986),
    !       USING A GAUSS-HERMITE APPROXIMATION SUGGESTED BY MILES(1997).
    !       IN THE PRESENT VERSION ONLY TWO HERMITE POLYNOMIALS ARE UTILISED
    !       IN THE EVALUATION OF THE PROBABILITY INTEGRAL. EXPLICITELY ONE THEN
    !       FINDS:
    !
    !             <GAMMA(X)> = 0.5*( GAMMA(X(1+SIG)) + GAMMA(X(1-SIG)) )
    !
    !       WHERE X IS THE FRICTION VELOCITY AND SIG IS THE RELATIVE GUSTINESS
    !       LEVEL.
    
    !**   INTERFACE.
    !     ----------
    
    !     *CALL* *SINPUT_ARD (NGST, LLSNEG, KIJS, KIJL, FL1,
    !    &                    WAVNUM, CINV, XK2CG,
    !    &                    WSWAVE, WDWAVE, UFRIC, Z0M,
    !    &                    COSWDIF, SINWDIF2,
    !    &                    RAORW, WSTAR, RNFAC,
    !    &                    FLD, SL, SPOS, XLLWS)
    !         *NGST* - IF = 1 THEN NO GUSTINESS PARAMETERISATION
    !                - IF = 2 THEN GUSTINESS PARAMETERISATION
    !         *LLSNEG- IF TRUE THEN THE NEGATIVE SINPUT (SWELL DAMPING) WILL BE COMPUTED
    !         *KIJS* - INDEX OF FIRST GRIDPOINT.
    !         *KIJL* - INDEX OF LAST GRIDPOINT.
    !          *FL1* - SPECTRUM.
    !       *WAVNUM* - WAVE NUMBER.
    !         *CINV* - INVERSE PHASE VELOCITY.
    !       *XK2CG*  - (WAVNUM)**2 * GROUP SPPED.
    !       *WDWAVE* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
    !                  NOTATION (POINTING ANGLE OF WIND VECTOR,
    !                  CLOCKWISE FROM NORTH).
    !        *UFRIC* - NEW FRICTION VELOCITY IN M/S.
    !        *Z0M* - ROUGHNESS LENGTH IN M.
    !      *COSWDIF* - COS(TH(K)-WDWAVE(IJ))
    !     *SINWDIF2* - SIN(TH(K)-WDWAVE(IJ))**2
    !        *RAORW* - RATIO AIR DENSITY TO WATER DENSITY.
    !        *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
    !        *RNFAC* - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
    !          *FLD* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
    !           *SL* - TOTAL SOURCE FUNCTION ARRAY.
    !         *SPOS* - POSITIVE SOURCE FUNCTION ARRAY.
    !       *XLLWS*  - = 1 WHERE SINPUT IS POSITIVE
    
    !     METHOD.
    !     -------
    
    !       SEE REFERENCE.
    
    
    ! ----------------------------------------------------------------------
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWSTAT, ONLY: IDAMPING_D
    USE YOWCOUP, ONLY: LLCAPCHNK_D, LLNORMAGAM_D
    USE YOWFRED, ONLY: FR_D, TH_D, DFIM_D, COSTH_D, SINTH_D, ZPIFR_D, DELTH_D
    USE YOWPARAM, ONLY: NANG_D, NFRE_D, NANG_PARAM
    USE YOWPCONS, ONLY: G_D, GM1_D, EPSMIN, EPSUS, ZPI_D
    USE YOWPHYS, ONLY: ZALP_D, TAUWSHELTER_D, XKAPPA, BETAMAXOXKAPPA2_D, RNU_D, RNUM_D, SWELLF, SWELLF2, SWELLF3, SWELLF4,  &
    & SWELLF5_D, SWELLF6, SWELLF7, SWELLF7M1, Z0RAT_D, Z0TUBMAX_D, ABMIN, ABMAX
    USE YOWTEST, ONLY: IU06
    USE YOWTABL, ONLY: IAB, SWELLFT_D
    USE YOWSTAT, ONLY: IDAMPING_D
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NGST
    LOGICAL, VALUE, INTENT(IN) :: LLSNEG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NFRE_LOKI_PARAM) :: WAVNUM, CINV, XK2CG
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: WDWAVE, WSWAVE, UFRIC, Z0M
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: RAORW, RNFAC
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: WSTAR
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(NANG_LOKI_PARAM) :: COSWDIF, SINWDIF2
    REAL(KIND=JWRB), INTENT(OUT), DEVICE, DIMENSION(NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FLD, SL, SPOS
    REAL(KIND=JWRB), INTENT(OUT), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: XLLWS
    
    
    INTEGER(KIND=JWIM) :: K, M, IND, IGST
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    REAL(KIND=JWRB) :: CONSTN
    REAL(KIND=JWRB) :: AVG_GST, ABS_TAUWSHELTER
    REAL(KIND=JWRB) :: CONST1
    REAL(KIND=JWRB) :: ZNZ
    REAL(KIND=JWRB) :: X1, X2, ZLOG, ZLOG1, ZLOG2, ZLOG2X, XV1, XV2, ZBETA1, ZBETA2
    REAL(KIND=JWRB) :: XI, X, DELI1, DELI2
    REAL(KIND=JWRB) :: FU, FUD, NU_AIR, SMOOTH, HFTSWELLF6, Z0TUB
    REAL(KIND=JWRB) :: FAC_NU_AIR, FACM1_NU_AIR
    REAL(KIND=JWRB) :: ARG, DELABM1
    REAL(KIND=JWRB) :: TAUPX, TAUPY
    REAL(KIND=JWRB) :: DSTAB2
    
    REAL(KIND=JWRB) :: SIG2, COEF, COEF5, DFIM_SIG2, COSLP
    
    REAL(KIND=JWRB) :: XNGAMCONST
    REAL(KIND=JWRB) :: CONSTF, CONST11, CONST22
    REAL(KIND=JWRB) :: Z0VIS, Z0NOZ, FWW
    REAL(KIND=JWRB) :: PVISC, PTURB
    REAL(KIND=JWRB) :: ZCN
    REAL(KIND=JWRB) :: SIG_N, UORBT, AORB, TEMP, RE, RE_C, ZORB
    REAL(KIND=JWRB) :: CNSN, SUMF, SUMFSIN2
    REAL(KIND=JWRB) :: CSTRNFAC
    REAL(KIND=JWRB) :: FLP_AVG, SLP_AVG
    REAL(KIND=JWRB) :: ROGOROAIR, AIRD_PVISC
    REAL(KIND=JWRB) :: DSTAB1, TEMP1, TEMP2
    
    REAL(KIND=JWRB), DIMENSION(2) :: XSTRESS, YSTRESS, FLP, SLP
    REAL(KIND=JWRB), DIMENSION(2) :: USG2, TAUX, TAUY, USTP, USTPM1, USDIRP, UCN
    REAL(KIND=JWRB), DIMENSION(2) :: UCNZALPD
    REAL(KIND=JWRB), DIMENSION(2) :: GAMNORMA    ! ! RENORMALISATION FACTOR OF THE GROWTH RATE
    REAL(KIND=JWRB), DIMENSION(2, NANG_PARAM) :: GAM0, DSTAB
    
    LOGICAL :: LTAUWSHELTER
    ! ----------------------------------------------------------------------
    
    
    AVG_GST = 1.0_JWRB / NGST
    CONST1 = BETAMAXOXKAPPA2_D
    CONSTN = DELTH_D / (XKAPPA*ZPI_D)
    
    ABS_TAUWSHELTER = ABS(TAUWSHELTER_D)
    IF (ABS_TAUWSHELTER == 0.0_JWRB) THEN
      LTAUWSHELTER = .false.
    ELSE
      LTAUWSHELTER = .true.
    END IF
    
    IF (NGST > 1) THEN
      CALL WSIGSTAR_CUF_PARAMETRISE(WSWAVE(IJ), UFRIC(IJ), Z0M(IJ), WSTAR(IJ), SIG_N)
    END IF
    
    
    IF (LLNORMAGAM_D) THEN
      CSTRNFAC = CONSTN*RNFAC / RAORW
    END IF
    
    
    !     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
    
    ! ----------------------------------------------------------------------
    IF (LLSNEG) THEN
      !!!!  only for the negative sinput
      NU_AIR = RNU_D
      FACM1_NU_AIR = 4.0_JWRB / NU_AIR
      
      FAC_NU_AIR = RNUM_D
      
      FU = ABS(SWELLF3)
      FUD = SWELLF2
      DELABM1 = REAL(IAB) / (ABMAX - ABMIN)
      
      
      !       computation of Uorb and Aorb
      UORBT = EPSMIN
      AORB = EPSMIN
      
      DO M=1,NFRE_D
        SIG2 = ZPIFR_D(M)**2
        DFIM_SIG2 = DFIM_D(M)*SIG2
        
        K = 1
        TEMP = FL1(IJ, K, M)
        DO K=2,NANG_D
          TEMP = TEMP + FL1(IJ, K, M)
        END DO
        
        UORBT = UORBT + DFIM_SIG2*TEMP
        AORB = AORB + DFIM_D(M)*TEMP
      END DO
      
      UORBT = 2.0_JWRB*SQRT(UORBT)        ! this is the significant orbital amplitude
      AORB = 2.0_JWRB*SQRT(AORB)        ! this 1/2 Hs
      RE = FACM1_NU_AIR*UORBT*AORB        ! this is the Reynolds number
      Z0VIS = FAC_NU_AIR / MAX(UFRIC(IJ), 0.0001_JWRB)
      Z0TUB = Z0RAT_D*MIN(Z0TUBMAX_D, Z0M(IJ))
      Z0NOZ = MAX(Z0VIS, Z0TUB)
      ZORB = AORB / Z0NOZ
      
      !         compute fww
      XI = (LOG10(MAX(ZORB, 3.0_JWRB)) - ABMIN)*DELABM1
      IND = MIN(IAB - 1, INT(XI))
      DELI1 = MIN(1.0_JWRB, XI - REAL(IND, kind=JWRB))
      DELI2 = 1.0_JWRB - DELI1
      FWW = SWELLFT_D(IND)*DELI2 + SWELLFT_D(IND + 1)*DELI1
      TEMP2 = FWW*UORBT
      
      !       Define the critical Reynolds number
      IF (SWELLF6 == 1.0_JWRB) THEN
        RE_C = SWELLF4
      ELSE
        HFTSWELLF6 = 1.0_JWRB - SWELLF6
        RE_C = SWELLF4*(2.0_JWRB / AORB)**HFTSWELLF6
      END IF
      
      !       Swell damping weight between viscous and turbulent boundary layer
      IF (SWELLF7 > 0.0_JWRB) THEN
        SMOOTH = 0.5_JWRB*TANH((RE - RE_C)*SWELLF7M1)
        PTURB = 0.5_JWRB + SMOOTH
        PVISC = 0.5_JWRB - SMOOTH
      ELSE
        IF (RE <= RE_C) THEN
          PTURB = 0.0_JWRB
          PVISC = 0.5_JWRB
        ELSE
          PTURB = 0.5_JWRB
          PVISC = 0.0_JWRB
        END IF
      END IF
      
      AIRD_PVISC = PVISC*RAORW
      
    END IF
    
    
    
    ! Initialisation
    
    IF (NGST == 1) THEN
      USTP(1) = UFRIC(IJ)
    ELSE
      USTP(1) = UFRIC(IJ)*(1.0_JWRB + SIG_N)
      USTP(2) = UFRIC(IJ)*(1.0_JWRB - SIG_N)
    END IF
    
    DO IGST=1,NGST
      USTPM1(IGST) = 1.0_JWRB / MAX(USTP(IGST), EPSUS)
    END DO
    
    IF (LTAUWSHELTER) THEN
      DO IGST=1,NGST
        XSTRESS(IGST) = 0.0_JWRB
        YSTRESS(IGST) = 0.0_JWRB
        USG2(IGST) = USTP(IGST)**2
        TAUX(IGST) = USG2(IGST)*SIN(WDWAVE(IJ))
        TAUY(IGST) = USG2(IGST)*COS(WDWAVE(IJ))
      END DO
      
      ROGOROAIR = G_D / RAORW
    END IF
    
    
    !*    2. MAIN LOOP OVER FREQUENCIES.
    !        ---------------------------
    
    IF (.not.LLNORMAGAM_D) THEN
      DO IGST=1,NGST
        GAMNORMA(IGST) = 1.0_JWRB
      END DO
    END IF
    
    IF (.not.LLSNEG) THEN
      DO K=1,NANG_D
        DO IGST=1,NGST
          DSTAB(IGST, K) = 0.0_JWRB
        END DO
      END DO
    END IF
    
    DO M=1,NFRE_D
      
      IF (LTAUWSHELTER) THEN
        DO IGST=1,NGST
          TAUPX = TAUX(IGST) - ABS_TAUWSHELTER*XSTRESS(IGST)
          TAUPY = TAUY(IGST) - ABS_TAUWSHELTER*YSTRESS(IGST)
          USDIRP(IGST) = ATAN2(TAUPX, TAUPY)
          USTP(IGST) = (TAUPX**2 + TAUPY**2)**0.25_JWRB
          USTPM1(IGST) = 1.0_JWRB / MAX(USTP(IGST), EPSUS)
        END DO
        
        CONSTF = ROGOROAIR*CINV(IJ, M)*DFIM_D(M)
      END IF
      
      
      !*      PRECALCULATE FREQUENCY DEPENDENCE.
      !       ----------------------------------
      
      DO IGST=1,NGST
        UCN(IGST) = USTP(IGST)*CINV(IJ, M)
        UCNZALPD(IGST) = XKAPPA / (UCN(IGST) + ZALP_D)
      END DO
      ZCN = LOG(WAVNUM(IJ, M)*Z0M(IJ))
      CNSN = ZPIFR_D(M)*CONST1*RAORW
      
      !*    2.1 LOOP OVER DIRECTIONS.
      !         ---------------------
      
      DO K=1,NANG_D
        XLLWS(IJ, K, M) = 0.0_JWRB
      END DO
      
      IF (LLSNEG) THEN
        !       SWELL DAMPING:
        
        SIG2 = ZPIFR_D(M)**2
        DFIM_SIG2 = DFIM_D(M)*SIG2
        
        COEF = -SWELLF*16._JWRB*SIG2 / G_D
        COEF5 = -SWELLF5_D*2._JWRB*SQRT(2._JWRB*NU_AIR*ZPIFR_D(M))
        
        DSTAB1 = COEF5*AIRD_PVISC*WAVNUM(IJ, M)
        TEMP1 = COEF*RAORW
      END IF
      
      DO K=1,NANG_D
        DO IGST=1,NGST
          
          SUMF = 0.0_JWRB
          SUMFSIN2 = 0.0_JWRB
          
          IF (LTAUWSHELTER) THEN
            COSLP = COS(TH_D(K) - USDIRP(IGST))
          ELSE
            COSLP = COSWDIF(K)
          END IF
          
          GAM0(IGST, K) = 0._JWRB
          IF (COSLP > 0.01_JWRB) THEN
            X = COSLP*UCN(IGST)
            ZLOG = ZCN + UCNZALPD(IGST) / COSLP
            IF (ZLOG < 0.0_JWRB) THEN
              ZLOG2X = ZLOG*ZLOG*X
              GAM0(IGST, K) = EXP(ZLOG)*ZLOG2X*ZLOG2X*CNSN
              XLLWS(IJ, K, M) = 1.0_JWRB
            END IF
          END IF
          
          IF (LLSNEG) THEN
            DSTAB2 = TEMP1*(TEMP2 + (FU + FUD*COSLP)*USTP(IGST))
            DSTAB(IGST, K) = DSTAB1 + PTURB*DSTAB2
          END IF
          
          SUMF = SUMF + GAM0(IGST, K)*FL1(IJ, K, M)
          SUMFSIN2 = SUMFSIN2 + GAM0(IGST, K)*FL1(IJ, K, M)*SINWDIF2(K)
        END DO
      END DO
      
      IF (LLNORMAGAM_D) THEN
        
        XNGAMCONST = CSTRNFAC*XK2CG(IJ, M)
        DO IGST=1,NGST
          ZNZ = XNGAMCONST*USTPM1(IGST)
          GAMNORMA(IGST) = (1.0_JWRB + ZNZ*SUMFSIN2) / (1.0_JWRB + ZNZ*SUMF)
        END DO
        
      END IF
      
      
      
      !*    2.2 UPDATE THE SHELTERING STRESS (in any),
      !         AND THEN ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
      !         ---------------------------------------------------------
      
      DO K=1,NANG_D
        
        DO IGST=1,NGST
          ! SLP: only the positive contributions
          SLP(IGST) = GAM0(IGST, K)*GAMNORMA(IGST)
          FLP(IGST) = SLP(IGST) + DSTAB(IGST, K)
        END DO
        
        DO IGST=1,NGST
          SLP(IGST) = SLP(IGST)*FL1(IJ, K, M)
        END DO
        
        IF (LTAUWSHELTER) THEN
          CONST11 = CONSTF*SINTH_D(K)
          CONST22 = CONSTF*COSTH_D(K)
          DO IGST=1,NGST
            XSTRESS(IGST) = XSTRESS(IGST) + SLP(IGST)*CONST11
            YSTRESS(IGST) = YSTRESS(IGST) + SLP(IGST)*CONST22
          END DO
        END IF
        
        IGST = 1
        SLP_AVG = SLP(IGST)
        FLP_AVG = FLP(IGST)
        DO IGST=2,NGST
          SLP_AVG = SLP_AVG + SLP(IGST)
          FLP_AVG = FLP_AVG + FLP(IGST)
        END DO
        
        SPOS(K, M) = AVG_GST*SLP_AVG
        FLD(K, M) = AVG_GST*FLP_AVG
        SL(K, M) = FLD(K, M)*FL1(IJ, K, M)
        
      END DO
      
    END DO
    ! END LOOP OVER FREQUENCIES
    
    
  END SUBROUTINE SINPUT_ARD_CUF_PARAMETRISE
  ATTRIBUTES(DEVICE) SUBROUTINE SINPUT_JAN_CUF_PARAMETRISE (NGST, LLSNEG, KIJS, KIJL, FL1, WAVNUM, CINV, XK2CG, WSWAVE, UFRIC,  &
  & Z0M, COSWDIF, SINWDIF2, RAORW, WSTAR, RNFAC, FLD, SL, SPOS, XLLWS, IJ)
    ! ----------------------------------------------------------------------
    
    !**** *SINPUT_JAN* - COMPUTATION OF INPUT SOURCE FUNCTION.
    
    !     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
    
    !     OPTIMIZED BY : H. GUENTHER
    
    !     MODIFIED BY :
    !       J-R BIDLOT NOVEMBER 1995
    !       J-R BIDLOT FEBRUARY 1996-97
    !       J-R BIDLOT FEBRUARY 1999 : INTRODUCE ICALL AND NCALL
    !       P.A.E.M. JANSSEN MAY 2000 : INTRODUCE GUSTINESS
    !       J-R BIDLOT FEBRUARY 2001 : MAKE IT FULLY IMPLICIT BY ONLY
    !                                  USING NEW STRESS AND ROUGHNESS.
    !       S. ABDALLA OCTOBER 2001:  INTRODUCTION OF VARIABLE AIR
    !                                 DENSITY AND STABILITY-DEPENDENT
    !                                 WIND GUSTINESS
    !       P.A.E.M. JANSSEN OCTOBER 2008: INTRODUCE DAMPING WHEN WAVES ARE
    !                                      RUNNING FASTER THAN THE WIND.
    !       J-R BIDLOT JANUARY 2013: SHALLOW WATER FORMULATION.
    
    !*    PURPOSE.
    !     ---------
    
    !       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
    !       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
    !       INPUT SOURCE FUNCTION.
    !
    !       GUSTINESS IS INTRODUCED FOLL0WING THE APPROACH OF JANSSEN(1986),
    !       USING A GAUSS-HERMITE APPROXIMATION SUGGESTED BY MILES(1997).
    !       IN THE PRESENT VERSION ONLY TWO HERMITE POLYNOMIALS ARE UTILISED
    !       IN THE EVALUATION OF THE PROBABILITY INTEGRAL. EXPLICITELY ONE THEN
    !       FINDS:
    !
    !             <GAMMA(X)> = 0.5*( GAMMA(X(1+SIG)) + GAMMA(X(1-SIG)) )
    !
    !       WHERE X IS THE FRICTION VELOCITY AND SIG IS THE RELATIVE GUSTINESS
    !       LEVEL.
    
    !**   INTERFACE.
    !     ----------
    
    !     *CALL* *SINPUT_JAN (NGST, LLSNEG, KIJS, KIJL, FL1,
    !    &                    WAVNUM, CINV, XK2CG,
    !    &                    WDWAVE, WSWAVE, UFRIC, Z0M,
    !    &                    COSWDIF, SINWDIF2,
    !    &                    RAORW, WSTAR, RNFAC,
    !    &                    FLD, SL, SPOS, XLLWS)
    !         *NGST* - IF = 1 THEN NO GUSTINESS PARAMETERISATION
    !                - IF = 2 THEN GUSTINESS PARAMETERISATION
    !         *LLSNEG- IF TRUE THEN THE NEGATIVE SINPUT (SWELL DAMPING) WILL BE COMPUTED
    !         *KIJS* - INDEX OF FIRST GRIDPOINT.
    !         *KIJL* - INDEX OF LAST GRIDPOINT.
    !          *FL1* - SPECTRUM.
    !       *WAVNUM* - WAVE NUMBER.
    !         *CINV* - INVERSE PHASE VELOCITY.
    !       *XK2CG*  - (WAVNUM)**2 * GROUP SPPED.
    !       *WDWAVE* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
    !                  NOTATION (POINTING ANGLE OF WIND VECTOR,
    !                  CLOCKWISE FROM NORTH).
    !        *UFRIC* - FRICTION VELOCITY IN M/S.
    !        *Z0M*   - ROUGHNESS LENGTH IN M.
    !      *COSWDIF* - COS(TH(K)-WDWAVE(IJ))
    !     *SINWDIF2* - SIN(TH(K)-WDWAVE(IJ))**2
    !        *RAORW* - RATIO AIR DENSITY TO WATER DENSITY
    !        *RNFAC* - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
    !        *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
    !          *FLD* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
    !           *SL* - TOTAL SOURCE FUNCTION ARRAY.
    !         *SPOS* - ONLY POSITIVE PART OF INPUT SOURCE FUNCTION ARRAY.
    !        *XLLWS* - 1 WHERE SINPUT IS POSITIVE.
    
    
    !     METHOD.
    !     -------
    
    !       SEE REFERENCE.
    
    !     EXTERNALS.
    !     ----------
    
    !       WSIGSTAR.
    
    !     MODIFICATIONS
    !     -------------
    
    !     - REMOVAL OF CALL TO CRAY SPECIFIC FUNCTIONS EXPHF AND ALOGHF
    !       BY THEIR STANDARD FORTRAN EQUIVALENT EXP and ALOGHF
    !     - MODIFIED TO MAKE INTEGRATION SCHEME FULLY IMPLICIT
    !     - INTRODUCTION OF VARIABLE AIR DENSITY
    !     - INTRODUCTION OF WIND GUSTINESS
    
    !     REFERENCE.
    !     ----------
    
    !       P. JANSSEN, J.P.O., 1989.
    !       P. JANSSEN, J.P.O., 1991
    
    ! ----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWCOUP, ONLY: LLNORMAGAM_D
    USE YOWFRED, ONLY: ZPIFR_D, DELTH_D, TH_D
    USE YOWFRED, ONLY: FR_D, TH_D, ZPIFR_D
    USE YOWPARAM, ONLY: NANG_D, NFRE_D, NANG_PARAM
    USE YOWPCONS, ONLY: G_D, GM1_D, ZPI_D, EPSUS
    USE YOWPHYS, ONLY: ZALP_D, XKAPPA, BETAMAXOXKAPPA2_D
    USE YOWSTAT, ONLY: IDAMPING_D
    USE YOWTEST, ONLY: IU06
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: NGST
    LOGICAL, VALUE, INTENT(IN) :: LLSNEG
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NFRE_LOKI_PARAM) :: WAVNUM, CINV, XK2CG
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: WSWAVE, UFRIC, Z0M
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(NANG_LOKI_PARAM) :: COSWDIF, SINWDIF2
    REAL(KIND=JWRB), INTENT(IN), DEVICE :: RAORW, RNFAC
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: WSTAR
    REAL(KIND=JWRB), INTENT(OUT), DEVICE, DIMENSION(NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FLD, SL, SPOS
    REAL(KIND=JWRB), INTENT(OUT), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: XLLWS
    
    
    INTEGER(KIND=JWIM) :: IG, K, M
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    INTEGER(KIND=JWIM) :: IGST
    
    REAL(KIND=JWRB) :: CONST1, CONST3, XKAPPAD
    REAL(KIND=JWRB) :: CONSTN
    REAL(KIND=JWRB) :: ZNZ
    REAL(KIND=JWRB) :: X, ZLOG, ZLOG2X, ZBETA
    REAL(KIND=JWRB) :: TEMPD
    
    REAL(KIND=JWRB), DIMENSION(2) :: WSIN
    REAL(KIND=JWRB) :: ZTANHKD
    REAL(KIND=JWRB) :: SIG_N
    REAL(KIND=JWRB) :: CNSN
    REAL(KIND=JWRB) :: SUMF, SUMFSIN2
    REAL(KIND=JWRB) :: CSTRNFAC
    REAL(KIND=JWRB) :: UFAC1, UFAC2
    REAL(KIND=JWRB), DIMENSION(2) :: GAMNORMA    ! ! RENORMALISATION FACTOR OF THE GROWTH RATE
    REAL(KIND=JWRB), DIMENSION(2) :: SIGDEV, US, Z0, UCN, ZCN
    REAL(KIND=JWRB), DIMENSION(2) :: USTPM1
    REAL(KIND=JWRB), DIMENSION(2) :: XVD, UCND, CONST3_UCN2
    REAL(KIND=JWRB), DIMENSION(2, NANG_PARAM) :: GAM0
    
    LOGICAL :: LZ
    
    ! ----------------------------------------------------------------------
    
    
    CONST1 = BETAMAXOXKAPPA2_D
    CONST3 = 2.0_JWRB*XKAPPA / CONST1      ! SEE IDAMPING
    XKAPPAD = 1.E0_JWRB / XKAPPA
    
    CONST3 = IDAMPING_D*CONST3
    
    CONSTN = DELTH_D / (XKAPPA*ZPI_D)
    
    !     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
    IF (NGST > 1) THEN
      CALL WSIGSTAR_CUF_PARAMETRISE(WSWAVE(IJ), UFRIC(IJ), Z0M(IJ), WSTAR(IJ), SIG_N)
    END IF
    
    !     DEFINE WHERE SINPUT WILL BE EVALUATED IN RELATIVE TERM WRT USTAR
    !     DEFINE ALSO THE RELATIVE WEIGHT OF EACH.
    
    IF (NGST == 1) THEN
      WSIN(1) = 1.0_JWRB
      SIGDEV(1) = 1.0_JWRB
    ELSE
      WSIN(1) = 0.5_JWRB
      WSIN(2) = 0.5_JWRB
      SIGDEV(1) = 1.0_JWRB - SIG_N
      SIGDEV(2) = 1.0_JWRB + SIG_N
    END IF
    
    
    IF (NGST == 1) THEN
      US(1) = UFRIC(IJ)
      Z0(1) = Z0M(IJ)
    ELSE
      DO IGST=1,NGST
        US(IGST) = UFRIC(IJ)*SIGDEV(IGST)
        Z0(IGST) = Z0M(IJ)
      END DO
    END IF
    
    DO IGST=1,NGST
      USTPM1(IGST) = 1.0_JWRB / MAX(US(IGST), EPSUS)
    END DO
    
    ! ----------------------------------------------------------------------
    
    !*    2. LOOP OVER FREQUENCIES.
    !        ----------------------
    
    DO M=1,NFRE_D
      
      !*      PRECALCULATE FREQUENCY DEPENDENCE.
      !       ----------------------------------
      
      ZTANHKD = ZPIFR_D(M)**2 / (G_D*WAVNUM(IJ, M))
      CNSN = CONST1*ZPIFR_D(M)*ZTANHKD*RAORW
      
      DO IGST=1,NGST
        UCN(IGST) = US(IGST)*CINV(IJ, M) + ZALP_D
        CONST3_UCN2(IGST) = CONST3*UCN(IGST)**2
        UCND(IGST) = 1.0_JWRB / UCN(IGST)
        ZCN(IGST) = LOG(WAVNUM(IJ, M)*Z0(IGST))
        XVD(IGST) = 1.0_JWRB / (-US(IGST)*XKAPPAD*ZCN(IGST)*CINV(IJ, M))
      END DO
      
      !*    2.1 LOOP OVER DIRECTIONS.
      !         ---------------------
      
      !       WIND INPUT:
      DO K=1,NANG_D
        XLLWS(IJ, K, M) = 0.0_JWRB
        
        DO IGST=1,NGST
          
          IF (COSWDIF(K) > 0.01_JWRB) THEN
            LZ = .true.
            TEMPD = XKAPPA / COSWDIF(K)
          ELSE
            LZ = .false.
            TEMPD = XKAPPA
          END IF
          
          GAM0(IGST, K) = 0.0_JWRB
          IF (LZ) THEN
            ZLOG = ZCN(IGST) + TEMPD*UCND(IGST)
            IF (ZLOG < 0.0_JWRB) THEN
              X = COSWDIF(K)*UCN(IGST)
              ZLOG2X = ZLOG*ZLOG*X
              GAM0(IGST, K) = ZLOG2X*ZLOG2X*EXP(ZLOG)*CNSN
              XLLWS(IJ, K, M) = 1.0_JWRB
            END IF
          END IF
        END DO
        
      END DO
      
      
      IF (LLNORMAGAM_D) THEN
        
        SUMF = 0.0_JWRB
        SUMFSIN2 = 0.0_JWRB
        DO K=1,NANG_D
          DO IGST=1,NGST
            SUMF = SUMF + GAM0(IGST, K)*FL1(IJ, K, M)
            SUMFSIN2 = SUMFSIN2 + GAM0(IGST, K)*FL1(IJ, K, M)*SINWDIF2(K)
          END DO
          
          CSTRNFAC = CONSTN*RNFAC / RAORW
          ZNZ = CSTRNFAC*XK2CG(IJ, M)*USTPM1(IGST)
          GAMNORMA(IGST) = (1.0_JWRB + ZNZ*SUMFSIN2) / (1.0_JWRB + ZNZ*SUMF)
          
        END DO
      ELSE
        DO IGST=1,NGST
          GAMNORMA(IGST) = 1.0_JWRB
        END DO
      END IF
      
      DO K=1,NANG_D
        UFAC1 = WSIN(1)*GAM0(1, K)*GAMNORMA(1)
        DO IGST=2,NGST
          UFAC1 = UFAC1 + WSIN(IGST)*GAM0(IGST, K)*GAMNORMA(IGST)
        END DO
        
        UFAC2 = 0.0_JWRB
        IF (LLSNEG) THEN
          !         SWELL DAMPING:
          ZBETA = CONST3_UCN2(1)*(COSWDIF(K) - XVD(1))
          UFAC2 = WSIN(1)*ZBETA
          DO IGST=2,NGST
            ZBETA = CONST3_UCN2(IGST)*(COSWDIF(K) - XVD(IGST))
            UFAC2 = UFAC2 + WSIN(IGST)*ZBETA
          END DO
        END IF
        
        FLD(K, M) = UFAC1 + UFAC2*CNSN
        SPOS(K, M) = UFAC1*FL1(IJ, K, M)
        SL(K, M) = FLD(K, M)*FL1(IJ, K, M)
      END DO
      
      !*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
      !         ------------------------------------------------
      
    END DO
    
    
  END SUBROUTINE SINPUT_JAN_CUF_PARAMETRISE
END MODULE SINPUT_ARD_CUF_PARAMETRISE_MOD
