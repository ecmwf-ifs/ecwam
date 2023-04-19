! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SINPUT_ARD_MOD
!CONTAINED SUBROUTINES:
! - WSIGSTAR
! - SINPUT_ARD
! - SINPUT_JAN
CONTAINS
      SUBROUTINE WSIGSTAR (WSWAVE, UFRIC, Z0M, WSTAR, SIG_N)
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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LLGCBZ0
      USE YOWPCONS , ONLY : G, EPSUS, ACDLIN, BCDLIN 
      USE YOWPHYS  , ONLY : XKAPPA, RNUM, ALPHAMIN, ALPHAMAX
      USE YOWWIND  , ONLY : WSPMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN) :: WSWAVE, UFRIC, Z0M, WSTAR
      REAL(KIND=JWRB), INTENT(OUT) :: SIG_N

      REAL(KIND=JWRB), PARAMETER :: BG_GUST = 0.0_JWRB ! NO BACKGROUND GUSTINESS (S0 12. IS NOT USED)
      REAL(KIND=JWRB), PARAMETER :: ONETHIRD = 1.0_JWRB/3.0_JWRB
      REAL(KIND=JWRB), PARAMETER :: SIG_NMAX = 0.9_JWRB ! MAX OF RELATIVE STANDARD DEVIATION OF USTAR 

      REAL(KIND=JWRB), PARAMETER :: LOG10 = LOG(10.0_JWRB)
      REAL(KIND=JWRB), PARAMETER :: C1 = 1.03E-3_JWRB
      REAL(KIND=JWRB), PARAMETER :: C2 = 0.04E-3_JWRB
      REAL(KIND=JWRB), PARAMETER :: P1 = 1.48_JWRB
      REAL(KIND=JWRB), PARAMETER :: P2 = -0.21_JWRB

      REAL(KIND=JWRB) :: ZCHAR, C_D, DC_DDU, SIG_CONV
      REAL(KIND=JWRB) :: XKAPPAD, U10, C2U10P1, U10P2
      REAL(KIND=JWRB) :: BCD, U10M1, ZN, Z0VIS
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WSIGSTAR',0,ZHOOK_HANDLE)

      IF (LLGCBZ0) THEN
        ZN = RNUM

        U10M1=1.0_JWRB/MAX(WSWAVE,WSPMIN)
        ! CHARNOCK:
        Z0VIS = ZN/MAX(UFRIC,EPSUS)
        ZCHAR=G*(Z0M-Z0VIS)/MAX(UFRIC**2,EPSUS)
        ZCHAR=MAX(MIN(ZCHAR,ALPHAMAX),ALPHAMIN)

        BCD = BCDLIN*SQRT(ZCHAR)
        C_D = ACDLIN + BCD * WSWAVE
        DC_DDU = BCD
        SIG_CONV = 1.0_JWRB + 0.5_JWRB*WSWAVE/C_D * DC_DDU
        SIG_N = MIN(SIG_NMAX, SIG_CONV * U10M1*(BG_GUST*UFRIC**3 + &
     &                  0.5_JWRB*XKAPPA*WSTAR**3)**ONETHIRD )
       ELSE
        ZN = 0.0_JWRB

!!! for consistency I have kept the old method, even though the new method above could be used,
!!! but until LLGCBZ0 is the default, keep the old scheme whe it is not...
!
!       IN THE FOLLOWING U10 IS ESTIMATED ASSUMING EVERYTHING IS
!       BASED ON U*
!
        XKAPPAD=1.0_JWRB/XKAPPA
        U10 = UFRIC*XKAPPAD*(LOG10-LOG(Z0M))
        U10 = MAX(U10,WSPMIN)
        U10M1=1.0_JWRB/U10
        C2U10P1=C2*U10**P1
        U10P2=U10**P2
        C_D = (C1 + C2U10P1)*U10P2
        DC_DDU = (P2*C1+(P1+P2)*C2U10P1)*U10P2*U10M1
        SIG_CONV = 1.0_JWRB + 0.5_JWRB*U10/C_D*DC_DDU
        SIG_N = MIN(SIG_NMAX, SIG_CONV * U10M1*(BG_GUST*UFRIC**3 + &
     &                  0.5_JWRB*XKAPPA*WSTAR**3)**ONETHIRD )
      ENDIF

      IF (LHOOK) CALL DR_HOOK('WSIGSTAR',1,ZHOOK_HANDLE)

      END SUBROUTINE WSIGSTAR
      SUBROUTINE SINPUT_ARD (NGST, LLSNEG, KIJS, KIJL, FL1, &
 &                     WAVNUM, CINV, XK2CG,           &
 &                     WDWAVE, WSWAVE, UFRIC, Z0M,    &
 &                     COSWDIF, SINWDIF2,             &
 &                     RAORW, WSTAR, RNFAC,           &
 &                     FLD, SL, SPOS, XLLWS)
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
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LLCAPCHNK,LLNORMAGAM
      USE YOWFRED  , ONLY : FR       ,TH       ,DFIM     ,COSTH  ,SINTH, ZPIFR, DELTH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,GM1      ,EPSMIN, EPSUS, ZPI
      USE YOWPHYS  , ONLY : ZALP     ,TAUWSHELTER, XKAPPA, BETAMAXOXKAPPA2,    &
     &                      RN1_RN, &
     &                      RNU      ,RNUM, &
     &                      SWELLF   ,SWELLF2  ,SWELLF3  ,SWELLF4  , SWELLF5, &
     &                      SWELLF6  ,SWELLF7  ,SWELLF7M1, Z0RAT   ,Z0TUBMAX , &
     &                      ABMIN  ,ABMAX
      USE YOWTEST  , ONLY : IU06
      USE YOWTABL  , ONLY : IAB      ,SWELLFT

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: NGST
      LOGICAL, INTENT(IN) :: LLSNEG
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM, CINV, XK2CG
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: WDWAVE, WSWAVE, UFRIC, Z0M
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: RAORW, WSTAR, RNFAC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG), INTENT(IN) :: COSWDIF, SINWDIF2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: FLD, SL, SPOS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      INTEGER(KIND=JWIM) :: IJ, K, M, IND, IGST

      REAL(KIND=JWRB) :: CONSTN 
      REAL(KIND=JWRB) :: AVG_GST, ABS_TAUWSHELTER 
      REAL(KIND=JWRB) :: CONST1
      REAL(KIND=JWRB) :: ZNZ
      REAL(KIND=JWRB) :: X1, X2, ZLOG, ZLOG1, ZLOG2, ZLOG2X, XV1, XV2, ZBETA1, ZBETA2
      REAL(KIND=JWRB) :: XI, X, DELI1, DELI2
      REAL(KIND=JWRB) :: FU, FUD, NU_AIR,SMOOTH, HFTSWELLF6, Z0TUB
      REAL(KIND=JWRB) :: FAC_NU_AIR, FACM1_NU_AIR
      REAL(KIND=JWRB) :: ARG, DELABM1
      REAL(KIND=JWRB) :: TAUPX, TAUPY
      REAL(KIND=JWRB) :: DSTAB2
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE) :: CONST, SIG, SIGM1, SIG2, COEF, COEF5, DFIM_SIG2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: CONSTF, CONST11, CONST22
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: Z0VIS, Z0NOZ, FWW
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: PVISC, PTURB
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ZCN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: SIG_N, UORBT, AORB, TEMP, RE, RE_C, ZORB
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: CNSN, SUMF, SUMFSIN2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: CSTRNFAC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: FLP_AVG, SLP_AVG
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ROGOROAIR, AIRD_PVISC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NGST) :: XSTRESS, YSTRESS, FLP, SLP
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NGST) :: USG2, TAUX, TAUY, USTP, USTPM1, USDIRP, UCN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NGST) :: UCNZALPD
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: XNGAMCONST
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NGST) :: GAMNORMA ! ! RENORMALISATION FACTOR OF THE GROWTH RATE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: DSTAB1, TEMP1, TEMP2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NGST) :: COSLP, GAM0, DSTAB

      LOGICAL :: LTAUWSHELTER
! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SINPUT_ARD',0,ZHOOK_HANDLE)

      AVG_GST = 1.0_JWRB/NGST
      CONST1  = BETAMAXOXKAPPA2
      CONSTN = DELTH/(XKAPPA*ZPI) 

      ABS_TAUWSHELTER=ABS(TAUWSHELTER)
      IF (ABS_TAUWSHELTER == 0.0_JWRB ) THEN
        LTAUWSHELTER = .FALSE.
      ELSE
        LTAUWSHELTER = .TRUE.
      ENDIF

      IF (LLNORMAGAM) THEN
        DO IJ=KIJS,KIJL
          CSTRNFAC(IJ) = CONSTN * RNFAC(IJ) / RAORW(IJ)
        ENDDO
        DO M=1,NFRE
          DO IJ=KIJS,KIJL
            XNGAMCONST(IJ,M) = CSTRNFAC(IJ)*XK2CG(IJ,M)
          ENDDO
        ENDDO

      ELSE
        XNGAMCONST(KIJS:KIJL,:) = 0.0_JWRB
      ENDIF


!     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
      DO IJ=KIJS,KIJL
        IF (NGST > 1) CALL WSIGSTAR (WSWAVE(IJ), UFRIC(IJ), Z0M(IJ), WSTAR(IJ), SIG_N(IJ))
      ENDDO

! ----------------------------------------------------------------------

      DO M=1,NFRE
        SIG(M) = ZPIFR(M)
        SIGM1(M) = 1.0_JWRB/SIG(M)
        SIG2(M) = SIG(M)**2
        DFIM_SIG2(M)=DFIM(M)*SIG2(M)
        CONST(M)=SIG(M)*CONST1
      ENDDO


      IF (LLSNEG) THEN
!!!!  only for the negative sinput
         NU_AIR = RNU
         FACM1_NU_AIR = 4.0_JWRB/NU_AIR

         FAC_NU_AIR = RNUM

         FU=ABS(SWELLF3)
         FUD=SWELLF2
         DELABM1= REAL(IAB)/(ABMAX-ABMIN)


!       computation of Uorb and Aorb
        DO IJ=KIJS,KIJL
          UORBT(IJ) = EPSMIN 
          AORB(IJ) = EPSMIN
        ENDDO

        DO M=1,NFRE
          COEF(M) =-SWELLF*16._JWRB*SIG2(M)/G
          COEF5(M)=-SWELLF5*2._JWRB*SQRT(2._JWRB*NU_AIR*SIG(M))

          K=1
          DO IJ=KIJS,KIJL
            TEMP(IJ) = FL1(IJ,K,M)
          ENDDO
          DO K=2,NANG
            DO IJ=KIJS,KIJL
              TEMP(IJ) = TEMP(IJ)+FL1(IJ,K,M)
            ENDDO
          ENDDO

          DO IJ=KIJS,KIJL
            UORBT(IJ) = UORBT(IJ)+DFIM_SIG2(M)*TEMP(IJ)
            AORB(IJ) = AORB(IJ)+DFIM(M)*TEMP(IJ)
          ENDDO
        ENDDO

        DO IJ=KIJS,KIJL
          UORBT(IJ) = 2.0_JWRB*SQRT(UORBT(IJ))  ! this is the significant orbital amplitude
          AORB(IJ)  = 2.0_JWRB*SQRT(AORB(IJ))   ! this 1/2 Hs
          RE(IJ)    = FACM1_NU_AIR*UORBT(IJ)*AORB(IJ) ! this is the Reynolds number 
          Z0VIS(IJ) = FAC_NU_AIR/MAX(UFRIC(IJ),0.0001_JWRB)
          Z0TUB = Z0RAT*MIN(Z0TUBMAX,Z0M(IJ))
          Z0NOZ(IJ) = MAX(Z0VIS(IJ),Z0TUB)
          ZORB(IJ)  = AORB(IJ)/Z0NOZ(IJ)

!         compute fww
          XI=(LOG10(MAX(ZORB(IJ),3.0_JWRB))-ABMIN)*DELABM1
          IND  = MIN (IAB-1, INT(XI))
          DELI1= MIN (1.0_JWRB ,XI-REAL(IND,JWRB))
          DELI2= 1.0_JWRB - DELI1
          FWW(IJ)= SWELLFT(IND)*DELI2+SWELLFT(IND+1)*DELI1
          TEMP2(IJ) = FWW(IJ)*UORBT(IJ)
        ENDDO

!       Define the critical Reynolds number
        IF ( SWELLF6 == 1.0_JWRB) THEN
          DO IJ=KIJS,KIJL
            RE_C(IJ) = SWELLF4
          ENDDO
        ELSE
          HFTSWELLF6=1.0_JWRB-SWELLF6
          DO IJ=KIJS,KIJL
            RE_C(IJ) = SWELLF4*(2.0_JWRB/AORB(IJ))**HFTSWELLF6
          ENDDO
        ENDIF

!       Swell damping weight between viscous and turbulent boundary layer
        IF (SWELLF7 > 0.0_JWRB) THEN
          DO IJ=KIJS,KIJL
            SMOOTH=0.5_JWRB*TANH((RE(IJ)-RE_C(IJ))*SWELLF7M1)
            PTURB(IJ)=0.5_JWRB+SMOOTH
            PVISC(IJ)=0.5_JWRB-SMOOTH
          ENDDO
        ELSE
          DO IJ=KIJS,KIJL
            IF (RE(IJ) <= RE_C(IJ)) THEN
              PTURB(IJ)=0.0_JWRB
              PVISC(IJ)=0.5_JWRB
            ELSE
              PTURB(IJ)=0.5_JWRB
              PVISC(IJ)=0.0_JWRB
            ENDIF
          ENDDO
        ENDIF

        DO IJ=KIJS,KIJL
          AIRD_PVISC(IJ) = PVISC(IJ)*RAORW(IJ)
        ENDDO

      ENDIF



! Initialisation

      IF (NGST == 1) THEN
        DO IJ=KIJS,KIJL
          USTP(IJ,1) = UFRIC(IJ)
        ENDDO
      ELSE IF (NGST == 2) THEN
        DO IJ=KIJS,KIJL
          USTP(IJ,1) = UFRIC(IJ)*(1.0_JWRB+SIG_N(IJ))
          USTP(IJ,2) = UFRIC(IJ)*(1.0_JWRB-SIG_N(IJ))
        ENDDO
      ELSE
         WRITE (IU06,*) '**************************************'
         WRITE (IU06,*) '*    FATAL ERROR                     *'
         WRITE (IU06,*) '*    ===========                     *'
         WRITE (IU06,*) '* IN SINPUT_ARD: NGST > 2            *'
         WRITE (IU06,*) '* NGST = ', NGST
         WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.  *'
         WRITE (IU06,*) '*                                    *'
         WRITE (IU06,*) '**************************************'
         CALL ABORT1
      ENDIF

      DO IGST=1,NGST
        DO IJ=KIJS,KIJL
          USTPM1(IJ,IGST) = 1.0_JWRB/MAX(USTP(IJ,IGST),EPSUS)
        ENDDO
      ENDDO

      IF (LTAUWSHELTER) THEN
        DO IGST=1,NGST
          DO IJ=KIJS,KIJL
            XSTRESS(IJ,IGST)=0.0_JWRB
            YSTRESS(IJ,IGST)=0.0_JWRB
            USG2(IJ,IGST)=USTP(IJ,IGST)**2
            TAUX(IJ,IGST)=USG2(IJ,IGST)*SIN(WDWAVE(IJ))
            TAUY(IJ,IGST)=USG2(IJ,IGST)*COS(WDWAVE(IJ))
          ENDDO
        ENDDO

        DO IJ=KIJS,KIJL
          ROGOROAIR(IJ) = G/RAORW(IJ)
        ENDDO

      ELSE
        DO IGST=1,NGST
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              COSLP(IJ,K,IGST) = COSWDIF(IJ,K)
            ENDDO
          ENDDO
        ENDDO
      ENDIF


!*    2. MAIN LOOP OVER FREQUENCIES.
!        ---------------------------

      IF ( .NOT. LLNORMAGAM) THEN
        GAMNORMA(KIJS:KIJL,:) = 1.0_JWRB
      ENDIF

      IF ( .NOT. LLSNEG) THEN
        DSTAB(KIJS:KIJL,:,:) = 0.0_JWRB
      ENDIF

      DO M=1,NFRE

        IF (LTAUWSHELTER) THEN
          DO IGST=1,NGST
            DO IJ=KIJS,KIJL
              TAUPX=TAUX(IJ,IGST)-ABS_TAUWSHELTER*XSTRESS(IJ,IGST)
              TAUPY=TAUY(IJ,IGST)-ABS_TAUWSHELTER*YSTRESS(IJ,IGST)
              USDIRP(IJ,IGST)=ATAN2(TAUPX,TAUPY)
              USTP(IJ,IGST)=(TAUPX**2+TAUPY**2)**0.25_JWRB
              USTPM1(IJ,IGST)=1.0_JWRB/MAX(USTP(IJ,IGST),EPSUS)
            ENDDO
          ENDDO

          DO IGST=1,NGST
            DO K=1,NANG
              DO IJ=KIJS,KIJL
                COSLP(IJ,K,IGST) = COS(TH(K)-USDIRP(IJ,IGST))
              ENDDO
            ENDDO
          ENDDO

          DO IJ=KIJS,KIJL
            CONSTF(IJ) = ROGOROAIR(IJ)*CINV(IJ,M)*DFIM(M)
          ENDDO
        ENDIF


!*      PRECALCULATE FREQUENCY DEPENDENCE.
!       ----------------------------------

        DO IGST=1,NGST
          DO IJ=KIJS,KIJL
            UCN(IJ,IGST) = USTP(IJ,IGST)*CINV(IJ,M)
            UCNZALPD(IJ,IGST) = XKAPPA/(UCN(IJ,IGST) + ZALP)
          ENDDO
        ENDDO
        DO IJ=KIJS,KIJL
          ZCN(IJ) = LOG(WAVNUM(IJ,M)*Z0M(IJ))
          CNSN(IJ) = CONST(M)*RAORW(IJ)
        ENDDO

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

        DO K=1,NANG
          DO IJ=KIJS,KIJL
            XLLWS(IJ,K,M)=0.0_JWRB
          ENDDO
        ENDDO

        DO IGST=1,NGST
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              IF (COSLP(IJ,K,IGST) > 0.01_JWRB) THEN
                X    = COSLP(IJ,K,IGST)*UCN(IJ,IGST)
                ZLOG = ZCN(IJ) + UCNZALPD(IJ,IGST)/COSLP(IJ,K,IGST)
                IF (ZLOG < 0.0_JWRB) THEN
                  ZLOG2X=ZLOG*ZLOG*X
                  GAM0(IJ,K,IGST) = EXP(ZLOG)*ZLOG2X*ZLOG2X * CNSN(IJ)
                  XLLWS(IJ,K,M) = 1.0_JWRB
                ELSE
                  GAM0(IJ,K,IGST) = 0.0_JWRB
                ENDIF
              ELSE
                GAM0(IJ,K,IGST) = 0.0_JWRB
              ENDIF
            ENDDO
          ENDDO
        ENDDO


        IF (LLNORMAGAM) THEN

          DO IGST=1,NGST

            SUMF(KIJS:KIJL) = 0.0_JWRB
            SUMFSIN2(KIJS:KIJL) = 0.0_JWRB
            DO K=1,NANG
              DO IJ=KIJS,KIJL
                SUMF(IJ) = SUMF(IJ) + GAM0(IJ,K,IGST)*FL1(IJ,K,M)
                SUMFSIN2(IJ) = SUMFSIN2(IJ) + GAM0(IJ,K,IGST)*FL1(IJ,K,M)*SINWDIF2(IJ,K)
              ENDDO
            ENDDO

            DO IJ=KIJS,KIJL
              ZNZ = XNGAMCONST(IJ,M)*USTPM1(IJ,IGST)
              GAMNORMA(IJ,IGST) = (1.0_JWRB + ZNZ*SUMFSIN2(IJ)) / (1.0_JWRB + ZNZ*SUMF(IJ))
            ENDDO

          ENDDO

        ENDIF


        IF (LLSNEG) THEN
!       SWELL DAMPING:
          DO IJ=KIJS,KIJL
            DSTAB1(IJ) = COEF5(M)*AIRD_PVISC(IJ)*WAVNUM(IJ,M)
            TEMP1(IJ) = COEF(M)*RAORW(IJ)
          ENDDO

          DO IGST=1,NGST
            DO K=1,NANG
              DO IJ=KIJS,KIJL
                DSTAB2 = TEMP1(IJ)*(TEMP2(IJ)+(FU+FUD*COSLP(IJ,K,IGST))*USTP(IJ,IGST))
                DSTAB(IJ,K,IGST) = DSTAB1(IJ)+PTURB(IJ)*DSTAB2
              ENDDO
            ENDDO
          ENDDO
        ENDIF


!*    2.2 UPDATE THE SHELTERING STRESS (in any),
!         AND THEN ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         ---------------------------------------------------------

        DO K=1,NANG

          DO IGST=1,NGST
            DO IJ=KIJS,KIJL
              ! SLP: only the positive contributions
              SLP(IJ,IGST) =  GAM0(IJ,K,IGST) * GAMNORMA(IJ,IGST)
              FLP(IJ,IGST) = SLP(IJ,IGST)+DSTAB(IJ,K,IGST)
            ENDDO
          ENDDO

          DO IGST=1,NGST
            DO IJ=KIJS,KIJL
              SLP(IJ,IGST) = SLP(IJ,IGST)*FL1(IJ,K,M)
            ENDDO
          ENDDO

          IF (LTAUWSHELTER) THEN
            DO IJ=KIJS,KIJL
              CONST11(IJ)=CONSTF(IJ)*SINTH(K)
              CONST22(IJ)=CONSTF(IJ)*COSTH(K)
            ENDDO
            DO IGST=1,NGST
              DO IJ=KIJS,KIJL
                XSTRESS(IJ,IGST)=XSTRESS(IJ,IGST)+SLP(IJ,IGST)*CONST11(IJ)
                YSTRESS(IJ,IGST)=YSTRESS(IJ,IGST)+SLP(IJ,IGST)*CONST22(IJ)
              ENDDO
            ENDDO
          ENDIF

          IGST=1
            DO IJ=KIJS,KIJL
              SLP_AVG(IJ) = SLP(IJ,IGST)
              FLP_AVG(IJ) = FLP(IJ,IGST)
            ENDDO
          DO IGST=2,NGST
            DO IJ=KIJS,KIJL
              SLP_AVG(IJ) = SLP_AVG(IJ)+SLP(IJ,IGST)
              FLP_AVG(IJ) = FLP_AVG(IJ)+FLP(IJ,IGST)
            ENDDO
          ENDDO

          DO IJ=KIJS,KIJL
            SPOS(IJ,K,M) = AVG_GST*SLP_AVG(IJ)
          ENDDO

          DO IJ=KIJS,KIJL
            FLD(IJ,K,M) = AVG_GST*FLP_AVG(IJ)
            SL(IJ,K,M) = FLD(IJ,K,M)*FL1(IJ,K,M)
          ENDDO

        ENDDO

      ENDDO ! END LOOP OVER FREQUENCIES

      IF (LHOOK) CALL DR_HOOK('SINPUT_ARD',1,ZHOOK_HANDLE)

      END SUBROUTINE SINPUT_ARD
      SUBROUTINE SINPUT_JAN (NGST, LLSNEG, KIJS, KIJL, FL1 , &
     &                       WAVNUM, CINV, XK2CG,            &
     &                       WSWAVE, UFRIC, Z0M,     &
     &                       COSWDIF, SINWDIF2,              &
     &                       RAORW, WSTAR, RNFAC,            &
                             FLD, SL, SPOS, XLLWS)
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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LLNORMAGAM
      USE YOWFRED  , ONLY : ZPIFR    ,DELTH    ,TH
      USE YOWFRED  , ONLY : FR       ,TH       , ZPIFR
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,GM1      ,ZPI  , EPSUS
      USE YOWPHYS  , ONLY : ZALP     ,XKAPPA, BETAMAXOXKAPPA2
      USE YOWSTAT  , ONLY : IDAMPING
      USE YOWTEST  , ONLY : IU06

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: NGST
      LOGICAL, INTENT(IN) :: LLSNEG
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM, CINV, XK2CG
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: WSWAVE, UFRIC, Z0M
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG), INTENT(IN) :: COSWDIF, SINWDIF2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: RAORW, WSTAR, RNFAC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: FLD, SL, SPOS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      INTEGER(KIND=JWIM) :: IJ, IG, K, M
      INTEGER(KIND=JWIM) :: IGST

      REAL(KIND=JWRB) :: CONST1, CONST3, XKAPPAD
      REAL(KIND=JWRB) :: CONSTN
      REAL(KIND=JWRB) :: ZNZ
      REAL(KIND=JWRB) :: X, ZLOG, ZLOG2X, ZBETA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NGST) :: WSIN
      REAL(KIND=JWRB), DIMENSION(NFRE) :: CONST
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ZTANHKD 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: SIG_N
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: CNSN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: SUMF, SUMFSIN2 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: CSTRNFAC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: XNGAMCONST
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NGST) :: GAMNORMA ! ! RENORMALISATION FACTOR OF THE GROWTH RATE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NGST) :: SIGDEV ,US, Z0, UCN, ZCN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NGST) :: USTPM1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NGST) :: XVD, UCND, CONST3_UCN2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: UFAC1, UFAC2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: TEMPD
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NGST) :: GAM0

      LOGICAL, DIMENSION(KIJS:KIJL,NANG) :: LZ

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SINPUT_JAN',0,ZHOOK_HANDLE)

      CONST1   = BETAMAXOXKAPPA2 
      CONST3   = 2.0_JWRB*XKAPPA/CONST1  ! SEE IDAMPING
      XKAPPAD  = 1.E0_JWRB/XKAPPA

      CONST3 = IDAMPING*CONST3

      CONSTN = DELTH/(XKAPPA*ZPI)

!*    1. PRECALCULATED ANGULAR DEPENDENCE.
!        ---------------------------------

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          IF (COSWDIF(IJ,K) > 0.01_JWRB) THEN
            LZ(IJ,K) = .TRUE.
            TEMPD(IJ,K) = XKAPPA/COSWDIF(IJ,K)
          ELSE
            LZ(IJ,K) = .FALSE.
            TEMPD(IJ,K) = XKAPPA 
          ENDIF
        ENDDO
      ENDDO

      IF (LLNORMAGAM) THEN
        DO IJ=KIJS,KIJL
          CSTRNFAC(IJ) = CONSTN * RNFAC(IJ) / RAORW(IJ)
        ENDDO
        DO M=1,NFRE
          DO IJ=KIJS,KIJL
            XNGAMCONST(IJ,M) = CSTRNFAC(IJ)*XK2CG(IJ,M)
          ENDDO
        ENDDO

      ELSE
        XNGAMCONST(KIJS:KIJL,:) = 0.0_JWRB
      ENDIF

!     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
      DO IJ=KIJS,KIJL
        IF (NGST > 1) CALL WSIGSTAR (WSWAVE(IJ), UFRIC(IJ), Z0M(IJ), WSTAR(IJ), SIG_N(IJ))
      ENDDO


!     DEFINE WHERE SINPUT WILL BE EVALUATED IN RELATIVE TERM WRT USTAR
!     DEFINE ALSO THE RELATIVE WEIGHT OF EACH.

      IF (NGST == 1) THEN
        WSIN(1) = 1.0_JWRB 
        DO IJ=KIJS,KIJL
          SIGDEV(IJ,1) = 1.0_JWRB
        ENDDO
      ELSE IF (NGST == 2) THEN
        WSIN(1) = 0.5_JWRB 
        WSIN(2) = 0.5_JWRB 
        DO IJ=KIJS,KIJL
          SIGDEV(IJ,1) = 1.0_JWRB-SIG_N(IJ)
          SIGDEV(IJ,2) = 1.0_JWRB+SIG_N(IJ)
        ENDDO
      ELSE
         WRITE (IU06,*) '**************************************'
         WRITE (IU06,*) '*    FATAL ERROR                     *'
         WRITE (IU06,*) '*    ===========                     *'
         WRITE (IU06,*) '* IN SINPUT_JAN: NGST > 2            *'
         WRITE (IU06,*) '* NGST = ', NGST
         WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.  *'
         WRITE (IU06,*) '*                                    *'
         WRITE (IU06,*) '**************************************'
         CALL ABORT1
      ENDIF


      IF (NGST == 1) THEN
        DO IJ=KIJS,KIJL
          US(IJ,1) = UFRIC(IJ)
          Z0(IJ,1) = Z0M(IJ)
        ENDDO
      ELSE
        DO IGST=1,NGST
          DO IJ=KIJS,KIJL
            US(IJ,IGST) = UFRIC(IJ)*SIGDEV(IJ,IGST)
            Z0(IJ,IGST) = Z0M(IJ)
          ENDDO
        ENDDO
      ENDIF

      DO IGST=1,NGST
        DO IJ=KIJS,KIJL
          USTPM1(IJ,IGST) = 1.0_JWRB/MAX(US(IJ,IGST),EPSUS)
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    2. LOOP OVER FREQUENCIES.
!        ----------------------

      IF ( .NOT. LLNORMAGAM) THEN
        GAMNORMA(KIJS:KIJL,:) = 1.0_JWRB
      ENDIF

      IF ( .NOT. LLSNEG) THEN
        UFAC2(KIJS:KIJL,:) = 0.0_JWRB
      ENDIF

      DO M=1,NFRE

        CONST(M)=ZPIFR(M)*CONST1

!*      PRECALCULATE FREQUENCY DEPENDENCE.
!       ----------------------------------

        DO IJ=KIJS,KIJL
          ZTANHKD(IJ) = ZPIFR(M)**2/(G*WAVNUM(IJ,M)) 
        ENDDO

        DO IJ=KIJS,KIJL
          CNSN(IJ) = CONST(M)*ZTANHKD(IJ)*RAORW(IJ)
        ENDDO

        DO IGST=1,NGST
          DO IJ=KIJS,KIJL
            UCN(IJ,IGST) = US(IJ,IGST)*CINV(IJ,M) + ZALP
            CONST3_UCN2(IJ,IGST) = CONST3*UCN(IJ,IGST)**2
            UCND(IJ,IGST) = 1.0_JWRB/ UCN(IJ,IGST)
            ZCN(IJ,IGST)  = LOG(WAVNUM(IJ,M)*Z0(IJ,IGST))
            XVD(IJ,IGST) = 1.0_JWRB/(-US(IJ,IGST)*XKAPPAD*ZCN(IJ,IGST)*CINV(IJ,M))
          ENDDO
        ENDDO

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------


!       WIND INPUT:
        DO K=1,NANG

          DO IJ=KIJS,KIJL
            XLLWS(IJ,K,M)= 0.0_JWRB
          ENDDO

          DO IGST=1,NGST
            DO IJ=KIJS,KIJL
              IF (LZ(IJ,K)) THEN
                ZLOG = ZCN(IJ,IGST) + TEMPD(IJ,K)*UCND(IJ,IGST)
                IF (ZLOG < 0.0_JWRB) THEN
                  X=COSWDIF(IJ,K)*UCN(IJ,IGST)
                  ZLOG2X=ZLOG*ZLOG*X
                  GAM0(IJ,K,IGST) = ZLOG2X*ZLOG2X*EXP(ZLOG) * CNSN(IJ)
                  XLLWS(IJ,K,M)= 1.0_JWRB
                ELSE
                  GAM0(IJ,K,IGST) = 0.0_JWRB
                ENDIF
              ELSE
                GAM0(IJ,K,IGST) = 0.0_JWRB
              ENDIF
            ENDDO
          ENDDO

        ENDDO


        IF (LLNORMAGAM) THEN

          DO IGST=1,NGST

            SUMF(KIJS:KIJL) = 0.0_JWRB
            SUMFSIN2(KIJS:KIJL) = 0.0_JWRB
            DO K=1,NANG
              DO IJ=KIJS,KIJL
                SUMF(IJ) = SUMF(IJ) + GAM0(IJ,K,IGST)*FL1(IJ,K,M)
                SUMFSIN2(IJ) = SUMFSIN2(IJ) + GAM0(IJ,K,IGST)*FL1(IJ,K,M)*SINWDIF2(IJ,K)
              ENDDO
            ENDDO

            DO IJ=KIJS,KIJL
              ZNZ = XNGAMCONST(IJ,M)*USTPM1(IJ,IGST)
              GAMNORMA(IJ,IGST) = (1.0_JWRB + ZNZ*SUMFSIN2(IJ)) / (1.0_JWRB + ZNZ*SUMF(IJ))
            ENDDO

          ENDDO

        ENDIF


        DO K=1,NANG
          DO IJ=KIJS,KIJL
            UFAC1(IJ,K) = WSIN(1)*GAM0(IJ,K,1)*GAMNORMA(IJ,1)
          ENDDO
          DO IGST=2,NGST
            DO IJ=KIJS,KIJL
              UFAC1(IJ,K) = UFAC1(IJ,K) + WSIN(IGST)*GAM0(IJ,K,IGST)*GAMNORMA(IJ,IGST)
            ENDDO
          ENDDO
        ENDDO

        IF (LLSNEG) THEN
!         SWELL DAMPING:
          DO K=1,NANG
            DO IGST=1,1
              DO IJ=KIJS,KIJL
                ZBETA = CONST3_UCN2(IJ,IGST)*(COSWDIF(IJ,K)-XVD(IJ,IGST))
                UFAC2(IJ,K) = WSIN(IGST)*ZBETA
              ENDDO
            ENDDO
            DO IGST=2,NGST
              DO IJ=KIJS,KIJL
                ZBETA = CONST3_UCN2(IJ,IGST)*(COSWDIF(IJ,K)-XVD(IJ,IGST))
                UFAC2(IJ,K) = UFAC2(IJ,K)+WSIN(IGST)*ZBETA
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         ------------------------------------------------

        DO K=1,NANG
          DO IJ=KIJS,KIJL
            FLD(IJ,K,M) = UFAC1(IJ,K) + UFAC2(IJ,K)*CNSN(IJ)
            SPOS(IJ,K,M) = UFAC1(IJ,K)*FL1(IJ,K,M)
            SL(IJ,K,M) = FLD(IJ,K,M)*FL1(IJ,K,M)
          ENDDO
        ENDDO

      ENDDO

      IF (LHOOK) CALL DR_HOOK('SINPUT_JAN',1,ZHOOK_HANDLE)

      END SUBROUTINE SINPUT_JAN
END MODULE SINPUT_ARD_MOD
