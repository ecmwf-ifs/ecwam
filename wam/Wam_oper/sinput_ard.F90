      SUBROUTINE SINPUT_ARD (F,FL,IJS,IJL,THWNEW,USNEW,Z0NEW,&
     &                       ROAIRN,WSTAR,SL,XLLWS)
! ----------------------------------------------------------------------

!**** *SINPUT* - COMPUTATION OF INPUT SOURCE FUNCTION.

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


!       L. AOUF    March 2011 : USE OF NEW DISSIPATION DEVELOPED BY ARDHUIN ET AL.2010

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

!     *CALL* *SINPUT (F, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,
!    &                   ROAIRN,WSTAR, SL, XLLWS)
!            *F* - SPECTRUM.
!           *FL* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!          *IJS* - INDEX OF FIRST GRIDPOINT.
!          *IJL* - INDEX OF LAST GRIDPOINT.
!       *THWNEW* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!        *USNEW* - NEW FRICTION VELOCITY IN M/S.
!        *Z0NEW* - ROUGHNESS LENGTH IN M.
!       *ROAIRN* - AIR DENSITY IN KG/M3
!        *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
!           *SL* - TOTAL SOURCE FUNCTION ARRAY.
!         *XLLWS*- 1 WHERE SINPUT IS POSITIVE

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

      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,TAUWSHELTER, XKAPPA
      USE YOWFRED  , ONLY : FR       ,TH       ,DFIM     ,COSTH  ,SINTH
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER
      USE YOWSHAL  , ONLY : TFAK     ,INDEP, DEPTH
      USE YOWSTAT  , ONLY : ISHALLO
      USE YOWTABL  , ONLY : IAB      ,SWELLFT
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS,IJL
      INTEGER(KIND=JWIM) :: IJ,K,M,IND,IGST

      REAL(KIND=JWRB), PARAMETER :: ABMIN = 0.3_JWRB
      REAL(KIND=JWRB), PARAMETER :: ABMAX = 8.0_JWRB 

      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: THWNEW, USNEW, Z0NEW, ROAIRN, WSTAR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL,SL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(OUT) :: XLLWS
      REAL(KIND=JWRB) :: CONST1
      REAL(KIND=JWRB) :: X1,X2,ZLOG,ZLOG1,ZLOG2,ZLOG2X,XV1,XV2,ZBETA1,ZBETA2
      REAL(KIND=JWRB) :: XI,X,DELI1,DELI2
      REAL(KIND=JWRB) :: FU,FUD,NU_AIR,SWELLFPAR,SWELLF,SWELLF2,SWELLF3,SWELLF4,SWELLF5
      REAL(KIND=JWRB) :: SWELLF7, SMOOTH
      REAL(KIND=JWRB) :: ARG, DELAB, CONST11, CONST22
      REAL(KIND=JWRB) :: TAUPX,TAUPY
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NANG) :: PP
      REAL(KIND=JWRB), DIMENSION(NFRE) :: FAC, CONST, SIG, CONSTF, COEF, COEF5
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: Z0VIS, Z0NOZ, FWW
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: PVISC, PTURB
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: ZCN, CM
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: SIG_N, UORBT, AORB, TEMP, RE, ZORB
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: CNSN
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,2) :: XSTRESS, YSTRESS, FLP, SLP
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,2) :: USG2, TAUX, TAUY, USTP, USDIRP, UCN
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE) :: XK
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG) :: DSTAB1
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,2) :: TEMP1, UFAC, DSTAB2, DSTAB

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SINPUT_ARD',0,ZHOOK_HANDLE)
#endif

      CONST1  = BETAMAX/XKAPPA**2 /ROWATER
      NU_AIR = 1.4E-05_JWRB
      SWELLFPAR = 3.0_JWRB
      SWELLF = 0.8_JWRB
      SWELLF2 = -0.018_JWRB
      SWELLF3 = 0.015_JWRB
      SWELLF4 = 1.0E05_JWRB
      SWELLF5 = 1.2_JWRB
      SWELLF7 = 2.3E05_JWRB

      FU=ABS(SWELLF3)
      FUD=SWELLF2
      DELAB= (ABMAX-ABMIN)/REAL(IAB)


!     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
      CALL WSIGSTAR (IJS, IJL, USNEW, Z0NEW, WSTAR, SIG_N)

! ----------------------------------------------------------------------
! computation of Uorb and Aorb

      DO IJ=IJS,IJL
        UORBT(IJ) = 0._JWRB
        AORB(IJ) = 0._JWRB
      ENDDO
      DO M=1,NFRE
        K=1
        SIG(M) = ZPI*FR(M)
        DO IJ=IJS,IJL
          TEMP(IJ) = F(IJ,K,M)
        ENDDO
        DO K=2,NANG
          DO IJ=IJS,IJL
            TEMP(IJ) = TEMP(IJ)+F(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          UORBT(IJ) = UORBT(IJ)+DFIM(M)*(SIG(M)**2)*TEMP(IJ)
          AORB(IJ) = AORB(IJ)+DFIM(M)*TEMP(IJ)
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        UORBT(IJ) = 2*SQRT(UORBT(IJ))  ! this is the significant orbital amplitude
        AORB(IJ) = 2*SQRT(AORB(IJ))
        Z0VIS(IJ) = 0.1_JWRB*NU_AIR/MAX(USNEW(IJ),0.0001_JWRB)
        Z0NOZ(IJ) = max(Z0VIS(IJ),0.04_JWRB*Z0NEW(IJ))
        ZORB(IJ) = AORB(IJ)/Z0NOZ(IJ)
        IF (UORBT(IJ).NE.0.0_JWRB) THEN
        RE(IJ) = 4*UORBT(IJ)*AORB(IJ) / NU_AIR ! this is the Reynolds number 
        ELSE
        RE(IJ) = 0.0_JWRB
        ENDIF
! calcul fww
        FU=ABS(SWELLF3)
        FUD=SWELLF2
        XI=(LOG10(MAX(ZORB(IJ),3.0_JWRB))-ABMIN)/DELAB
        IND  = MIN (IAB-1, INT(XI))
        DELI1= MIN (1.0_JWRB ,XI-FLOAT(IND))
        DELI2= 1.0_JWRB - DELI1
        FWW(IJ) =SWELLFT(IND)*DELI2+SWELLFT(IND+1)*DELI1
!!!!!!!!!!!!!!!!

      END DO

      DO IGST=1,2
        DO IJ=IJS,IJL
          XSTRESS(IJ,IGST)=0.0_JWRB
          YSTRESS(IJ,IGST)=0.0_JWRB
        ENDDO
      ENDDO
      DO IJ=IJS,IJL
        USG2(IJ,1)=(USNEW(IJ)*(1.0_JWRB+SIG_N(IJ)))**2
        USG2(IJ,2)=(USNEW(IJ)*(1.0_JWRB-SIG_N(IJ)))**2
      ENDDO

      DO IGST=1,2
        DO IJ=IJS,IJL
          TAUX(IJ,IGST)=USG2(IJ,IGST)*SIN(THWNEW(IJ))
          TAUY(IJ,IGST)=USG2(IJ,IGST)*COS(THWNEW(IJ))
        ENDDO
      ENDDO

!*    2. LOOP OVER FREQUENCIES.
!        ----------------------

      DO M=1,NFRE

        DO IGST=1,2
          DO IJ=IJS,IJL
            TAUPX=TAUX(IJ,IGST)-ABS(TAUWSHELTER)*XSTRESS(IJ,IGST)
            TAUPY=TAUY(IJ,IGST)-ABS(TAUWSHELTER)*YSTRESS(IJ,IGST)
            USTP(IJ,IGST)=(TAUPX**2+TAUPY**2)**0.25_JWRB
            USDIRP(IJ,IGST)=ATAN2(TAUPX,TAUPY)
          ENDDO
        ENDDO

        CONSTF(M) =ZPI*ROWATER*FR(M)*DFIM(M)
        FAC(M) = ZPI*FR(M)
        CONST(M)=FAC(M)*CONST1
        COEF(M) =-SWELLF*16._JWRB*SIG(M)**2/(G*ROWATER)
        COEF5(M)=-SWELLF5*2._JWRB*SQRT(2._JWRB*NU_AIR*SIG(M))/ROWATER

        DO IGST=1,2
          DO K=1,NANG
            DO IJ=IJS,IJL
              TEMP1(IJ,K,IGST) = COS(TH(K)-USDIRP(IJ,IGST))
            ENDDO
          ENDDO
        ENDDO

!*      INVERSE OF PHASE VELOCITIES AND WAVE NUMBER.
!       -------------------------------------------

        IF (ISHALLO.EQ.1) THEN
          DO IJ=IJS,IJL
            CM(IJ) = FAC(M)/G
            XK(IJ,M) = ((ZPI*FR(M))**2)/G
          ENDDO
        ELSE
          DO IJ=IJS,IJL
            CM(IJ) = TFAK(INDEP(IJ),M)/FAC(M)
            XK(IJ,M) = TFAK(INDEP(IJ),M)
          ENDDO
        ENDIF

!*      PRECALCULATE FREQUENCY DEPENDENCE.
!       ----------------------------------

        DO IGST=1,2
          DO IJ=IJS,IJL
            UCN(IJ,IGST) = USTP(IJ,IGST)*CM(IJ) + ZALP
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          ZCN(IJ) = LOG(G*Z0NEW(IJ)*CM(IJ)**2)
          CNSN(IJ) = CONST(M) * ROAIRN(IJ)
        ENDDO

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

        DO K=1,NANG
          DO IJ=IJS,IJL
            XLLWS(IJ,K,M)= 0.0_JWRB
          ENDDO
        ENDDO

        DO IGST=1,2
          DO K=1,NANG
            DO IJ=IJS,IJL
              IF (TEMP1(IJ,K,IGST).GT.0.01_JWRB) THEN
                X    = TEMP1(IJ,K,IGST)*UCN(IJ,IGST)
                ZLOG = ZCN(IJ) + XKAPPA/X
                IF (ZLOG.LT.0.0_JWRB) THEN
                  ZLOG2X=ZLOG*ZLOG*X
                  UFAC(IJ,K,IGST) = EXP(ZLOG)*ZLOG2X*ZLOG2X
                  XLLWS(IJ,K,M)= 1.0_JWRB
                ELSE
                  UFAC(IJ,K,IGST) = 0.0_JWRB
                ENDIF
              ELSE
                UFAC(IJ,K,IGST) = 0.0_JWRB
              ENDIF
            ENDDO
          ENDDO
        ENDDO

!       SWELL DAMPING:
        DO K=1,NANG
          PP(K) = 1.
          DO IJ=IJS,IJL
            DSTAB1(IJ,K) = COEF5(M)*ROAIRN(IJ)*XK(IJ,M)*PP(K)
          END DO

          DO IGST=1,2
            DO IJ=IJS,IJL
              DSTAB2(IJ,K,IGST) = COEF(M)*ROAIRN(IJ)*(FWW(IJ)*UORBT(IJ)+(FU+FUD*TEMP1(IJ,K,IGST))*USTP(IJ,IGST))
            END DO
          ENDDO

        END DO

        IF (SWELLF7.GT.0.0_JWRB) THEN
          DO IJ=IJS,IJL
            SMOOTH=0.5_JWRB*TANH((RE(IJ)-SWELLF4)/SWELLF7)
            PTURB(IJ)=0.5_JWRB+SMOOTH
            PVISC(IJ)=0.5_JWRB-SMOOTH
          ENDDO
          DO IGST=1,2
            DO K=1,NANG
              DO IJ=IJS,IJL
                DSTAB(IJ,K,IGST) = PVISC(IJ)*DSTAB1(IJ,K)+PTURB(IJ)*DSTAB2(IJ,K,IGST)
              ENDDO
            END DO
          ENDDO
        ELSE
          DO IJ=IJS,IJL
            PTURB(IJ)=0.5_JWRB
            PVISC(IJ)=0.5_JWRB
          ENDDO
          IF (RE(IJ).LE.SWELLF4) THEN
            DO IGST=1,2
              DO K=1,NANG
                DO IJ=IJS,IJL
                  DSTAB(IJ,K,IGST) = PVISC(IJ)*DSTAB1(IJ,K)
                ENDDO
              END DO
            ENDDO
          ELSE
            DO IGST=1,2
              DO K=1,NANG
                DO IJ=IJS,IJL
                DSTAB(IJ,K,IGST) = PTURB(IJ)*DSTAB2(IJ,K,IGST)
                ENDDO
              END DO
            ENDDO
          ENDIF
        ENDIF

!*    2.2 UPDATE THE SHELTERING STRESS,
!         AND THEN ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         ---------------------------------------------------------

        DO K=1,NANG
          CONST11=CONSTF(M)*SINTH(K)
          CONST22=CONSTF(M)*COSTH(K)
          DO IGST=1,2
            DO IJ=IJS,IJL
              FLP(IJ,IGST) = CNSN(IJ)*UFAC(IJ,K,IGST)+DSTAB(IJ,K,IGST)
              SLP(IJ,IGST) = FLP(IJ,IGST)*F(IJ,K,M)
              XSTRESS(IJ,IGST)=XSTRESS(IJ,IGST)+SLP(IJ,IGST)*CONST11/MAX(ROAIRN(IJ),1.0_JWRB)
              YSTRESS(IJ,IGST)=YSTRESS(IJ,IGST)+SLP(IJ,IGST)*CONST22/MAX(ROAIRN(IJ),1.0_JWRB)
            ENDDO
          ENDDO
          DO IJ=IJS,IJL
            FL(IJ,K,M) = 0.5_JWRB*(FLP(IJ,1)+FLP(IJ,2))
            SL(IJ,K,M) = FL(IJ,K,M)*F(IJ,K,M)
          ENDDO
        ENDDO


      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SINPUT_ARD',1,ZHOOK_HANDLE)
#endif

      RETURN
      END SUBROUTINE SINPUT_ARD
